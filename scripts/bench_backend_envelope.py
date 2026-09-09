#!/usr/bin/env python3
"""Backend and device envelope for the PINNACLE CP-PFD benchmark (research question R2).

Measures where GPU/Torch execution wins and where it loses, separately for the
two workload shapes the coupled problem contains:

  cp      the dislocation-density CP constitutive update over a quadrature
          batch -- many small 3x3 einsums, the shape ADR 0015 warns about;
  stack   the state marshalling inside that update (list-of-per-QP-dicts ->
          stacked tensors), isolated so dispatch overhead can be separated
          from arithmetic;
  pinn    a float64 PINN training step with double-backward (create_graph),
          the shape the damage-as-PINN-field route trains on.

Protocol follows the platform-benchmark discipline: capability check before
speed, one configuration per process, time-budgeted warmup, median of N
synchronised reps with min-max spread reported, bounded loops, threads pinned
in-process. A configuration whose first warmup iteration exceeds --cap-s is
recorded as capped rather than run to completion, so a pathological
configuration cannot hang the job.

Usage (one configuration per invocation):
    python bench_backend_envelope.py --bench cp --backend torch --device cuda \
        --n 4096 --out results.jsonl
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import time

WARMUP_BUDGET_S = 3.0
MIN_WARMUP_ITERS = 8
N_REP = 11

AL2219_MATERIAL = {
    # Independently extracted Al2219 set, micrometre-MPa-second units.
    # Provenance and the (k0, y_c) -> (k1, k2) derivation are Phase 1 work;
    # these are used here only to make the timing representative of the real
    # problem's stiffness, not to produce a physics result.
    "type": "dislocation_density_cp",
    "C11": 108.0e3,
    "C12": 61.3e3,
    "C44": 28.5e3,
    "mu": 25.7e3,
    "gamma0": 1.0e-5,
    "rate_exponent_n": 1.0 / 0.12,
    "tau_c0": 0.112,
    "q_lat": 1.4,
    "burgers": 0.000286,
    "rho0": 220.0,
    "k1": 0.3432,
    "k2": 2.0 * 0.00268,
}
AL2219_CRYSTAL = {"structure": "fcc", "euler_angles": [0.0, 0.0, 0.0]}


def timed(fn, sync=None):
    """Median / min / max wall time of N_REP individually synchronised reps."""
    t0 = time.perf_counter()
    iters = 0
    while iters < MIN_WARMUP_ITERS or (time.perf_counter() - t0) < WARMUP_BUDGET_S:
        fn()
        iters += 1
        if iters > 2000:
            break
    if sync:
        sync()
    samples = []
    for _ in range(N_REP):
        t = time.perf_counter()
        fn()
        if sync:
            sync()
        samples.append(time.perf_counter() - t)
    med = statistics.median(samples)
    spread = (max(samples) - min(samples)) / med * 100.0 if med > 0 else float("nan")
    return {
        "median_s": med,
        "min_s": min(samples),
        "max_s": max(samples),
        "spread_pct": spread,
        "warmup_iters": iters,
    }


def make_backend(backend_name, device):
    from pinnacle.backend import get_backend

    return get_backend(backend_name, precision="float64", device=device, seed=0)


def capability_check(backend_name, device):
    """Allocate a tiny tensor on the target device before measuring anything."""
    if backend_name != "torch" or device != "cuda":
        return True, "n/a"
    import torch

    if not torch.cuda.is_available():
        return False, "torch.cuda.is_available() is False"
    x = torch.ones(4, device="cuda", dtype=torch.float64)
    torch.cuda.synchronize()
    if float(x.sum()) != 4.0:
        return False, "tiny allocation returned the wrong value"
    return True, torch.cuda.get_device_name(0)


def build_cp(backend, n_qp):
    from pinnacle import registry
    from pinnacle.backend import ops
    from pinnacle.plugins import discover_components

    discover_components()

    mat = registry.create(
        "material",
        "dislocation_density_cp",
        backend=backend,
        config={"material": AL2219_MATERIAL, "crystal": AL2219_CRYSTAL},
    )
    eye = ops.eye(backend, 3)
    state0 = mat.initial_state()
    state_batch = [{k: v for k, v in state0.items()} for _ in range(n_qp)]
    F_batch = backend.stack([eye for _ in range(n_qp)], axis=0)
    return mat, F_batch, state_batch


def bench_cp(backend, n_qp, sync):
    mat, F_batch, state_batch = build_cp(backend, n_qp)
    return lambda: mat.compute_batch(F_batch, state_batch, dt=1.0e-3), sync


def bench_stack(backend, n_qp, sync):
    """Isolate the state-marshalling step inside compute_batch."""
    _, _, state_batch = build_cp(backend, n_qp)

    def run():
        backend.stack([s["Fp"] for s in state_batch], axis=0)
        backend.stack([s["rho"] for s in state_batch], axis=0)

    return run, sync


def bench_pinn(backend, n_colloc, sync):
    """float64 MLP forward plus double-backward, the PINN residual shape."""
    import torch

    dev = "cuda" if backend.device.startswith("cuda") else "cpu"
    net = torch.nn.Sequential(
        torch.nn.Linear(3, 64),
        torch.nn.Tanh(),
        torch.nn.Linear(64, 64),
        torch.nn.Tanh(),
        torch.nn.Linear(64, 1),
    ).double().to(dev)
    opt = torch.optim.Adam(net.parameters(), lr=1e-3)
    x = torch.rand(n_colloc, 3, device=dev, dtype=torch.float64, requires_grad=True)

    def run():
        opt.zero_grad(set_to_none=True)
        d = net(x)
        g = torch.autograd.grad(d.sum(), x, create_graph=True)[0]
        lap = torch.autograd.grad(g.sum(), x, create_graph=True)[0].sum(dim=1)
        loss = ((lap + d.squeeze(-1)) ** 2).mean()
        loss.backward()
        opt.step()

    return run, sync


BENCHES = {"cp": bench_cp, "stack": bench_stack, "pinn": bench_pinn}


def run_one(args):
    ok, note = capability_check(args.backend, args.device)
    rec = {
        "bench": args.bench,
        "backend": args.backend,
        "device": args.device,
        "n": args.n,
        "precision": "float64",
        "threads": args.threads,
        "host": platform.node() or "masked",
        "capability_ok": ok,
        "capability_note": note,
    }
    if not ok:
        rec["status"] = "device_unavailable"
        return rec

    if args.backend == "torch":
        import torch

        torch.set_num_threads(args.threads)
        rec["torch_version"] = torch.__version__
        sync = torch.cuda.synchronize if args.device == "cuda" else None
    else:
        sync = None

    backend = make_backend(args.backend, args.device)
    rec["backend_device"] = backend.device

    if args.bench == "pinn" and args.backend != "torch":
        rec["status"] = "not_applicable"
        rec["reason"] = "the NumPy backend cannot train a PINN (no autodiff)"
        return rec

    build_t0 = time.perf_counter()
    fn, sync = BENCHES[args.bench](backend, args.n, sync)
    rec["build_s"] = time.perf_counter() - build_t0

    # Bounded single-shot probe before any loop.
    probe_t0 = time.perf_counter()
    fn()
    if sync:
        sync()
    probe = time.perf_counter() - probe_t0
    rec["probe_s"] = probe
    if probe > args.cap_s:
        rec["status"] = "capped"
        rec["reason"] = f"single iteration {probe:.1f}s exceeds cap {args.cap_s}s"
        return rec

    rec.update(timed(fn, sync))
    rec["status"] = "ok"
    if args.bench in ("cp", "stack"):
        rec["per_qp_us"] = rec["median_s"] / args.n * 1e6
    else:
        rec["per_point_us"] = rec["median_s"] / args.n * 1e6
    return rec


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bench", choices=sorted(BENCHES), required=True)
    p.add_argument("--backend", choices=["numpy", "torch"], required=True)
    p.add_argument("--device", choices=["cpu", "cuda"], default="cpu")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    p.add_argument("--cap-s", type=float, default=30.0)
    p.add_argument("--out", default="bench_backend_envelope.jsonl")
    args = p.parse_args()

    try:
        rec = run_one(args)
    except Exception as exc:  # a failed configuration must not lose the others
        rec = {
            "bench": args.bench,
            "backend": args.backend,
            "device": args.device,
            "n": args.n,
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }
    with open(args.out, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(json.dumps(rec))
    return 0


if __name__ == "__main__":
    sys.exit(main())
