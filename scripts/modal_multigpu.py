#!/usr/bin/env python3
"""Part V: scaling of the PINNACLE device pool over several NVIDIA GPUs, run on Modal.

usage: modal run scripts/modal_multigpu.py [--gpu A10G] [--counts 1,2,4] [--steps 200] [--mesh 5] [--no-mechanics] [--profile]

Run with the client of `.venv-modal/` (the author's token). The container asks for 8 CPU
cores: with Modal's default share of a core the host-side work of a step (services,
Python) hides the GPUs (first run of 2026-11-25). The script sends the PINNACLE source and the
decks of this repository to Modal and runs them there; nothing else leaves the machine.

What it measures, on one container with the largest count of `--counts` GPUs:

  explicit   the single track with melt-pool flow (decks/lpbf/in718/k1/cases/flow_h10/
             case0.i, mesh refined to `--mesh` um by an override), `--steps` steps, with
             `[Global] devices` set to the first 1, 2, 4, ... GPUs: seconds per step
             and the largest difference of the final temperature against one GPU
  mechanics  the tensile cube of case K3 (k3_cube_n6_seed1.i refined to 16 x 16 x 16
             elements, 32,768 quadrature points, 3 load steps) with the crystal-plasticity
             points on the same device lists: wall time

It prints one JSON document; copy it into results/lpbf/hardware/ and add a ledger row.
Cost: the container bills all its GPUs for the whole run; with the defaults (4 x A10G,
about 15 minutes) a few US dollars. Check the current prices of Modal before a larger run.
"""

from __future__ import annotations

import json
from pathlib import Path

import modal

REPO = Path(__file__).resolve().parent.parent
PINNACLE = REPO.parent / "pinnacle"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch>=2.4", "numpy>=1.26", "scipy>=1.11", "jsonschema>=4.0", "pyyaml", "h5py")
    .add_local_dir(str(PINNACLE / "src"), "/work/pinnacle/src", copy=True)
    .add_local_dir(str(PINNACLE / "schemas"), "/work/pinnacle/schemas", copy=True)
    .add_local_dir(str(REPO / "decks" / "lpbf"), "/work/bench/decks/lpbf", copy=True)
    .add_local_dir(str(REPO / "scripts"), "/work/bench/scripts", copy=True)
    .add_local_dir(str(REPO / "data" / "lpbf" / "in718"), "/work/bench/data/lpbf/in718", copy=True)
)
app = modal.App("pinnacle-device-pool-scaling", image=image)


def _run(script: str, deck: str, overrides: list[str], out: str, profile: bool = False) -> dict:
    import os
    import subprocess
    import time

    env = dict(os.environ, PYTHONPATH="/work/pinnacle/src")
    tic = time.time()
    proc = subprocess.run(
        ["python", *(["-m", "cProfile", "-o", f"{out}.prof"] if profile else []),
         f"/work/bench/scripts/{script}", *deck.split(), *overrides, f"global.output_dir={out}"],
        cwd="/work/bench", env=env, capture_output=True, text=True,
    )
    wall = time.time() - tic
    report = {"wall_with_setup_s": round(wall, 1), "returncode": proc.returncode}
    wall_file = Path(out) / "wall.json"
    if wall_file.exists():
        report.update(json.loads(wall_file.read_text()))
    else:
        report["stderr_tail"] = proc.stderr[-2000:]
    if profile and Path(f"{out}.prof").exists():
        import io
        import pstats

        text = io.StringIO()
        pstats.Stats(f"{out}.prof", stream=text).sort_stats("cumulative").print_stats(70)
        report["profile_main_thread"] = text.getvalue()[-12000:]
    return report


def _measure(counts: list[int], steps: int, mesh_um: float, mechanics: bool = True,
             profile: bool = False) -> dict:
    import numpy as np
    import torch

    names = [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]
    report: dict = {"gpus": names, "torch": str(torch.__version__), "explicit": {}, "mechanics": {}}
    n = [int(round(v / mesh_um)) for v in (700.0, 300.0, 350.0)]
    dt = 5.2869e-07 * (mesh_um / 10.0) ** 2
    track = "decks/lpbf/in718/k1/cases/flow_h10/case0.i"
    cube = "run decks/lpbf/in718/k3/k3_cube_n6_seed1.i"
    reference = None
    for count in counts:
        devices = json.dumps(" ".join(f"cuda:{i}" for i in range(count)))
        out = f"/tmp/explicit_{count}"
        row = _run(
            "lpbf_run.py", track,
            [f"mesh.divisions={json.dumps(n)}", f"executioner.dt={dt}",
             f"executioner.end_time={dt * steps}", "executioner.snapshots_every=0",
             f"kernels.carried.mesh_size={mesh_um * 1e-6}", f"kernels.momentum.mesh_size={mesh_um * 1e-6}",
             f"kernels.surface.band_thickness={mesh_um * 1e-6}",
             f"kernels.marangoni.band_thickness={mesh_um * 1e-6}",
             f"global.devices={devices}", "global.pool_min_batch=1", "global.num_threads=8"],
            out, profile,
        )
        fields = Path(out) / "fields.npz"
        if fields.exists():
            temperature = np.load(fields)["temperature"].ravel()
            if reference is None:
                reference = temperature
            row["max_abs_difference_to_one_gpu_K"] = float(np.abs(temperature - reference).max())
        report["explicit"][str(count)] = row
        if not mechanics:
            continue
        report["mechanics"][str(count)] = _run(
            "lpbf_k3.py", cube,
            ["global.use_gpu=true", f"global.devices={devices}", "global.pool_min_batch=1",
             "mesh.divisions=[16,16,16]", "executioner.n_steps=3", "executioner.t_end=6.0",
             "bcs.pull.values=[7.2e-7]", "global.num_threads=8"],
            f"/tmp/mechanics_{count}",
        )
    return report


@app.function(gpu="A10G:4", cpu=8.0, memory=32768, timeout=3600)
def measure_a10g(counts: list[int], steps: int, mesh_um: float, mechanics: bool = True,
                 profile: bool = False) -> dict:
    return _measure(counts, steps, mesh_um, mechanics, profile)


@app.function(gpu="L4:4", cpu=8.0, memory=32768, timeout=3600)
def measure_l4(counts: list[int], steps: int, mesh_um: float, mechanics: bool = True,
                 profile: bool = False) -> dict:
    return _measure(counts, steps, mesh_um, mechanics, profile)


@app.function(gpu="H100:8", cpu=8.0, memory=32768, timeout=1800)
def measure_h100(counts: list[int], steps: int, mesh_um: float, mechanics: bool = True,
                 profile: bool = False) -> dict:
    return _measure(counts, steps, mesh_um, mechanics, profile)


@app.function(gpu="A10G:4", cpu=8.0, memory=32768, timeout=900)
def diagnose() -> dict:
    """Where do the arrays of a pool live, and do four GPUs work at the same time at all?"""
    import sys
    import time
    from concurrent.futures import ThreadPoolExecutor

    import torch

    sys.path.insert(0, "/work/pinnacle/src")
    from pinnacle.backend import get_backend
    from pinnacle.backend.device_pool import DevicePool

    backend = get_backend("torch", precision="float64", device="cuda:0")
    devices = [f"cuda:{i}" for i in range(torch.cuda.device_count())]
    pool = DevicePool(backend, devices)
    report: dict = {
        "pool_devices": devices,
        "backend_tensor_devices": [str(b.tensor([0.0]).device) for b in pool.backends],
        "namespace_asarray_devices": [str(ns.asarray([0.0]).device) for ns in pool.namespaces],
        "namespace_index_devices": [str(ns.asarray([0], integer=True).device) for ns in pool.namespaces],
    }

    def work(device: str, n: int) -> float:
        x = torch.rand(n, 8, 3, dtype=torch.float64, device=device)
        total = torch.zeros((), dtype=torch.float64, device=device)
        for _ in range(300):
            total = total + (x * x).sum(dim=1).exp().sum()
        torch.cuda.synchronize(device)
        return float(total)

    n = 2_000_000
    for device in devices:
        work(device, 1000)
    tic = time.time(); work("cuda:0", n); one = time.time() - tic
    tic = time.time()
    for device in devices:
        work(device, n // len(devices))
    serial = time.time() - tic
    tic = time.time()
    with ThreadPoolExecutor(len(devices)) as ex:
        list(ex.map(lambda d: work(d, n // len(devices)), devices))
    threaded = time.time() - tic
    report["micro_benchmark_s"] = {"all_on_one_gpu": one, "quarters_one_after_another": serial,
                                   "quarters_in_threads": threaded}
    return report


@app.local_entrypoint()
def check() -> None:
    print(json.dumps(diagnose.remote(), indent=1))


@app.local_entrypoint()
def main(gpu: str = "A10G", counts: str = "1,2,4", steps: int = 200, mesh: float = 5.0,
         mechanics: bool = True, profile: bool = False) -> None:
    wanted = [int(v) for v in counts.split(",")]
    function = {"A10G": measure_a10g, "L4": measure_l4, "H100": measure_h100}[gpu]
    print(json.dumps(function.remote(wanted, steps, mesh, mechanics, profile), indent=1))
