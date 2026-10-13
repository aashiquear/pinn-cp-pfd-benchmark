#!/usr/bin/env python3
"""Can a smooth network reach the strong-form residual of the FE solution on the plate?

Diagnostic for the physics-trained predictor (PINNACLE ADR 0028): on the
polycrystal plate its scaled residual loss stays O(1). This script separates
"the optimiser does not find it" from "the target is not reachable":

1. run the plate deck on the torch backend for ``--steps`` load steps (FE,
   Newton), keeping the live quadrature-point state;
2. fit the predictor's network (same architecture, same plane-Dirichlet
   ansatz) to the FE nodal displacement of the last step by REGRESSION only;
3. evaluate the predictor's physics loss (scaled strong-form residual at the
   quadrature points with the frozen state) on that regressed field, and
   the FE residual / reaction ratio of the regressed field through the
   assembler.

If a network that reproduces the FE field to a small relative error still
carries an O(1) strong-form loss, the collocation target itself is the
obstacle (the equilibrium field of a 30-grain plate has strain kinks at the
grain boundaries that a global smooth network cannot represent), not the
training.

Usage: strong_form_floor.py results/cp_study/cal/gn2/base/cpfe.i --steps 4 --epochs 3000
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("deck")
    ap.add_argument("--steps", type=int, default=4)
    ap.add_argument("--epochs", type=int, default=3000)
    ap.add_argument("--neurons", type=int, default=64)
    ap.add_argument("--layers", type=int, default=3)
    ap.add_argument("--device", default="auto")
    a = ap.parse_args()
    from pinnacle.backend import get_backend
    from pinnacle.backend.device import resolve_device_spec
    from pinnacle.backend.optimizers import make_optimizer, optimizer_step, zero_grad
    from pinnacle.components.predictors.load_step_predictor import PlaneDirichletAnsatz
    from pinnacle.io.loader import load_input
    from pinnacle.materials.cp_field_callback import CPFieldCallback
    from pinnacle.physics.cp_residual import CPStateResidual
    from pinnacle.pinn.networks.mlp import MLP
    from pinnacle.plugins import discover_components
    from pinnacle.solvers.hybrid import build_fem_boundary_conditions
    from pinnacle.system.executioner.cp_quasi_static import _quadrature_coordinates, _residual_norms
    from pinnacle.system.executioner.fem_steady import _element_for_dim
    from pinnacle.system.problem import build_composed_problem
    from pinnacle.fem.assembly import Assembler

    discover_components()
    deck = Path(a.deck).resolve()
    cfg = load_input(str(deck))
    cfg["global"]["backend"] = "torch"
    cfg["global"]["use_gpu"] = a.device != "cpu"
    cfg["global"]["output_dir"] = str(deck.parent / "strong_form_floor")
    cfg["executioner"]["n_steps"] = a.steps
    cfg["executioner"]["t_end"] = float(cfg["executioner"]["t_end"]) * a.steps / 256.0
    cfg["executioner"]["snapshots_every"] = 0
    cfg["output"]["formats"] = ["json"]
    device = resolve_device_spec(None if a.device == "auto" else a.device, prefer_gpu=a.device != "cpu", backend_name="torch")
    backend = get_backend("torch", precision="float64", device=device, seed=0)
    import os
    os.chdir(deck.parent)
    problem = build_composed_problem(backend, cfg)
    t0 = time.perf_counter()
    result = problem.run()
    print(f"FE: {a.steps} steps on torch/{device} in {time.perf_counter() - t0:.1f} s")
    kernel = next(k for k in problem.kernels if k.__class__.__name__ == "StressDivergenceCPKernel")
    cb = kernel._material_cache[2]
    mesh = problem.mesh
    dim = int(mesh.dim)
    u_fe = np.asarray(backend.to_numpy(result["displacement"]), dtype=float).reshape(-1, dim)
    rows = result["history"]["cp_quasi_static"]
    lam = float(rows[-1]["load_factor"])
    t_end = float(cfg["executioner"]["t_end"])
    # Boundary conditions of the last step (same construction as the executioner).
    fem_var = next(v for v in problem.variables.values() if v.family == "fem")
    element = _element_for_dim(dim, order=fem_var.order)
    dirichlet_specs, neumann_specs = problem.normalized_bc_specs(dim, fem_var.name)
    scaled = []
    for spec in dirichlet_specs:
        fn = problem.functions.get(str(spec.get("function", "") or ""))
        scale = float(fn.value(lam * t_end)) if fn is not None else lam
        scaled.append({**spec, "values": [float(v) * scale for v in spec["values"]]})
    bcs = build_fem_boundary_conditions(mesh, element, backend, scaled, neumann_specs, traction_scale=lam)
    assembler = Assembler(backend, mesh, element, cb, finite_strain=bool(kernel.finite_strain))
    r_fe = _residual_norms(assembler, bcs, backend, u_fe.reshape(-1), dim)
    print(f"FE field: residual {r_fe[0]:.3e}, reaction {r_fe[1]:.3e}, ratio {r_fe[0] / r_fe[1]:.2e}")
    # Network + ansatz, regression to the FE field.
    nodes = np.asarray(mesh.nodes, dtype=float)[:, :dim]
    low, high = nodes.min(axis=0), nodes.max(axis=0); span = np.where(high - low > 0, high - low, 1.0)
    ansatz = PlaneDirichletAnsatz(nodes, dim, bcs)
    print(f"ansatz: {sum(len(v) for v in ansatz.planes.values())} plane conditions, {len(ansatz.soft)} soft")
    net = MLP(backend, input_dim=dim + 1, output_dim=dim, hidden_dims=(a.neurons,) * a.layers, activation="tanh", final_activation=None)
    opt = make_optimizer(backend, "adam", net.parameters(), lr=2e-3)
    u_scale = float(np.max(np.abs(u_fe)))
    low_t, span_t = backend.tensor(low), backend.tensor(span)
    x_nodes = backend.tensor(nodes)
    y = backend.tensor(u_fe)

    def u_fn(x):
        xn = 2.0 * (x - low_t) / span_t - 1.0
        feats = backend.concatenate([xn, backend.tensor(np.full((int(x.shape[0]), 1), lam))], axis=1)
        return ansatz.apply(backend, x, net(feats) * u_scale)

    t0 = time.perf_counter()
    for ep in range(a.epochs):
        def closure():
            zero_grad(backend, opt)
            d = (u_fn(x_nodes) - y) / u_scale
            return backend.mean(d * d)
        loss = optimizer_step(backend, opt, closure)
    u_net = np.asarray(backend.to_numpy(u_fn(x_nodes)), dtype=float)
    rel = np.linalg.norm(u_net - u_fe) / np.linalg.norm(u_fe)
    print(f"regression: {a.epochs} epochs in {time.perf_counter() - t0:.1f} s, loss {float(backend.to_numpy(loss)):.3e}, rel L2 error vs FE {rel:.3e}")
    r_net = _residual_norms(assembler, bcs, backend, u_net.reshape(-1), dim)
    print(f"regressed field: FE residual {r_net[0]:.3e}, reaction {r_net[1]:.3e}, ratio {r_net[0] / r_net[1]:.3f}")
    # Strong-form loss of the regressed field, exactly as the predictor computes it.
    xq = _quadrature_coordinates(mesh, element)
    twin = type(cb.material)(backend, cb.material.config)
    helper = CPFieldCallback(backend, twin, dim, xq.shape[0], finite_strain=bool(kernel.finite_strain), orientations=cb._orientations, batched=True)
    res = CPStateResidual(backend, twin, dim, bool(kernel.finite_strain))
    stacks = cb._current_stacks()
    x_q = backend.tensor(xq[:, :dim], requires_grad=True)
    r = res.residual(u_fn, x_q, stacks, helper._stacked_orientation())
    mu = float(getattr(twin, "mu", 1.0)); L = float(np.max(high - low)); scale = mu / L
    r_np = np.asarray(backend.to_numpy(r)) / scale
    print(f"strong-form loss of the regressed field (scale mu/L = {scale:.1f}): mean r^2 = {np.mean(r_np**2):.3f}, "
          f"p50 |r| {np.percentile(np.abs(r_np), 50):.3f}, p90 {np.percentile(np.abs(r_np), 90):.3f}, max {np.max(np.abs(r_np)):.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
