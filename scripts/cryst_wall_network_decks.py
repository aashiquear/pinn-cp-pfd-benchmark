#!/usr/bin/env python3
"""Case B1, network routes: decks derived from the FE decks of scripts/cryst_wall_decks.py.

Usage: python3 scripts/cryst_wall_network_decks.py
Writes into decks/crystallization/walls/, for Model I (`b1_wall_model1_dT040`) and Model II
(`b1_wall_model2_dT040_nucleus`):
  *_pinn_seed{0,1,2}.i     `pf_pinn_stationary residual_form = fem`: a network of (x, y) trained
                           on the assembled weak residual (wall flux and symmetry plane included)
  *_hybrid_seed{0,1,2}.i   `pf_pinn_fem_hybrid`: the same network, trained briefly, starts a
                           PLAIN Newton solve (no volume constraint)
  *_newton.i               control: the plain Newton solve from the same circular cap
and for Model II at 20, 60 and 80 K the hybrid (seed 0) and the control; at 20 K also the
hybrid with 150 training iterations (`*_hybrid150_seed0.i`).
Mesh, initial condition, materials, kernels, boundary conditions and postprocessors are those
of the FE deck, so the three routes are compared on one discretisation.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "decks/crystallization/walls"

PINN = """  type               = pf_pinn_stationary
  residual_form      = fem         # loss: discrete L2 norm of the assembled weak residual
  hidden_layers      = 4
  neurons            = 48
  pretrain_epochs    = 3000        # regression on the circular cap (the initial condition)
  lbfgs_iters        = {iters}
  deflation_power    = 1.0         # the state without a nucleus is not a minimum of the loss
{deflation}  collapse_tol       = 0.1
  seed               = {seed}
"""
NEWTON = """  type            = pf_fem_transient
  scheme          = bdf1
  n_steps         = 1
  dt              = 1.0e12        # one backward-Euler step: plain Newton on the stationary equation
  max_newton_iter = 60
  newton_tol      = 1.0e-9
  max_cutbacks    = 0
  snapshots_every = 1
"""
HOW = {
    "pinn": """# HOW IT IS SOLVED (PINN route, `pf_pinn_stationary residual_form = fem`). A network
# phi_theta(x, y) (4 x 48, tanh) is fitted to the circular cap of the classical radius and then
# trained (L-BFGS) on the discrete L2 norm of the weak residual the FE arm assembles for its
# nodal values on the mesh of the FE deck: loss = sum_i R_i^2 / m_i / sum_i m_i (m_i: nodal
# measure), with the exact adjoint of the consistent Jacobian. The wall condition, the symmetry
# plane and the Dirichlet values are inside R; no strong form and no boundary penalty is used.
# The nucleus is a saddle of F but a zero, hence a minimum, of the loss; the state without a
# nucleus is deflated ({deflated}).
#
# EXPECTED RESULT. The free energy of the FE deck ({fe_deck}) within 1e-3 and the field
# within about 0.02 of it; values and seed scatter in the report.
""",
    "hybrid": """# HOW IT IS SOLVED (PINN+FEM route, `pf_pinn_fem_hybrid`). Predictor: the network of the
# PINN route, trained for 40 L-BFGS iterations only. Corrector: one backward-Euler step with a
# very large time step, i.e. PLAIN Newton on the weak stationary equation (no volume
# constraint, unlike the FE deck), started from the network's nodal field. The result is the
# corrector's. The run reports the distance between prediction and corrected field.
#
# EXPECTED RESULT. The state of the FE deck ({fe_deck}, found there by the constrained
# saddle search) to the Newton tolerance. The control deck `{control}` shows what plain Newton
# does from the circular cap without the network.
""",
    "newton": """# HOW IT IS SOLVED (control). PLAIN Newton on the weak stationary equation (one backward-Euler
# step with a very large time step) from the circular cap of the classical radius: no volume
# constraint and no network. This is the solve the hybrid deck corrects with.
#
# EXPECTED RESULT. Recorded in the report: for Model I the nucleus of the FE deck
# ({fe_deck}); for Model II another stationary state or no convergence, which is why the FE
# deck uses `pf_fem_saddle`.
""",
}


def derive(fe_name: str, route: str, seed: int = 0) -> tuple[str, str]:
    text = (DIR / f"{fe_name}.i").read_text()
    stem = fe_name[: -len("_fem")]
    name = f"{stem}_{route}" + ("" if route == "newton" else f"_seed{seed}")
    model2 = "model2" in fe_name
    head, body = text.split("[Global]", 1)
    body = "[Global]" + body
    head = head.replace("(FE route)", {"pinn": "(PINN route)", "hybrid": "(PINN+FEM route)",
                                        "newton": "(FE control, plain Newton)"}[route])
    cut = head.index("# HOW IT IS SOLVED")
    head = head[:cut] + HOW[route].format(
        fe_deck=f"{fe_name}.i", control=f"{stem}_newton.i",
        deflated="the wall layer, `deflation_function = layer`" if model2 else "the liquid phi = 0",
    ) + "\n"
    body = body.replace(f"walls/{fe_name}", f"walls/{name}")
    if route != "newton":
        body = body.replace("backend    = numpy", "backend    = torch")
        body = re.sub(r"(  seed       = )\d+", rf"\g<1>{seed}", body)
    if model2 and route != "newton":
        body = body.replace("  # window of the constrained amount", "  # the wall layer of the "
                            "melting point: the stationary state without a nucleus (deflated)\n"
                            "  [./layer]\n    kind       = parsed\n"
                            "    expression = 1/(1+0.5*exp(y))\n  []\n"
                            "  # window of the constrained amount (used by the FE deck only)", 1)
    deflation = "  deflation_function = layer\n" if model2 else ""
    start, end = body.index("[Executioner]"), body.index("[Postprocessors]")
    if route == "pinn":
        block = "[Executioner]\n" + PINN.format(iters=150, deflation=deflation, seed=seed) + "[]\n\n"
    elif route == "hybrid":
        inner = "".join("  " + line + "\n" for line in PINN.format(
            iters=40, deflation=deflation, seed=seed).splitlines())
        fem = "".join("  " + line + "\n" for line in NEWTON.splitlines())
        block = ("[Executioner]\n  type = pf_pinn_fem_hybrid\n  [./pinn]\n" + inner + "  []\n"
                 "  [./fem]\n" + fem + "  []\n[]\n\n")
    else:
        block = "[Executioner]\n" + NEWTON + "[]\n\n"
    return name, head + body[:start] + block + body[end:]


def main() -> None:
    written = []
    for fe_name in ("b1_wall_model1_dT040_fem", "b1_wall_model2_dT040_nucleus_fem"):
        for route in ("pinn", "hybrid"):
            for seed in (0, 1, 2):
                written.append(derive(fe_name, route, seed))
        written.append(derive(fe_name, "newton"))
    for d_t in (20, 60, 80):
        fe_name = f"b1_wall_model2_dT{d_t:03d}_nucleus_fem"
        written.append(derive(fe_name, "hybrid", 0))
        written.append(derive(fe_name, "newton"))
    # 20 K: the 40-iteration predictor leaves Newton outside the basin; the same deck with the
    # training length of the PINN route
    name, text = derive("b1_wall_model2_dT020_nucleus_fem", "hybrid", 0)
    written.append((name.replace("_hybrid_", "_hybrid150_"),
                    text.replace(name, name.replace("_hybrid_", "_hybrid150_"))
                    .replace("lbfgs_iters        = 40", "lbfgs_iters        = 150")
                    .replace("trained for 40 L-BFGS iterations only", "trained for 150 L-BFGS iterations")))
    for name, text in written:
        (DIR / f"{name}.i").write_text(text)
    print("wrote", len(written), "decks:", " ".join(n for n, _ in written))


if __name__ == "__main__":
    main()
