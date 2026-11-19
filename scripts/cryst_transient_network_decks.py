#!/usr/bin/env python3
"""Cases A2 and A5, network routes: decks derived from the FE decks of the single-field model.

Usage: python3 scripts/cryst_transient_network_decks.py
For each of a2_front_velocity_standard, a5_seed_subcritical, a5_seed_supercritical it writes
into decks/crystallization/single_field/:
  <case>_window_fem.i        FE route on the window the network routes use (shorter domain and
                             horizon than the *_fem.i deck; same mesh size and time step)
  <case>_pinn_seed{0,1,2}.i  `pf_pinn`: one network phi_theta(x, t) on the strong residual
  <case>_hybrid_seed0.i      `pf_pinn_fem_hybrid`: that network starts every Newton step of the
                             FE march of the window deck
Materials, kernels, initial and boundary conditions are those of the FE deck.
"""
from __future__ import annotations

import re
from pathlib import Path

DIR = Path(__file__).resolve().parents[1] / "decks/crystallization/single_field"

# case: (domain end, divisions, steps of dt = 0.02, replacements in the body)
CASES = {
    "a2_front_velocity_standard": (10.0, 1000, 300, {"position  = 5.0": "position  = 2.5"}),
    "a5_seed_subcritical": (6.0, 600, 150, {}),
    "a5_seed_supercritical": (6.0, 600, 200, {}),
}
PINN = """  type                 = pf_pinn
  residual_form        = strong
  n_steps              = {steps}
  dt                   = 0.02        # horizon t = {horizon:g}; the network is continuous in time
  slabs                = {slabs}           # time-slab curriculum, one shared network
  continuity           = data
  slab_epochs          = 1000        # Adam epochs per slab
  adam_lr              = 2.0e-3
  lbfgs_iters          = 30          # then L-BFGS per slab
  causal               = false
  weighting            = fixed
  hidden_layers        = 4
  neurons              = 48
  transform_bounded    = true        # phi in [0, 1] by construction
  bound_lower          = 0.0
  bound_upper          = 1.0
  transform_ic_exact   = true        # the initial condition holds exactly
  collocation_interior = 2000
  n_grid               = 201
  seed                 = {seed}
"""
FEM = """  type            = pf_fem_transient
  scheme          = bdf2
  n_steps         = {steps}
  dt              = 0.02       # to t = {horizon:g}
  newton_tol      = 1.0e-10
  snapshots_every = 5
"""
HOW = {
    "window_fem": """# HOW IT IS SOLVED (FE route on the window of the network routes). The deck {fe_deck} on
# a shorter domain ([0, {end:g}] nm) and horizon (t = {horizon:g}), with the same mesh size
# (0.01 nm) and time step: the reference the PINN and PINN+FEM decks of this case are
# compared with, and the corrector of the PINN+FEM deck.
""",
    "pinn": """# HOW IT IS SOLVED (PINN route, `pf_pinn`). One network phi_theta({coord}, t) (4 x 48, tanh)
# is trained on the strong residual of the equation above at 2000 collocation points over
# [0, {end:g}] nm x [0, {horizon:g}], in {slabs} time slabs of one time unit (1000 Adam epochs and 30
# L-BFGS iterations each). The initial condition and the range 0 <= phi <= 1 hold by
# construction of the ansatz; the boundary values are penalised.{radial}
#
# EXPECTED RESULT. Compared with {window} (same window, FE): values in the report. The
# solid amount at the slab ends is in summary.json (`invariants`, `mass`).
""",
    "hybrid": """# HOW IT IS SOLVED (PINN+FEM route, `pf_pinn_fem_hybrid`). Predictor: the network of the
# PINN deck, trained over the horizon. Corrector: the FE march of {window} (BDF2,
# dt = 0.02); the Newton iteration of every step starts from the network's field at the end
# time of that step instead of from the previous state. Every state is the corrector's.
# solver_history.csv records per step `newton_iters`, `predictor_rms` (distance of the
# prediction from the converged state) and `predictor_mass` (the solid amount of the
# predicted field).
#
# EXPECTED RESULT. The fields of {window} to the Newton tolerance. What the network
# changes is the number of Newton iterations per step; values in the report.
""",
}


def derive(case: str, route: str, seed: int = 0) -> tuple[str, str]:
    end, divisions, steps, swaps = CASES[case]
    horizon = steps * 0.02
    slabs = int(round(horizon))
    text = (DIR / f"{case}_fem.i").read_text()
    name = f"{case}_{route}" + ("" if route == "window_fem" else f"_seed{seed}")
    head, body = text.split("[Global]", 1)
    body = "[Global]" + body
    spherical = "spherical" in body
    head = head.replace("(FE route)", {"window_fem": "(FE route, window of the network routes)",
                                        "pinn": "(PINN route)", "hybrid": "(PINN+FEM route)"}[route])
    cut = head.index("# UNITS.")
    tail = head[cut:]
    tail = tail[: tail.index("# EXPECTED RESULT.")]
    head = head[:cut] + HOW[route].format(
        fe_deck=f"{case}_fem.i", window=f"{case}_window_fem.i", end=end, horizon=horizon, slabs=slabs,
        coord="r" if spherical else "x",
        radial=("\n# On the spherical mesh the network reads r^2 (the field is even in r) and the strong "
                "form\n# carries the term (2/r) phi_r." if spherical else "")) + "#\n" + tail.rstrip("\n# ") + "\n\n"
    body = body.replace(f"single_field/{case}_fem", f"single_field/{name}")
    body = re.sub(r"(domain\s*=\s*)'0\.0 [0-9.]+'", rf"\g<1>'0.0 {end:.1f}'", body)
    body = re.sub(r"(divisions\s*=\s*)\d+(\s*#.*)?", rf"\g<1>{divisions}             # h = 0.01 nm", body)
    for old, new in swaps.items():
        assert old in body, old
        body = body.replace(old, new)
    if route != "window_fem":
        body = body.replace("backend    = numpy", "backend    = torch")
        body = re.sub(r"(  seed       = )\d+", rf"\g<1>{seed}", body)
    start, end_ = body.index("[Executioner]"), body.index("[Postprocessors]")
    pinn = PINN.format(steps=steps, horizon=horizon, slabs=slabs, seed=seed)
    fem = FEM.format(steps=steps, horizon=horizon)
    if route == "pinn":
        body = body.replace("family     = fem\n    components = 1\n    order      = 2",
                            "family     = pinn\n    components = 1\n    network    = mlp")
        start = body.index("[Executioner]")
        # the FE postprocessor needs the FE field; the network run reports its own invariants
        body = body[:start] + "[Executioner]\n" + pinn + "[]\n\n" + body[body.index("[Output]"):]
    elif route == "hybrid":
        ind = lambda t: "".join("  " + ln + "\n" for ln in t.splitlines())  # noqa: E731
        block = ("[Executioner]\n  type = pf_pinn_fem_hybrid\n  [./pinn]\n" + ind(pinn) + "  []\n  [./fem]\n"
                 + ind(fem) + "  []\n[]\n\n")
        body = body[:start] + block + body[end_:]
    else:
        body = body[:start] + "[Executioner]\n" + fem + "[]\n\n" + body[end_:]
    return name, head + body


def main() -> None:
    written = []
    for case in CASES:
        written.append(derive(case, "window_fem"))
        written += [derive(case, "pinn", s) for s in (0, 1, 2)]
        written.append(derive(case, "hybrid", 0))
    for name, text in written:
        (DIR / f"{name}.i").write_text(text)
    print("wrote", len(written), "decks:", " ".join(n for n, _ in written))


if __name__ == "__main__":
    main()
