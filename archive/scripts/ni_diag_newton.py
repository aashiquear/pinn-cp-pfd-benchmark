"""Run a PINNACLE deck with per-iteration Newton diagnostics (docs/12 §8.7, 2026-11-04).

Diagnoses the softening stall of the W_c = 0 baseline (L.3: the residual of a vanishing
increment levels off just above newton_tol and Newton walks away from it). Logged into
the run's own solver/pinnacle.log through the logger ``pinnacle.diag``, so the lines sit
in order with the step messages:

- ``newton``: every global Newton iteration, its residual norm, the line-search step
  taken and whether that step lowered the residual (Armijo) or was the fallback.
- ``branch``: every stress evaluation of the implicit CP update with damage. Points with
  g(d) < 0.5 that sit within 1e-3 of Je = 1, and how many of them changed side of
  Je = 1 since the previous evaluation. degrade_pk2 degrades the volumetric stress only
  for Je >= 1, so the volumetric tangent jumps from g K to K there.

- ``check`` (2026-11-04, with NI_DIAG_CHECK="15 40"): at those Newton iterations of a
  solve, whether the Newton step is consistent with the residual: globally,
  (R(u + eps du) - R(u))/eps + R(u) -> 0 since K du = -R, and per quadrature point
  sigma(H + eps dH) - sigma(H) against eps C:dH; and whether R(u) repeats when the
  same point is evaluated again (warm-started and cold local solves).

Usage: [NI_DIAG_CHECK="15 40"] ni_diag_newton.py DECK.i   (run from the deck's folder)
"""
from __future__ import annotations

import logging
import os
import sys

import numpy as np

from pinnacle.cli.main import main as pinnacle_main
from pinnacle.materials import cp_field_callback as cpf
from pinnacle.solvers import nonlinear

LOG = logging.getLogger("pinnacle.diag")
STATE: dict = {"side": None, "it": 0}
# NI_DIAG_CHECK="15 40": run the tangent check at these Newton iterations of every solve
CHECK_AT = {int(v) for v in os.environ.get("NI_DIAG_CHECK", "").split()}


def _patch_newton() -> None:
    orig_solve = nonlinear.NewtonSolver.solve
    orig_backtrack = nonlinear.NewtonSolver._backtrack
    orig_probe = nonlinear.NewtonSolver._probe_norm

    def solve(self, u0):
        STATE["it"] = 0
        LOG.info("newton start tol %.3e", self.tol)
        return orig_solve(self, u0)

    def probe(self, trial, ladder=False):
        norm = orig_probe(self, trial, ladder)
        STATE["probes"].append(norm)
        return norm

    def backtrack(self, u, du, norm_r):
        if CHECK_AT and STATE["it"] in CHECK_AT:
            try:
                _tangent_check(self, u, du)
            except Exception as exc:  # diagnostics must not stop the run
                LOG.info("check unavailable: %s", exc)
        STATE["probes"] = []
        try:
            step = orig_backtrack(self, u, du, norm_r)
        except nonlinear.ConvergenceError:
            LOG.info("newton it %d |R| %.4e no feasible trial", STATE["it"], norm_r)
            raise
        alpha = float(np.linalg.norm(step) / max(np.linalg.norm(du), 1e-300))
        probes = STATE["probes"]
        armijo = any(p is not None and p <= (1.0 - 1.0e-4 * a) * norm_r
                     for p, a in zip(probes, [2.0 ** -k for k in range(len(probes))]))
        best = min((p for p in probes if p is not None), default=float("nan"))
        LOG.info("newton it %d |R| %.4e alpha %.4g %s best probe %.4e |du| %.3e",
                 STATE["it"], norm_r, alpha, "armijo" if armijo else "FALLBACK", best,
                 float(np.abs(du).max()))
        STATE["it"] += 1
        return step

    nonlinear.NewtonSolver.solve = solve
    nonlinear.NewtonSolver._probe_norm = probe
    nonlinear.NewtonSolver._backtrack = backtrack


def _free_residual(solver, u):
    be = solver.backend
    R = np.asarray(be.to_numpy(solver.assembler.residual(be.tensor(u, dtype="float64"))),
                   dtype=float).copy()
    solver._residual_norm_at(u)          # builds the fixed-DOF mask
    R[solver._fixed_mask] = 0.0
    return R


def _tangent_check(solver, u, du):
    """Is the Newton step consistent with the residual?  K du = -R on the free DOFs, so
    (R(u + eps du) - R(u))/eps + R(u) -> 0 as eps -> 0 when the tangent is the residual's
    derivative.  Global, then per quadrature point on the material: sigma(H + eps dH) -
    sigma(H) against eps C:dH."""
    asm, be = solver.assembler, solver.backend
    mat = asm.material
    R0 = _free_residual(solver, u)
    n0 = float(np.linalg.norm(R0))
    for eps in (1e-2, 1e-4, 1e-6):
        Re = _free_residual(solver, solver.bcs.enforce_values(u + eps * du, dim=asm.dim))
        mis = (Re - R0) / eps + R0
        LOG.info("check global eps %.0e |R0| %.4e |mismatch|/|R0| %.4e |dR/eps| %.4e",
                 eps, n0, float(np.linalg.norm(mis)) / n0, float(np.linalg.norm(Re - R0)) / eps)
    # the assembled tangent itself: K du against -R (the linear solve) and against the
    # finite-difference change of R (the tangent)
    u_t = be.tensor(u, dtype="float64")
    asm.residual(u_t)                    # the tangent is taken at the stress of u
    K = asm.tangent(u_t, method=solver.tangent_method, format="coo").to_scipy().tocsr()
    Kdu = np.asarray(K @ du, dtype=float)
    Kdu[solver._fixed_mask] = 0.0
    Re = _free_residual(solver, solver.bcs.enforce_values(u + 1e-4 * du, dim=asm.dim))
    LOG.info("check K |Kdu+R0|/|R0| %.4e |dR/eps-Kdu|/|R0| %.4e |Kdu| %.4e",
             float(np.linalg.norm(Kdu + R0)) / n0,
             float(np.linalg.norm((Re - R0) / 1e-4 - Kdu)) / n0, float(np.linalg.norm(Kdu)))
    # repeatability: the same point twice (warm start from itself, then cold)
    Ra = _free_residual(solver, u)
    mat._warm_enabled = False
    try:
        Rc = _free_residual(solver, u)
    finally:
        mat._warm_enabled = True
    LOG.info("check repeat |R(u)-R(u)| warm %.3e cold %.3e",
             float(np.linalg.norm(Ra - R0)), float(np.linalg.norm(Rc - R0)))
    # per quadrature point
    eps = 1e-4
    H0 = asm._compute_strain(asm._u_to_nodes(be.tensor(u, dtype="float64")))
    H1 = asm._compute_strain(asm._u_to_nodes(be.tensor(u + eps * du, dtype="float64")))
    s0 = np.asarray(be.to_numpy(mat.stress(H0)), dtype=float)
    dg0 = None if mat._last_dgamma is None else np.asarray(mat._last_dgamma).copy()
    C = np.asarray(be.to_numpy(mat.tangent(H0)), dtype=float)
    s1 = np.asarray(be.to_numpy(mat.stress(H1)), dtype=float)
    dH = (np.asarray(be.to_numpy(H1), dtype=float) - np.asarray(be.to_numpy(H0), dtype=float))
    lin = np.einsum("eqacbd,eqdb->eqac", C, dH)
    lin2 = np.einsum("eqacbd,eqbd->eqac", C, dH)
    fd = s1 - s0
    err = np.linalg.norm((fd - lin).reshape(-1, 9), axis=1)
    err2 = np.linalg.norm((fd - lin2).reshape(-1, 9), axis=1)
    LOG.info("check contraction C:dH[d,b] sum err %.4e, C:dH[b,d] sum err %.4e",
             float(err.sum()), float(err2.sum()))
    if err2.sum() < err.sum():
        err, lin = err2, lin2
    ref = np.linalg.norm(fd.reshape(-1, 9), axis=1)
    g = mat._g_per_qp()
    g = np.ones_like(err) if g is None else np.asarray(be.to_numpy(g), dtype=float).reshape(-1)
    LOG.info("check qp sum|fd-C:dH| %.4e sum|fd| %.4e points err>0.1|fd| %d of %d",
             float(err.sum()), float(ref.sum()), int(np.count_nonzero(err > 0.1 * ref + 1e-30)),
             err.size)
    order = np.argsort(err)[::-1][:12]
    for q in order:
        act = -1 if dg0 is None else int(np.count_nonzero(np.abs(dg0[q]) > 1e-12))
        mx = float("nan") if dg0 is None else float(np.abs(dg0[q]).max())
        LOG.info("check qp %d err %.3e |fd| %.3e |C:dH| %.3e g %.4f active %d max|dgamma| %.3e",
                 int(q), float(err[q]), float(ref[q]),
                 float(np.linalg.norm(lin.reshape(-1, 9)[q])), float(g[q]), act, mx)
    for lo, hi in ((0.0, 0.1), (0.1, 0.5), (0.5, 0.9), (0.9, 1.01)):
        m = (g >= lo) & (g < hi)
        LOG.info("check g in [%.1f,%.2f): points %d sum err %.3e sum |fd| %.3e", lo, hi,
                 int(m.sum()), float(err[m].sum()), float(ref[m].sum()))


def _patch_branch() -> None:
    orig_integrate = cpf.CPFieldCallback._integrate_step

    def integrate(self, F3, stacks0, dt, gb, g_qp, F0, *a, **k):
        out = orig_integrate(self, F3, stacks0, dt, gb, g_qp, F0, *a, **k)
        if g_qp is None:
            return out
        try:
            be = self.backend
            g = np.asarray(be.to_numpy(g_qp), dtype=float).reshape(-1)
            F = np.asarray(be.to_numpy(F3), dtype=float).reshape(-1, 3, 3)
            Fp = np.asarray(be.to_numpy(out["Fp"]), dtype=float).reshape(-1, 3, 3)
            Je = np.linalg.det(F) / np.linalg.det(Fp)
        except Exception as exc:  # diagnostics must not stop the run
            LOG.info("branch unavailable: %s", exc)
            return out
        dmg = g < 0.5
        side = Je >= 1.0
        prev = STATE["side"]
        flips = -1
        if prev is not None and prev.shape == side.shape:
            flips = int(np.count_nonzero((prev != side) & dmg))
        STATE["side"] = side
        if np.any(dmg):
            LOG.info("branch damaged %d near(1e-3) %d compressed %d flips %d min|Je-1| %.2e",
                     int(dmg.sum()), int(np.count_nonzero(dmg & (np.abs(Je - 1.0) < 1e-3))),
                     int(np.count_nonzero(dmg & ~side)), flips,
                     float(np.abs(Je[dmg] - 1.0).min()))
        return out

    cpf.CPFieldCallback._integrate_step = integrate


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    _patch_newton()
    _patch_branch()
    sys.exit(pinnacle_main(["run", sys.argv[1]]))
