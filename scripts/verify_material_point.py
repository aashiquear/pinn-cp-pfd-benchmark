#!/usr/bin/env python3
"""Material-point verification of the dislocation-density CP law, and the C.15 test.

Three groups of checks, all at a single quadrature point so the constitutive law
is exercised without any mesh, solver, or damage coupling in the way:

  paths      stress update under monotonic, unload-reload, biaxial and
             nonproportional (L-shaped) deformation histories;
  numerics   time-step refinement (the update is explicit Euler, so first-order
             convergence is the expectation) and the update's Jacobian against
             central finite differences;
  c15        Taylor (isostrain) polycrystal average over random orientations,
             which discriminates the two readings of the Taylor-term density
             against the digitized Al2219 calibration curve.

Units throughout: micrometre, megapascal, second (the Al2219 paper's system).
Backend: NumPy, per the measured policy in docs/02_backend_envelope.md.

Usage:
    python verify_material_point.py --group all --out ../results/
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# Al2219, independently extracted (data/extracted/benchmark_dataset.csv).
# k1 = k_m and k2 = 2*y_c per the Eq. (12) derivation; alpha0 = 0.4 per Table 2.
AL2219 = {
    "type": "dislocation_density_cp",
    "C11": 108.0e3,
    "C12": 61.3e3,
    "C44": 28.5e3,
    "mu": 25.7e3,
    "gamma0": 1.0e-5,
    "rate_exponent_n": 1.0 / 0.12,
    "tau_c0": 0.112,
    "q_lat": 1.4,
    "alpha0": 0.4,
    "burgers": 0.000286,
    "rho0": 220.0,
    "k1": 0.3432,
    "annihilation_diameter": 0.00268,
}


def make_material(backend, euler=(0.0, 0.0, 0.0), **overrides):
    from pinnacle import registry
    from pinnacle.plugins import discover_components

    discover_components()
    mat = {**AL2219, **overrides}
    return registry.create(
        "material",
        "dislocation_density_cp",
        backend=backend,
        config={
            "material": mat,
            "crystal": {"structure": "fcc", "euler_angles": list(euler)},
        },
    )


def integrate(mat, backend, F_of_t, t_end, n_steps, state=None):
    """Explicit-Euler drive of the material point along a prescribed F(t)."""
    st = {k: v for k, v in (state or mat.initial_state()).items()}
    dt = t_end / n_steps
    out = {"t": [], "F11": [], "sig11": [], "sig22": [], "sig33": [],
           "rho_mean": [], "tau_c_max": [], "gdot_absmax": []}
    for i in range(n_steps + 1):
        t = i * dt
        F = backend.tensor(np.ascontiguousarray(F_of_t(t)))
        step_dt = None if i == n_steps else dt
        cauchy, st_new, info = mat.compute(F, st, dt=step_dt)
        sig = backend.to_numpy(cauchy)
        out["t"].append(t)
        out["F11"].append(float(np.asarray(F_of_t(t))[0, 0]))
        out["sig11"].append(float(sig[0, 0]))
        out["sig22"].append(float(sig[1, 1]))
        out["sig33"].append(float(sig[2, 2]))
        out["rho_mean"].append(float(np.mean(backend.to_numpy(st["rho"]))))
        out["tau_c_max"].append(float(np.max(backend.to_numpy(info["tau_c"]))))
        out["gdot_absmax"].append(float(np.max(np.abs(backend.to_numpy(info["gamma_dot"])))))
        st = st_new
    return {k: np.asarray(v) for k, v in out.items()}, st


def F_isochoric_uniaxial(e_of_t):
    """Volume-preserving uniaxial extension: the standard large-strain Taylor drive."""

    def F(t):
        lam = 1.0 + e_of_t(t)
        lat = 1.0 / np.sqrt(lam)
        return np.diag([lam, lat, lat])

    return F


def group_paths(backend, edot=1.0e-3, e_max=0.02, n_steps=400):
    """Monotonic, unload-reload, biaxial, and nonproportional histories."""
    t_end = e_max / edot
    res = {}

    res["monotonic"] = integrate(
        make_material(backend), backend,
        F_isochoric_uniaxial(lambda t: edot * t), t_end, n_steps)[0]

    def unload_reload(t):
        # load to e_max, unload to e_max/2, reload to e_max
        f = t / t_end
        if f <= 0.5:
            return e_max * (f / 0.5)
        if f <= 0.75:
            return e_max * (1.0 - 2.0 * (f - 0.5))
        return e_max * (0.5 + 2.0 * (f - 0.75))

    res["unload_reload"] = integrate(
        make_material(backend), backend,
        F_isochoric_uniaxial(unload_reload), t_end, n_steps)[0]

    def F_biaxial(t):
        lam = 1.0 + edot * t
        return np.diag([lam, lam, 1.0 / lam**2])

    res["biaxial"] = integrate(make_material(backend), backend, F_biaxial, t_end, n_steps)[0]

    def F_nonproportional(t):
        # axial extension first, then shear at fixed axial stretch (an L-path)
        f = t / t_end
        lam = 1.0 + edot * t_end * min(f, 0.5) * 2.0 * 0.5
        g = 0.0 if f <= 0.5 else 0.5 * e_max * (f - 0.5) / 0.5
        Fm = np.diag([lam, 1.0 / np.sqrt(lam), 1.0 / np.sqrt(lam)])
        Fm[0, 1] = g
        return Fm

    res["nonproportional"] = integrate(
        make_material(backend), backend, F_nonproportional, t_end, n_steps)[0]
    return res


def group_numerics(backend, edot=1.0e-3, e_max=0.01):
    """Time-step refinement and a central-difference Jacobian check."""
    t_end = e_max / edot
    ref = None
    conv = []
    for n in (100, 200, 400, 800, 1600, 3200):
        r, st = integrate(make_material(backend), backend,
                          F_isochoric_uniaxial(lambda t: edot * t), t_end, n)
        val = (r["sig11"][-1], r["rho_mean"][-1])
        if ref is not None:
            conv.append({"n_steps": n, "sig11": val[0], "rho_mean": val[1]})
        ref = val
    # observed order from successive halvings of dt against the finest solution
    fine = conv[-1]["sig11"]
    for c in conv:
        c["sig11_err_vs_finest"] = abs(c["sig11"] - fine)
    orders = []
    for a, b in zip(conv[:-2], conv[1:-1]):
        if b["sig11_err_vs_finest"] > 0:
            orders.append(np.log2(a["sig11_err_vs_finest"] / b["sig11_err_vs_finest"]))

    # Jacobian d sigma_11 / d F_ij by central differences, at a plastically
    # active state, against the same quantity from a doubled perturbation:
    # agreement of the two step sizes is the check that the update is smooth
    # and differentiable, which is what the PINN routes require of it.
    mat = make_material(backend)
    _, st = integrate(mat, backend, F_isochoric_uniaxial(lambda t: edot * t), t_end, 400)
    F0 = np.asarray(F_isochoric_uniaxial(lambda t: edot * t)(t_end))

    def sig11(Fm):
        c, _, _ = mat.compute(backend.tensor(np.ascontiguousarray(Fm)), st, dt=None)
        return float(backend.to_numpy(c)[0, 0])

    jac = {}
    for h in (1.0e-6, 2.0e-6):
        J = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                Fp, Fm_ = F0.copy(), F0.copy()
                Fp[i, j] += h
                Fm_[i, j] -= h
                J[i, j] = (sig11(Fp) - sig11(Fm_)) / (2.0 * h)
        jac[h] = J
    d = np.max(np.abs(jac[1.0e-6] - jac[2.0e-6]))
    scale = np.max(np.abs(jac[1.0e-6]))
    return {
        "refinement": conv,
        "observed_order_sig11": orders,
        "jacobian_max_abs": float(scale),
        "jacobian_step_doubling_max_diff": float(d),
        "jacobian_relative_step_sensitivity": float(d / scale) if scale else float("nan"),
    }


def group_c15(backend, n_grains=100, edot=1.0e-3, e_max=0.02, n_steps=150, seed=0):
    """Taylor average under true uniaxial stress, for both readings of Eq. (10).

    A prescribed isochoric stretch is *not* a uniaxial-stress test: its apparent
    modulus is not Young's modulus, so a 0.2 %-offset yield taken from it is not
    comparable to a measured tensile curve. Here the macroscopic lateral stretch
    is instead solved at every increment so the aggregate transverse stress
    vanishes, which is the standard Taylor model with the transverse constraint
    relaxed. Recovering Young's modulus for aluminum from the elastic segment is
    the check that the drive is right.

    The two readings differ only in the interaction matrix entering the Taylor
    term:
      latent     a^{ab} = 1 on the diagonal, q off-diagonal -- Eq. (10) with rho
                 from Eq. (11), which is what the framework implements;
      self_only  a^{ab} = identity (q_lat = 0), so the Taylor term is
                 alpha_0 mu b sqrt(rho_SSD^a) -- the companion sheet's 43.6 MPa
                 reading.

    Taylor isostrain is an upper bound on aggregate strength, so a finite-element
    polycrystal of the same material yields at or below these values. That makes
    the bound one-sided: a reading whose Taylor yield falls *below* the published
    finite-element result is excluded, while one above it is admissible.
    """
    from pinnacle.materials.orientations import sample_orientations

    # 'mackenzie' is the Haar-uniform sampler; 'uniform' samples Euler angles
    # uniformly, which is not orientation-uniform and would bias the average.
    eul = np.asarray(sample_orientations({"type": "mackenzie"}, n_grains, seed=seed))
    dt = (e_max / edot) / n_steps
    out = {}

    for name, over in (("latent", {}), ("self_only", {"q_lat": 0.0})):
        mats = [make_material(backend, euler=tuple(eul[g]), **over) for g in range(n_grains)]
        states = [m.initial_state() for m in mats]

        def aggregate(lam_a, lam_t, sts):
            """Volume-average Cauchy stress over the aggregate at a trial F."""
            Fm = backend.tensor(np.ascontiguousarray(np.diag([lam_a, lam_t, lam_t])))
            s = np.zeros((3, 3))
            for m, st in zip(mats, sts):
                c, _, _ = m.compute(Fm, st, dt=None)
                s += backend.to_numpy(c)
            return s / len(mats)

        e_hist, s_hist, lat_hist = [], [], []
        lam_t = 1.0
        for i in range(n_steps + 1):
            e = edot * (i * dt)
            lam_a = 1.0 + e
            # Secant iteration on the lateral stretch so that the aggregate
            # transverse stress vanishes: one unknown, because a Haar-uniform
            # texture is transversely isotropic on average.
            x0 = lam_t
            f0 = 0.5 * (aggregate(lam_a, x0, states)[1, 1] + aggregate(lam_a, x0, states)[2, 2])
            x1 = x0 - 1.0e-4
            for _ in range(12):
                f1 = 0.5 * (aggregate(lam_a, x1, states)[1, 1]
                            + aggregate(lam_a, x1, states)[2, 2])
                if abs(f1) < 1.0e-6 * max(1.0, abs(f0)) or abs(f1 - f0) < 1.0e-30:
                    break
                x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
                x0, f0, x1 = x1, f1, x2
            lam_t = x1
            sig = aggregate(lam_a, lam_t, states)
            e_hist.append(e)
            s_hist.append(float(sig[0, 0]))
            lat_hist.append(lam_t)
            if i < n_steps:
                Fm = backend.tensor(np.ascontiguousarray(np.diag([lam_a, lam_t, lam_t])))
                states = [m.compute(Fm, st, dt=dt)[1] for m, st in zip(mats, states)]

        e = np.asarray(e_hist)
        sig = np.asarray(s_hist)
        el = (e > 2.0e-4) & (e < 1.0e-3)
        E = np.polyfit(e[el], sig[el], 1)[0]
        off = sig - E * (e - 0.002)
        idx = np.where((off < 0) & (e > 0.002))[0]
        m0 = make_material(backend, **over)
        out[name] = {
            "tau_c_initial_MPa": float(np.max(backend.to_numpy(
                m0.hardening.critical_stress(m0.initial_state()["rho"])))),
            "E_fit_GPa": float(E / 1000.0),
            "nu_effective": float(-(lat_hist[len(lat_hist) // 20] - 1.0)
                                  / max(e[len(e) // 20], 1e-12)),
            "sigma_y_02_MPa": float(sig[idx[0]]) if len(idx) else float("nan"),
            "sigma_at_2pct_MPa": float(np.interp(0.02, e, sig)),
            "max_abs_transverse_stress_MPa": float(abs(
                aggregate(1.0 + e_max, lam_t, states)[1, 1])),
            "strain": e.tolist(),
            "stress_MPa": sig.tolist(),
        }

    out["_meta"] = {
        "n_grains": n_grains, "strain_rate_per_s": edot, "seed": seed, "n_steps": n_steps,
        "taylor_model": ("isostrain aggregate with the transverse constraint relaxed "
                         "(lateral stretch solved for zero aggregate transverse stress); "
                         "an upper bound on aggregate strength"),
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", choices=["paths", "numerics", "c15", "all"], default="all")
    ap.add_argument("--n-grains", type=int, default=200)
    ap.add_argument("--strain-rate", type=float, default=1.0e-3)
    ap.add_argument("--out", default=".")
    args = ap.parse_args()

    from pinnacle.backend import get_backend

    backend = get_backend("numpy", precision="float64", seed=0)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    want = ["paths", "numerics", "c15"] if args.group == "all" else [args.group]

    if "paths" in want:
        r = group_paths(backend, edot=args.strain_rate)
        np.savez(outdir / "mp_paths.npz", **{f"{k}_{f}": v[f] for k, v in r.items() for f in v})
        for k, v in r.items():
            print(f"path {k:16s} sig11 final {v['sig11'][-1]:9.3f} MPa   "
                  f"rho_mean {v['rho_mean'][-1]:9.2f} um^-2   "
                  f"tau_c_max {v['tau_c_max'][-1]:8.3f} MPa   "
                  f"monotone_rho={bool(np.all(np.diff(v['rho_mean']) >= -1e-12))}")

    if "numerics" in want:
        r = group_numerics(backend, edot=args.strain_rate)
        (outdir / "mp_numerics.json").write_text(json.dumps(r, indent=1))
        print("\nrefinement (sig11 error vs finest):")
        for c in r["refinement"]:
            print(f"  n_steps {c['n_steps']:5d}  sig11 {c['sig11']:9.4f}  "
                  f"err {c['sig11_err_vs_finest']:.3e}")
        print("  observed order:", [round(o, 3) for o in r["observed_order_sig11"]])
        print(f"  jacobian max |dsig11/dF| = {r['jacobian_max_abs']:.4g} MPa, "
              f"step-doubling relative sensitivity = {r['jacobian_relative_step_sensitivity']:.3e}")

    if "c15" in want:
        r = group_c15(backend, n_grains=args.n_grains, edot=args.strain_rate)
        (outdir / "mp_c15.json").write_text(json.dumps(r, indent=1))
        print(f"\nC.15 Taylor average, {r['_meta']['n_grains']} random orientations, "
              f"edot = {r['_meta']['strain_rate_per_s']:g}/s")
        for k in ("latent", "self_only"):
            d = r[k]
            print(f"  {k:10s} tau_c(0) {d['tau_c_initial_MPa']:8.3f} MPa  "
                  f"E {d['E_fit_GPa']:6.1f} GPa  sigma_y(0.2%) {d['sigma_y_02_MPa']:8.2f} MPa  "
                  f"sigma(2%) {d['sigma_at_2pct_MPa']:8.2f} MPa")
        print("  digitized Fig. 3: simulation sigma_y(0.2%) = 363.4 MPa, "
              "experiment = 343.6 MPa, E = 72-76 GPa")
    return 0


if __name__ == "__main__":
    sys.exit(main())
