#!/usr/bin/env python3
"""Independent reference for stage C: binary phase-field model, explicit finite differences.

No PINNACLE code is used. The model, the dimensionless equations and every parameter with its
provenance are in docs/14_phase_field_crystallization_plan.md, section 7 (model of
(granasy 2005 spherulites), Eq. 2.3-2.4; phi = 1 liquid, phi = 0 crystal; c = Cu fraction):

    phi_t = m_phi { lap(phi) + div[(s^2 - 1) grad phi] - d_x(s s' phi_y) + d_y(s s' phi_x)
                    - a_w g'(phi) - p'(phi) b(c) } + noise
    c_t   = div{ p(phi) [ grad c - c (1 - c) (q_f p'(phi) + q_w g'(phi)) grad phi ] }

with the orientation uniform (single crystal). Square periodic grid, forward Euler; the
isotropic part of the gradient term uses the five-point Laplacian, the anisotropic remainder
central differences; the solute flux is evaluated on cell faces, so the solute is conserved to
round-off.

Usage: python3 scripts/cryst_binary_reference.py <deck.i>
The deck gives grid, steps, supersaturation, anisotropy, seed and output cadence; the script
writes <out_dir>/frames.npz (phi and c at the saved steps) and history.csv.
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from cryst_reference import read_deck  # noqa: E402

# ---- model constants (plan section 7.2) ------------------------------------------------
R_GAS = 8.314462618
T = 1574.0
XI = 2.1e-6  # m
DELTA = 20.8e-9  # m, tanh parameter of the interface after the x20.8 rescaling
V_M = 7.4e-6  # m^3/mol
NI = {"Tf": 1728.0, "L": 2350e6, "gamma": 0.315 / 6.0}
CU = {"Tf": 1358.0, "L": 1728e6, "gamma": 0.247 / 6.0}


def constants() -> dict:
    out = {}
    for name, m in (("ni", NI), ("cu", CU)):
        out[f"eps2_{name}"] = 6.0 * np.sqrt(2.0) * m["gamma"] * DELTA / m["Tf"]
        out[f"w_{name}"] = 12.0 * m["gamma"] / (np.sqrt(2.0) * DELTA * m["Tf"])
        out[f"df_{name}"] = m["L"] * (m["Tf"] - T) / m["Tf"]
    eps2 = 0.5 * (out["eps2_ni"] + out["eps2_cu"])
    out["a_w_ni"] = XI**2 * out["w_ni"] / eps2
    out["a_w_cu"] = XI**2 * out["w_cu"] / eps2
    out["b_ni"] = XI**2 * out["df_ni"] / (eps2 * T)
    out["b_cu"] = XI**2 * out["df_cu"] / (eps2 * T)
    out["q_f"] = V_M * (out["df_ni"] - out["df_cu"]) / (R_GAS * T)
    out["q_w"] = V_M * (out["w_ni"] - out["w_cu"]) / R_GAS
    # common tangent of the ideal solution
    e_cu = np.exp(V_M * out["df_cu"] / (R_GAS * T))
    e_ni = np.exp(V_M * out["df_ni"] / (R_GAS * T))
    c_l = (1.0 - e_ni) / (e_cu - e_ni)
    out["c_l"], out["c_s"] = c_l, c_l * e_cu
    out["delta_scaled"] = DELTA / XI
    return out


def g1(u):
    return 0.5 * u * (1.0 - u) * (1.0 - 2.0 * u)


def p0(u):
    return u**3 * (10.0 - 15.0 * u + 6.0 * u**2)


def p1(u):
    return 30.0 * u**2 * (1.0 - u) ** 2


def run(deck: dict, out_dir: Path) -> None:
    k = constants()
    grid = deck["Grid"]
    n = int(grid["n"])
    dx = float(grid["dx"])
    dt = float(grid["dt"])
    steps = int(grid["steps"])
    save_every = int(grid["save_every"])
    case = deck["Case"]
    s_sat = float(case["supersaturation"])
    s0 = float(case.get("s0", 0.0))
    kfold = int(case.get("k", 4))
    theta0 = float(case.get("crystal_angle_deg", 0.0)) * np.pi / 180.0
    r_seed = float(case["seed_radius_dx"]) * dx
    noise_phi = float(case.get("noise_phi", 0.0))
    rng = np.random.default_rng(int(case.get("rng_seed", 0)))
    c_inf = k["c_l"] - s_sat * (k["c_l"] - k["c_s"])

    x = (np.arange(n) - n / 2 + 0.5) * dx
    xx, yy = np.meshgrid(x, x, indexing="xy")
    rr = np.hypot(xx, yy)
    width = 2.0 * np.sqrt(2.0) * k["delta_scaled"]  # phi = (1 + tanh(x / width)) / 2
    phi = 0.5 * (1.0 + np.tanh((rr - r_seed) / width))
    c = k["c_s"] + (c_inf - k["c_s"]) * phi
    solute0 = float(c.sum())

    def dxc(a):  # central differences, periodic
        return (np.roll(a, -1, 1) - np.roll(a, 1, 1)) / (2 * dx)

    def dyc(a):
        return (np.roll(a, -1, 0) - np.roll(a, 1, 0)) / (2 * dx)

    def lap(a):
        return (np.roll(a, 1, 0) + np.roll(a, -1, 0) + np.roll(a, 1, 1) + np.roll(a, -1, 1) - 4 * a) / dx**2

    out_dir.mkdir(parents=True, exist_ok=True)
    frames_phi, frames_c, frame_steps = [phi.astype(np.float32)], [c.astype(np.float32)], [0]
    rows = []
    t0 = time.time()
    for step in range(1, steps + 1):
        # ---- phase field
        rhs = lap(phi)
        if s0 != 0.0:
            px, py = dxc(phi), dyc(phi)
            ang = np.arctan2(py, px)
            s = 1.0 + s0 * np.cos(kfold * (ang - theta0))
            ds = -s0 * kfold * np.sin(kfold * (ang - theta0))
            fx = (s * s - 1.0) * px - s * ds * py
            fy = (s * s - 1.0) * py + s * ds * px
            rhs = rhs + dxc(fx) + dyc(fy)
        a_w = (1.0 - c) * k["a_w_ni"] + c * k["a_w_cu"]
        b = (1.0 - c) * k["b_ni"] + c * k["b_cu"]
        rhs = rhs - a_w * g1(phi) - p1(phi) * b
        # ---- concentration: fluxes on the faces between i and i+1
        q = k["q_f"] * p1(phi) + k["q_w"] * g1(phi)
        new_c = c.copy()
        for axis in (0, 1):
            phi_n, c_n, q_n = np.roll(phi, -1, axis), np.roll(c, -1, axis), np.roll(q, -1, axis)
            lam = p0(0.5 * (phi + phi_n))
            cf = 0.5 * (c + c_n)
            flux = lam * ((c_n - c) - cf * (1.0 - cf) * 0.5 * (q + q_n) * (phi_n - phi)) / dx
            new_c += dt / dx * (flux - np.roll(flux, 1, axis))
        new_phi = phi + dt * rhs
        if noise_phi > 0.0:
            new_phi += noise_phi * p0(np.clip(phi, 0, 1)) * rng.standard_normal(phi.shape)
        phi, c = new_phi, new_c
        if step % save_every == 0 or step == steps:
            solid = float(np.mean(1.0 - phi))
            rows.append({"step": step, "t": step * dt, "solid_fraction": solid,
                         "equivalent_radius_dx": float(np.sqrt(solid * n * n / np.pi)),
                         "solute_drift": float(c.sum() / solute0 - 1.0),
                         "phi_min": float(phi.min()), "phi_max": float(phi.max()),
                         "c_min": float(c.min()), "c_max": float(c.max())})
            frames_phi.append(phi.astype(np.float32))
            frames_c.append(c.astype(np.float32))
            frame_steps.append(step)
            if not np.isfinite(phi).all():
                raise SystemExit(f"non-finite field at step {step}")
    np.savez_compressed(out_dir / "frames.npz", phi=np.array(frames_phi), c=np.array(frames_c),
                        steps=np.array(frame_steps), dx=dx, dt=dt, c_inf=c_inf, c_l=k["c_l"], c_s=k["c_s"])
    with (out_dir / "history.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    print(f"{steps} steps on {n}x{n} in {time.time() - t0:.0f} s; c_l {k['c_l']:.6f} c_s {k['c_s']:.6f} "
          f"c_inf {c_inf:.6f}; last: {rows[-1]}")


def main() -> None:
    deck_path = Path(sys.argv[1]).resolve()
    deck = read_deck(deck_path)
    if len(sys.argv) > 2 and sys.argv[2] == "constants":
        for key, value in constants().items():
            print(f"{key:14s} {value:.6g}")
        return
    run(deck, ROOT / deck["Case"]["out_dir"])


if __name__ == "__main__":
    main()
