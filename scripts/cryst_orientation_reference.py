#!/usr/bin/env python3
"""Independent reference for stage C: binary + orientation-field model, Torch finite differences.

Model C of (granasy 2004 modelling), appendix A.3, as used in (granasy 2005 spherulites);
phi = 1 liquid, phi = 0 crystal; theta in [0, 1) is the orientation in units of 2 pi / k.
Dimensionless equations (docs/14_phase_field_crystallization_plan.md, sections 7 and 8):

    phi_t   = m_phi(vartheta, theta) { lap(phi) + div[(s^2-1) grad phi] - d_x(s s' phi_y) + d_y(s s' phi_x)
                                       - a_w g'(phi) - p'(phi) [ b(c) - h F(u) ] } + noise
    c_t     = div{ p(phi) [ grad c - c (1-c) q_f p'(phi) grad phi ] }
    theta_t = chi_l p(phi) div{ [1 - p(phi)] pi [x F0~ m + (1-x) F1~ n] grad theta / |grad theta| } + noise
    F  = x F0 + (1-x) F1,  u = xi0 |grad theta|  (xi0 = dx)
    s  = 1 + s0 cos(k vartheta - 2 pi theta),   m_phi = 1 + delta0 cos(k vartheta - 2 pi theta)

No PINNACLE code is used. Assumed where the papers give nothing (section 8.1): the scale
h = xi^2 H / (2 xi0 eps^2) of the orientational energy in the phi equation; differences of
theta wrapped to [-1/2, 1/2]; the modulus regularized as sqrt(|grad theta|^2 + eps^2); the
noise amplitudes as standard deviations of the per-step increment; the seed.

Usage: python cryst_orientation_reference.py <deck.i>     (Torch; GPU if available)
Writes <out_dir>/frames.npz (phi, c, theta at the saved steps) and history.csv.
"""
from __future__ import annotations

import csv
import math
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
A_W, B_NI, B_CU, Q_F = 10193.2, 109551.0, -143771.0, 0.273838
C_L, C_S = 0.466219, 0.399112
DELTA = 9.90476e-3


def read_deck(path: Path) -> dict:
    out: dict = {}
    block = None
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            name = line[1:-1].strip()
            block = None if name == "" else out.setdefault(name, {})
        elif "=" in line and block is not None:
            k, v = (s.strip() for s in line.split("=", 1))
            block[k] = v.strip("'\"")
    return out


def p0(u):
    return u**3 * (10.0 - 15.0 * u + 6.0 * u * u)


def p1(u):
    return 30.0 * u * u * (1.0 - u) ** 2


def g1(u):
    return 0.5 * u * (1.0 - u) * (1.0 - 2.0 * u)


def main() -> None:
    deck = read_deck(Path(sys.argv[1]).resolve())
    grid, case = deck["Grid"], deck["Case"]
    n, dx, dt = int(grid["n"]), float(grid["dx"]), float(grid["dt"])
    steps, save_every = int(grid["steps"]), int(grid["save_every"])
    sub = int(grid.get("theta_substeps", 20))
    dtype = torch.float64 if grid.get("precision", "float64") == "float64" else torch.float32
    device = "cuda" if torch.cuda.is_available() and grid.get("device", "auto") != "cpu" else "cpu"
    s_sat = float(case["supersaturation"])
    s0, delta0, kfold = float(case.get("s0", 0.0)), float(case.get("delta0", 0.0)), int(case.get("k", 4))
    m_branch, n_branch, x_branch = float(case.get("m", 1.0)), float(case.get("n", 0.5)), float(case.get("x", 0.0))
    chi_l = float(case.get("chi_liquid", 360.0))
    h_ori = float(case["h_orientation"])
    noise_phi, noise_theta = float(case.get("noise_phi", 0.0025)), float(case.get("noise_theta", 0.0375))
    seed_radius = float(case["seed_radius_dx"]) * dx
    seed_angle = float(case.get("seed_theta", 0.0))
    gen = torch.Generator(device=device).manual_seed(int(case.get("rng_seed", 0)))
    c_inf = C_L - s_sat * (C_L - C_S)

    xs = (torch.arange(n, dtype=dtype) - n / 2 + 0.5) * dx
    rr = torch.sqrt(xs[None, :] ** 2 + xs[:, None] ** 2).to(device)
    phi = 0.5 * (1.0 + torch.tanh((rr - seed_radius) / (2.0 * math.sqrt(2.0) * DELTA)))
    c = C_S + (c_inf - C_S) * phi
    theta = torch.rand((n, n), generator=gen, dtype=dtype, device=device)
    theta = torch.where(rr < seed_radius, torch.full_like(theta, seed_angle), theta)
    eps_reg = float(case.get("theta_regularization", 1e-3)) / dx  # in units of |grad theta|

    def roll(a, shift, axis):
        return torch.roll(a, shift, axis)

    def wrap(d):
        return d - torch.round(d)

    out_dir = ROOT / case["out_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    frames = {"phi": [phi.cpu().numpy().astype(np.float32)], "c": [c.cpu().numpy().astype(np.float32)],
              "theta": [theta.cpu().numpy().astype(np.float32)]}
    frame_steps, rows = [0], []
    solute0 = float(c.sum())
    t_start = time.time()
    two_pi = 2.0 * math.pi
    with torch.no_grad():
        for step in range(1, steps + 1):
            # ---------------- phase field
            px = (roll(phi, -1, 1) - roll(phi, 1, 1)) / (2 * dx)
            py = (roll(phi, -1, 0) - roll(phi, 1, 0)) / (2 * dx)
            ang = kfold * torch.atan2(py, px) - two_pi * theta
            rhs = (roll(phi, 1, 0) + roll(phi, -1, 0) + roll(phi, 1, 1) + roll(phi, -1, 1) - 4 * phi) / dx**2
            if s0 != 0.0:
                s = 1.0 + s0 * torch.cos(ang)
                ds = -s0 * kfold * torch.sin(ang)
                fx = (s * s - 1.0) * px - s * ds * py
                fy = (s * s - 1.0) * py + s * ds * px
                rhs = rhs + (roll(fx, -1, 1) - roll(fx, 1, 1) + roll(fy, -1, 0) - roll(fy, 1, 0)) / (2 * dx)
            # orientational energy at the node: u = xi0 |grad theta| from wrapped central differences
            tx = wrap(roll(theta, -1, 1) - roll(theta, 1, 1)) / 2.0
            ty = wrap(roll(theta, -1, 0) - roll(theta, 1, 0)) / 2.0
            u = torch.sqrt(tx * tx + ty * ty)  # xi0 = dx
            f0 = torch.where(u < 0.75 / m_branch, torch.abs(torch.sin(two_pi * m_branch * u)), torch.ones_like(u))
            f1 = torch.where(u < 0.25 / n_branch, torch.abs(torch.sin(two_pi * n_branch * u)), torch.ones_like(u))
            f_ori = x_branch * f0 + (1.0 - x_branch) * f1
            b = (1.0 - c) * B_NI + c * B_CU
            rhs = rhs - A_W * g1(phi) - p1(phi) * (b - h_ori * f_ori)
            m_phi = 1.0 + delta0 * torch.cos(ang) if delta0 != 0.0 else 1.0
            new_phi = phi + dt * m_phi * rhs
            if noise_phi > 0.0:
                new_phi = new_phi + noise_phi * p0(phi.clamp(0, 1)) * torch.randn(phi.shape, generator=gen, dtype=dtype, device=device)
            # ---------------- composition
            q = Q_F * p1(phi)
            new_c = c.clone()
            for axis in (0, 1):
                phi_n, c_n, q_n = roll(phi, -1, axis), roll(c, -1, axis), roll(q, -1, axis)
                lam = p0(0.5 * (phi + phi_n))
                cf = 0.5 * (c + c_n)
                flux = lam * ((c_n - c) - cf * (1.0 - cf) * 0.5 * (q + q_n) * (phi_n - phi)) / dx
                new_c = new_c + dt / dx * (flux - roll(flux, 1, axis))
            # ---------------- orientation: sub-steps with phi frozen
            pphi = p0(phi.clamp(0, 1))
            weight = 1.0 - pphi
            chi = chi_l * pphi
            for _ in range(sub):
                div = torch.zeros_like(theta)
                for axis in (0, 1):
                    other = 1 - axis
                    d = wrap(roll(theta, -1, axis) - theta)
                    tc = wrap(roll(theta, -1, other) - roll(theta, 1, other))
                    tt = 0.25 * (tc + roll(tc, -1, axis))
                    mod = torch.sqrt(d * d + tt * tt + (eps_reg * dx) ** 2)  # xi0 |grad theta| at the face
                    f0t = torch.where(mod < 0.75 / m_branch,
                                      torch.sign(torch.sin(two_pi * m_branch * mod)) * torch.cos(two_pi * m_branch * mod),
                                      torch.zeros_like(mod))
                    f1t = torch.where(mod < 0.25 / n_branch,
                                      torch.sign(torch.sin(two_pi * n_branch * mod)) * torch.cos(two_pi * n_branch * mod),
                                      torch.zeros_like(mod))
                    coef = math.pi * (x_branch * f0t * m_branch + (1.0 - x_branch) * f1t * n_branch)
                    flux = 0.5 * (weight + roll(weight, -1, axis)) * coef * d / mod  # grad theta / |grad theta|
                    div = div + (flux - roll(flux, 1, axis)) / dx
                theta = theta + (dt / sub) * chi * div
                if noise_theta > 0.0:
                    theta = theta + (noise_theta / math.sqrt(sub)) * pphi * torch.randn(theta.shape, generator=gen, dtype=dtype, device=device)
                theta = theta - torch.floor(theta)
            phi, c = new_phi, new_c
            if device == "cuda":
                # rocm-server: a long run of queued GPU operations without a synchronization
                # point ends in a segmentation fault (measured: a loop of two-line updates
                # dies before 40,000 iterations, and runs 200,000 with a sync every 2,000)
                torch.cuda.synchronize()
            if step % save_every == 0 or step == steps:
                if not torch.isfinite(phi).all():
                    raise SystemExit(f"non-finite field at step {step}")
                solid = float((1.0 - phi).mean())
                rows.append({"step": step, "t": step * dt, "solid_fraction": solid,
                             "equivalent_radius_dx": math.sqrt(solid * n * n / math.pi),
                             "solute_drift": float(c.sum()) / solute0 - 1.0,
                             "phi_min": float(phi.min()), "phi_max": float(phi.max())})
                frames["phi"].append(phi.cpu().numpy().astype(np.float32))
                frames["c"].append(c.cpu().numpy().astype(np.float32))
                frames["theta"].append(theta.cpu().numpy().astype(np.float32))
                frame_steps.append(step)
    np.savez_compressed(out_dir / "frames.npz", phi=np.array(frames["phi"]), c=np.array(frames["c"]),
                        theta=np.array(frames["theta"]), steps=np.array(frame_steps), dx=dx, dt=dt,
                        c_inf=c_inf, c_l=C_L, c_s=C_S)
    with (out_dir / "history.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    print(f"{steps} steps on {n}x{n} ({device}, {sub} theta sub-steps) in {time.time() - t_start:.0f} s; last {rows[-1]}")


if __name__ == "__main__":
    main()
