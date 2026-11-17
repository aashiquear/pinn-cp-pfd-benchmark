#!/usr/bin/env python3
"""Time one explicit step of the binary + orientation-field model in Torch (CPU and GPU).

Purpose: estimate the wall time of the large cases of (granasy 2005 spherulites) on a given
machine before any of them is run. The step is the finite-difference form of the model
(scripts/cryst_binary_reference.py) plus an orientation sub-step with the singular diffusion
term div{[1 - p(phi)] grad(theta) / |grad(theta)|} (wrapped differences, regularized modulus),
which the papers say needs 1/20 to 1/50 of the time step of the other fields.

Usage: python cryst_gpu_benchmark.py [sizes...]      (default: 500 1000 2000 4000)
Prints seconds per step for the phase + composition step and per orientation sub-step, for
every available device and float32 / float64. Each timing is bounded (a fixed number of steps)
and synchronized. Standalone: needs torch only.
"""
from __future__ import annotations

import sys
import time

import torch

A_W, B_NI, B_CU, Q_F = 10193.0, 109551.0, -143771.0, 0.273838
DX, DT = 6.25e-3, 4.75e-6


def p0(u):
    return u**3 * (10.0 - 15.0 * u + 6.0 * u * u)


def p1(u):
    return 30.0 * u * u * (1.0 - u) ** 2


def g1(u):
    return 0.5 * u * (1.0 - u) * (1.0 - 2.0 * u)


def step_phi_c(phi, c):
    lap = (torch.roll(phi, 1, 0) + torch.roll(phi, -1, 0) + torch.roll(phi, 1, 1)
           + torch.roll(phi, -1, 1) - 4.0 * phi) / DX**2
    px = (torch.roll(phi, -1, 1) - torch.roll(phi, 1, 1)) / (2 * DX)
    py = (torch.roll(phi, -1, 0) - torch.roll(phi, 1, 0)) / (2 * DX)
    ang = torch.atan2(py, px)
    s = 1.0 + 0.1 * torch.cos(4.0 * ang)
    ds = -0.4 * torch.sin(4.0 * ang)
    fx = (s * s - 1.0) * px - s * ds * py
    fy = (s * s - 1.0) * py + s * ds * px
    div = (torch.roll(fx, -1, 1) - torch.roll(fx, 1, 1) + torch.roll(fy, -1, 0) - torch.roll(fy, 1, 0)) / (2 * DX)
    b = (1.0 - c) * B_NI + c * B_CU
    rhs = lap + div - A_W * g1(phi) - p1(phi) * b
    q = Q_F * p1(phi)
    new_c = c.clone()
    for axis in (0, 1):
        phi_n, c_n, q_n = torch.roll(phi, -1, axis), torch.roll(c, -1, axis), torch.roll(q, -1, axis)
        lam = p0(0.5 * (phi + phi_n))
        cf = 0.5 * (c + c_n)
        flux = lam * ((c_n - c) - cf * (1.0 - cf) * 0.5 * (q + q_n) * (phi_n - phi)) / DX
        new_c = new_c + DT / DX * (flux - torch.roll(flux, 1, axis))
    return phi + DT * rhs, new_c


def substep_theta(theta, phi, dt_sub):
    """One explicit sub-step of theta_t = chi div{(1 - p(phi)) grad theta / |grad theta|}."""
    out = theta
    weight = 1.0 - p0(phi)
    total = torch.zeros_like(theta)
    for axis in (0, 1):
        other = 1 - axis
        d = torch.roll(theta, -1, axis) - theta
        d = d - torch.round(d)  # wrapped difference (theta is an angle in units of 2 pi / k)
        # transverse difference at the face, for the modulus
        t = 0.25 * ((torch.roll(theta, -1, other) - torch.roll(theta, 1, other))
                    + torch.roll(torch.roll(theta, -1, other) - torch.roll(theta, 1, other), -1, axis))
        t = t - torch.round(t)
        mod = torch.sqrt(d * d + t * t + 1e-12) / DX
        flux = 0.5 * (weight + torch.roll(weight, -1, axis)) * (d / DX) / mod
        total = total + (flux - torch.roll(flux, 1, axis)) / DX
    return out + dt_sub * 360.0 * total


def bench(n: int, device: str, dtype: torch.dtype, steps: int) -> tuple[float, float]:
    g = torch.Generator().manual_seed(0)
    x = torch.linspace(-1, 1, n, dtype=dtype)
    r = torch.sqrt(x[None, :] ** 2 + x[:, None] ** 2)
    phi = (0.5 * (1.0 + torch.tanh((r - 0.3) / 0.03))).to(device)
    c = (0.399112 + 0.0134 * phi.cpu()).to(device)
    theta = torch.rand((n, n), generator=g, dtype=dtype).to(device)

    def sync():
        if device != "cpu":
            torch.cuda.synchronize()

    with torch.no_grad():
        for _ in range(3):  # warm-up
            phi, c = step_phi_c(phi, c)
            theta = substep_theta(theta, phi, DT / 50)
        sync()
        t0 = time.perf_counter()
        for _ in range(steps):
            phi, c = step_phi_c(phi, c)
        sync()
        t_pc = (time.perf_counter() - t0) / steps
        t0 = time.perf_counter()
        for _ in range(steps):
            theta = substep_theta(theta, phi, DT / 50)
        sync()
        t_th = (time.perf_counter() - t0) / steps
    assert torch.isfinite(phi).all()
    return t_pc, t_th


def main() -> None:
    sizes = [int(a) for a in sys.argv[1:]] or [500, 1000, 2000, 4000]
    devices = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])
    print(f"torch {torch.__version__}; threads {torch.get_num_threads()}; "
          f"gpu {torch.cuda.get_device_name(0) if len(devices) > 1 else 'none'}")
    print("device dtype    n   s/step(phi,c)  s/substep(theta)  s/step with 20 / 50 theta substeps")
    for device in devices:
        for dtype in (torch.float32, torch.float64):
            for n in sizes:
                steps = max(3, min(40, int(4e7 / (n * n))))
                try:
                    t_pc, t_th = bench(n, device, dtype, steps)
                except RuntimeError as exc:  # out of memory on the device
                    print(f"{device:5s} {str(dtype)[6:]:8s} {n:5d}  failed: {str(exc)[:60]}")
                    continue
                print(f"{device:5s} {str(dtype)[6:]:8s} {n:5d}  {t_pc:10.4f}  {t_th:12.4f}      "
                      f"{t_pc + 20 * t_th:8.3f} / {t_pc + 50 * t_th:8.3f}", flush=True)


if __name__ == "__main__":
    main()
