# Stress-concentrator shape study (docs/09 §6; author's go-ahead 2026-10-14:
# run all three series and compare with the circular hole): the calibrated
# CP+PFD plate (Gc 1106.12 MPa um, l0 5 um, W_c 0) on the h40 elliptical-hole
# meshes of decks/pfd_study/mesh/h40_ell_*.i. The circle is the base case.
#   series A (transverse extent 5 um): circle, x2.5_y5, x1.25_y5
#   series B (hole area of the circle): circle, x3.54_y7.07, x2.5_y10
#   series C (5 x 2.5 ellipse rotated): x2.5_y5 (0), a45, a90
[Study]
  base       = base_h40.i
  out_root   = ../../../results/pfd_study/cp/shape_h40
  concurrent = 3
[]
[Cases]
  [./A_x2.5_y5]
    mesh = ../../../results/pfd_study/mesh/h40_ell_x2.5_y5/mesh.e
  []
  [./A_x1.25_y5]
    mesh = ../../../results/pfd_study/mesh/h40_ell_x1.25_y5/mesh.e
  []
  [./B_x3.54_y7.07]
    mesh = ../../../results/pfd_study/mesh/h40_ell_x3.54_y7.07/mesh.e
  []
  [./B_x2.5_y10]
    mesh = ../../../results/pfd_study/mesh/h40_ell_x2.5_y10/mesh.e
  []
  [./C_a45]
    mesh = ../../../results/pfd_study/mesh/h40_ell_x2.5_y5_a45/mesh.e
  []
  [./C_a90]
    mesh = ../../../results/pfd_study/mesh/h40_ell_x2.5_y5_a90/mesh.e
  []
[]
