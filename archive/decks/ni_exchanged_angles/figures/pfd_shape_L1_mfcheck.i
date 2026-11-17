# Multi-fidelity check (docs/12 §8.5): fields at failure of L.1 on h020 and h040, same carried set;
# color = exx. Read by scripts/ni_pfd_shape_figures.py.
[Figures]
  out         = ../../../figures/ni/pfd/shape_L1_mfcheck.png
  panel       = ../../../data/digitized/ni_fig5_panels/L.1_exp.png
  panel_label = 'L.1: experiment, DIC axial strain at the onset of failure (Rezwan 2026, Fig. 5)'
  exx_range   = '0.23 0.98'
  x_range     = '5.0 20.4'
  crack       = 0.85
  ridge       = true
  ridge_x     = '8.7 16.7'
  color       = exx
  dpi         = 200
[]
[Runs]
  [./h020]
    dir   = ../../../results/ni_study/pfd/wc_L1_h020/chi0_wc75/l0.6/gc0.1
    label = 'h020 (calibration run)'
  []
  [./h040]
    dir   = ../../../results/ni_study/pfd/mfcheck_L1_h040
    label = 'h040 (same set)'
  []
[]
