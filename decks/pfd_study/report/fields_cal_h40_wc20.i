# Field maps of the calibrated J2+PFD (W_c 20, the carried configuration) and CP+PFD plates at h40:
# damage and equivalent plastic strain at 5, 10, 15, 17.2 (initiation) and
# 19% nominal strain, and the crack path at the end (elements with mean
# damage >= 0.8). Read by scripts/pfd_fields.py.
[Fields]
  figures    = ../../../figures/pfd
  prefix     = cal_h40_wc20
  strains    = '5 10 15 17.2 19'
  strain_max = 0.20
  threshold  = 0.8
  field      = equivalent_plastic_strain
  vmax       = 1.0
  deform     = 1.0
[]
[Runs]
  [./cp]
    dir   = ../../../results/pfd_study/cp/cal_h40/l5/gc1106.12
    label = 'CP+PFD, l0 = 5 um'
  []
  [./j2]
    dir   = ../../../results/pfd_study/j2/cal_h40_wc20/l5/gc2003.06
    label = 'J2+PFD (W_c 20), l0 = 5 um'
  []
[]
