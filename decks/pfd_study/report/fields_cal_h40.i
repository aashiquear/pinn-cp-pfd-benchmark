# Report field maps of the calibrated plates (docs/10 §6): damage and
# equivalent plastic strain at fixed strains, the crack path with the
# d >= 0.8 elements solid over the plastic strain, and the GIF of damage
# initiation and propagation (author's request 2026-10-14). The CP run is the
# calibrated case as rerun by the neural-route study (identical parameters,
# the series carries the CP state fields); the J2 run is the base at the W_c
# chosen by the W_c scan (docs/10 §4.10).
[Fields]
  figures    = ../../../figures/pfd
  prefix     = cal_h40
  strains    = '5 10 15 17.2 19'
  threshold  = 0.8
  field      = equivalent_plastic_strain
  vmax       = 0.5
  deform     = 1.0
  gif        = true
  gif_every  = 1
[]
[Runs]
  [./cp]
    dir   = ../../../results/pfd_study/cp/neural_h40/base
    label = 'CP+PFD, l0 = 5 um'
  []
  [./j2]
    dir   = ../../../results/pfd_study/j2/sens_h40_wc200/base
    label = 'J2+PFD, l0 = 5 um'
  []
[]
