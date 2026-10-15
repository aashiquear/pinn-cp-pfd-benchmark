# W_c scan of the J2+PFD plate at l0 = 5 um (docs/10 §4.10): the calibrated
# curve at every threshold against the experiment markers and the CP+PFD
# reference; panels of Gc, softening onset and peak stress against W_c.
[Figures]
  figures       = ../../../figures/pfd
  out_dir       = ../../../results/pfd_study/j2
  prefix        = j2_wc_scan_h40
  l0            = 5
  markers       = ../../../data/digitized/al2219_fig3_experiment_markers.csv
  reference     = ../../../results/pfd_study/cp/neural_h40/base
  target_strain = 0.172
[]
[Roots]
  [./wc0]
    W_c = 0
    dir = ../../../results/pfd_study/j2/cal_h40
  []
  [./wc20]
    W_c = 20
    dir = ../../../results/pfd_study/j2/cal_h40_wc20
  []
  [./wc200]
    W_c = 200
    dir = ../../../results/pfd_study/j2/cal_h40_wc200
  []
  [./wc400]
    W_c = 400
    dir = ../../../results/pfd_study/j2/cal_h40_wc400
  []
  [./wc600]
    W_c = 600
    dir = ../../../results/pfd_study/j2/cal_h40_wc600
  []
  [./wc800]
    W_c = 800
    dir = ../../../results/pfd_study/j2/cal_h40_wc800
  []
[]
