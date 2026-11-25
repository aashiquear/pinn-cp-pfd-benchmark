# Part V, network route on case K1: design of the training runs of the parametric
# surrogate of the single-track response.
#
# Read by scripts/lpbf_k1_surrogate.py (design | decks | run | targets). Thirty tracks are
# drawn by Latin hypercube sampling over laser power, scan speed and beam diameter, in a
# box that contains the seven NIST cases (245-325 W, 800-1200 mm/s, 49-82 um) with a
# margin. Each is one PINNACLE run of the K1 template at 10 um with the calibrated source
# (23 s on the GPU host). The seven NIST cases are NOT training data: they are the rows
# the trained surrogate is asked for.

[Design]
  template       = ../k1_template.i
  samples        = 30
  seed           = 20261125
  power_W        = '230.0 340.0'
  speed_mm_s     = '750.0 1250.0'
  d4sigma_um     = '45.0 86.0'
  mesh_um        = 10.0
  backend        = torch
  use_gpu        = true
  absorptivity   = 0.7724
  penetration_um = 45.01
  end_time       = 1.2e-3
  case_dir       = ../cases/surrogate_h10
  out_dir        = ../../../../../results/lpbf/in718/k1/surrogate_h10
  nist_targets   = ../../../../../data/lpbf/in718/k1_single_tracks.csv
  profile_points = 24            # peak temperature along the depth below the track axis
  profile_depth_um = 230.0
[]
