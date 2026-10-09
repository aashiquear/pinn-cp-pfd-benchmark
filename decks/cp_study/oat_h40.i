# One-at-a-time CP parameter study on the h40 plate (proposed levels, 2026-10-09).
# Each case changes ONE value of decks/cp_study/h40_base.i; the base case is
# the reference. Factors of 2 bracket the published values; tau_c0 is raised
# in absolute steps because its published value (0.112 MPa) is negligible
# against the Taylor term (176 MPa with the latent-hardened sum).
#   rho0   [um^-2]: 55, 110, 440      (base 220)   -> initial slip resistance and early hardening
#   k1     [-]    : 0.1716, 0.6864    (base 0.3432) -> hardening rate (storage)
#   y_c    [um]   : 0.00134, 0.00536  (base 0.00268) -> saturation (recovery)
#   tau_c0 [MPa]  : 10, 50            (base 0.112)  -> uniform yield offset
[Study]
  base       = h40_base.i
  out_root   = ../../results/cp_study
  markers    = ../../data/digitized/al2219_fig3_experiment_markers.csv
  concurrent = 5
[]
[Sweep]
  rho0   = '55 110 440'
  k1     = '0.1716 0.6864'
  y_c    = '0.00134 0.00536'
  tau_c0 = '10 50'
[]
