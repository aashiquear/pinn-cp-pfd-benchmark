# Second sensitivity set on the calibrated J2+PFD plate at h40: the crack
# density functional (AT1, with an elastic threshold, against the AT2 base)
# and the plastic-work prefactor (zeta 0.5 against 1), both at the calibrated
# Gc and l0 = 5 um, to quantify the early diffuse softening seen in the
# calibrated curves (docs/10 §3.2). Run after sens_h40.
[Study]
  base       = base_h40.i
  out_root   = ../../../results/pfd_study/j2/sens2_h40
  concurrent = 2
[]
[Sweep]
  damage_model = 'at1'
  zeta         = '0.5'
[]
