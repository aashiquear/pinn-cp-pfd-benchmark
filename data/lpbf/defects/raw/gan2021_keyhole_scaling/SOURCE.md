# Gan et al. 2021, Nature Communications: keyhole scaling laws, Supplementary Data 1 and 2

- Title: Universal scaling laws of keyhole stability and porosity in 3D printing of metals
- Authors: Z. Gan, O. L. Kafka, N. Parab, C. Zhao, L. Fang, O. Heinonen, T. Sun, W. K. Liu
- Nature Communications 12 (2021) 2379, DOI 10.1038/s41467-021-22704-0, https://www.nature.com/articles/s41467-021-22704-0
- Files from https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-22704-0/MediaObjects/41467_2021_22704_MOESM{1,3}_ESM.pdf and MOESM{4,5}_ESM.xlsx
- Date accessed: 2026-11-25
- Licence: article is open access under CC BY 4.0 (stated on the article page)

## Content (checked by opening the workbooks)

- `MOESM4_ESM.xlsx` = Supplementary Data 1 ("all process parameters, material properties, and measured keyhole dimensions"): one sheet, 87 data rows, 16 columns:
  P (103-522), eta, V (0.3-1.2), d, r0, K, rho, cp, alpha, Tv-T0, Lv, Tl-T0, Lm, e, `Length (µm)` (27 rows empty), `tan(theta)`.
  Groups (label in the first column, first row of each block only): `Ti64, spot size 95 μm` 48 rows, `Ti6, spot size 140 μm` 23 rows (label as written),
  `Al6061, spot size 88 μm` 9 rows, `SS316, spot size 110 μm` 7 rows.
  Units are not written in the file; values are consistent with SI (P in W, V in m/s, d and r0 in m; e = keyhole depth in m per the paper's notation, range 1.5e-5 to 5.3e-4).
  The exact definition of `Length (µm)` and `tan(theta)` is not stated in the file.
- `MOESM5_ESM.xlsx` = Supplementary Data 2 ("transient powers with different process parameters and material properties"): one sheet, 28 data rows, 21 columns:
  eta(P-P0), V, d, r0, material properties, e, eta, eta_reflect, eta_convect, eta_rediate, eta_evap, eta_spatt. Groups: Al6061 88 um 8 rows, Ti64 95 um 15 rows, Ti64 140 um 5 rows (all marked EXP).
- `MOESM1_ESM.pdf` Supplementary Information (Supplementary Tables 1-3: material properties, simulation constants, dimensionless numbers); `MOESM3_ESM.pdf` description of the additional files.

## Caveats

- These are bare-plate/powder synchrotron X-ray imaging experiments (single tracks), not a commercial LPBF machine; no layer thickness column.
- Neither workbook has a porosity column or a pore/no-pore label. The keyhole-porosity threshold of the paper (Keyhole number) is in the text/figures only.
- Supplementary Movies 1-11 (mp4) were not downloaded.

## Files (downloaded 2026-11-25)

| file | bytes | sha256 |
|---|---:|---|
| `41467_2021_22704_MOESM1_ESM.pdf` | 1519413 | `22a8d2dc870bbd83f6bd3330899437911733f9335d52758d95a328bd6993d6e5` |
| `41467_2021_22704_MOESM3_ESM.pdf` | 336575 | `1be86f67ed5d0b333f0c72b05861f7073656db3a17a69f11b982bd109476c1e3` |
| `41467_2021_22704_MOESM4_ESM.xlsx` | 18260 | `b52d9173eb3a63983dd87c56b2720abd419a200e9594b863218de13da6a7597a` |
| `41467_2021_22704_MOESM5_ESM.xlsx` | 13859 | `8040341d67d3c9390a6e1eff209cee059d4e922fa6fc1519e94941a1de84c763` |
