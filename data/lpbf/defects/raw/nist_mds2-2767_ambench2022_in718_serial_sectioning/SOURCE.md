# NIST mds2-2767: AM Bench 2022 IN718 serial sectioning and XRCT (documentation only)

- Title: AM Bench 2022: IN718 Serial Sectioning and X-ray Computed Tomography Measurement Data
- Authors: E. J. Schwalbach, M. G. Chapman, M. N. Shah, M. D. Uchic, N. Hrabe, O. Kafka, N. Moser, B. Lane, R. Carson, J. Belak, L. E. Levine
- DOI: 10.18434/mds2-2767 ; listing used: https://data.nist.gov/rmm/records/mds2-2767 ; files from https://data.nist.gov/od/ds/mds2-2767/<path>
- Date accessed: 2026-11-25
- Licence (as stated in the record metadata): https://www.nist.gov/open/license

## Finding: no pore table and no pore segmentation

The record holds 8465 files, 2.04 TB. Only the README and the CONTEXT pdf were downloaded. Neither contains the words pore, porosity or void.
The only csv in the record is `Reconstructed Data with Analysis/9 - EBSD_1um - Feature_Statistics_2deg.csv` (152 MB), which per the README is a
grain (EBSD feature) statistics table, not a pore table. The segmentation pipeline `9 - EBSD Segmentation.json` segments grains.

XRCT data are grey-value only (alloy IN718, AMB2022-01 specimen, 2.5 mm leg region; cubic voxel 2.0025 um per the README):

| file(s) | count | bytes |
|---|---:|---:|
| `Section Data/XRCT_Sections/sliver_in718_stitch_16bit_0000.tif` ... `_3811.tif` | 3812 | 7 570 639 624 |
| `Reconstructed Data/0 - XRCT_stack.dream3d` (+ .json, .xdmf) | 1 | 2 418 619 146 |
| `Reconstructed Data/2 - XRCT_stack_BuildCoordSystem.dream3d` (+ .json, .xdmf) | 1 | 2 225 109 315 |
| `Reconstructed Data/1 - CT_2_Build.dream3d` (+ `1 - CT_to_Build.json`, .xdmf) | 1 | 452 770 029 |

A pore population would have to be segmented by the user from these stacks. One specimen, one process condition: no process-parameter variation.
None of these was downloaded (each XRCT product exceeds or approaches the 2 GB limit).

## Files (downloaded 2026-11-25)

| file | bytes | sha256 |
|---|---:|---|
| `CONTEXT_2767.pdf` | 12032873 | `2f5a4835e2a7aaa278383a5d307378f9a6601d6a702ed422e45bf011292c2909` |
| `README_2767_IN718SerialSectioning.txt` | 21356 | `588587f9a55a7424beacd78e17268eea9c43232e20b46f8840acdd1104636e0e` |
