# NIST mds2-3083: IN718 as-built process-structure-property data (partial download)

- Title: Process-structure-properties investigations for laser powder bed fused IN718 in the as-built condition
- Authors: D. S. Watring, J. T. Benzing, O. L. Kafka, L.-A. Liew, N. H. Moser, N. Hrabe, A. D. Spear
- DOI: 10.18434/mds2-3083 ; landing page https://data.nist.gov/od/id/mds2-3083
- File listing used: https://data.nist.gov/rmm/records/mds2-3083 ; files from `https://data.nist.gov/od/ds/mds2-3083/<url-encoded path>`
- Related papers: 10.1016/j.msea.2019.06.003 (Paper1), 10.1016/j.addma.2020.101425 (Paper2), 10.1016/j.actamat.2021.117464 (Paper3)
- Date accessed: 2026-11-25
- Licence (as stated in the record metadata): https://www.nist.gov/open/license

## What was downloaded

Only small csv/xlsx/txt/m files (66 files, about 1.6 MB) of a 932-file, 48.3 GB record. Folder layout is kept.

## Content (checked by opening the tables)

- Alloy IN718, laser powder bed fusion, as-built. The machine is not named in the downloaded files (see the papers).
- `Paper1/Surface Roughness/sr.csv`: 1 header row + 1 units row + 269 data rows, 22 columns. Per roughness scan:
  Label (P1..P12 only, although the readme speaks of 25 parameter sets), Sample, Side, Scan, Ra, Rms (um), three columns whose
  headers are garbled (`?`, `?`, `?s`; meaning unclear), contour power/speed (115 W, 625 mm/s, constant), fill power
  (115-465 W, 7 values), fill scan speed (620-1770 mm/s, 7 values), layer thickness (30 or 60 um), build orientation
  (file values 0/45/90; the readme says 0/60/90, unresolved), spot size 50 um, fatigue stress (330 MPa mean, R 0.1) and fatigue life (cycles).
  It contains NO porosity column.
- `Paper2/Micro-CT data/Data/P{1,2,4}_{LOF,gas}.csv`: per-pore lists from X-ray CT for three build conditions only
  (per `Paper2/paper2_readme.txt`: P1 = 60 deg, 62 J/mm3; P2 = 0 deg, 62 J/mm3; P4 = 0 deg, 38 J/mm3, 6.91 % porosity by volume).
  No header row; 3 numeric columns. Rows: P1_LOF 57, P1_gas 580, P2_LOF 120, P2_gas 532, P4_LOF 196, P4_gas 546.
  Column meaning is NOT documented in the record. Column 3 equals the equivalent-sphere diameter of column 2 to within 0.2 %
  (so column 2 is probably a pore volume in voxels and column 3 a diameter in voxel units; column 1 probably a label or slice index). This is an inference.
- `Paper2/Tensile data/*.txt`, `ANOVA/*.TXT`, specimen list xlsx: tensile properties per specimen for P1/P2/P4.
- `Paper3/Meso-Load-Strain/*`: load-strain curves of six meso tensile specimens (P1, P2). `Paper3/Simulated Dataset/pores/{fracture,vdf_values}.csv`:
  simulated void-descriptor-function data (120 synthetic specimens per the readme); no headers.
- `ReadMe_ScanDetails*.txt`: XCT scan settings.

## Caveat for defect benchmarking

Porosity versus process parameters exists for three conditions only (two energy densities). The pore-size tables are per-pore lists, not porosity fractions.
The power/speed/layer of P1/P2/P4 in the Paper2 naming must be taken from the paper; do not assume the `Label` P1..P12 of `sr.csv` is the same numbering (not verified).

## Not downloaded (listed in the record)

| folder | files | bytes | content |
|---|---:|---:|---|
| Paper2/Micro-CT data/P1, P2, P4 | 12 | 7 019 388 007 | XCT tif stacks (raw and processed) |
| Paper3/Micro-CT data/P1, P2, P4 | 61 (minus the readme txt files) | 11 885 804 769 | XCT tif stacks, pre and post test (13 tif_stack zip files of 89-250 MB among them) |
| Paper2/EBSD data | 130 | 28 298 385 057 | .osc/.oim/.ang EBSD, maps, EBSD_ReadMe.docx, reference directions pdf/pptx |
| Paper2/Fractography, SEM, EDS | about 300 | about 0.6 GB | SEM images, EDS line scans |
| Paper1/SEM/Powder, Surface Roughness/Images and Tests | 363 | about 0.38 GB | powder SEM, profilometry images and .dat height maps |
| Paper3/IN718meso.inp, Simulated Dataset zip files, Optical Images | - | about 0.11 GB | FE input deck, zipped simulation inputs/outputs, optical images |
| Paper2/Tensile data/raw data.zip | 1 | 4 921 344 | raw tensile data (zip not downloaded) |

## Files (downloaded 2026-11-25)

| file | bytes | sha256 |
|---|---:|---|
| `readme.txt` | 1908 | `6275f2e52856e8e8be6faaeb8165f929c0ca455c8fe4d2fd53961edd89ef6f5b` |
| `Paper1/paper1_readme.txt` | 2726 | `bff215b73379275d619e202aa7790eff46b67af275ce7022a76545ed3fefba6c` |
| `Paper1/Surface Roughness/ReadMe.txt` | 405 | `e1d729ed529623c5b905da9d1b6758a3f229f24e70dbb582b14fa06ea97f09a5` |
| `Paper1/Surface Roughness/sr.csv` | 27548 | `733aa044b9e83b245916a7b970899f4fa31bea80ad7b22b11555ba005225e844` |
| `Paper2/paper2_readme.txt` | 3972 | `fe6bd85fa12e51f8bb194e86b17351981eb1a8f9de304346fbb9791ca5a23677` |
| `Paper2/Micro-CT data/Data/P1_LOF.csv` | 1165 | `f2ab64a2a2b0337b3404325619a39e0de6cb8467b8767d59fdea70f5e5fd0070` |
| `Paper2/Micro-CT data/Data/P1_gas.csv` | 8430 | `0378d4b1c41596b5ba180634054ee2c57bdd62174c68f12c24364920d90549ff` |
| `Paper2/Micro-CT data/Data/P2_LOF.csv` | 2519 | `187582804d03b03a77bafb762604828918e2bd0abf62a641ffa5a1d16b5a1f51` |
| `Paper2/Micro-CT data/Data/P2_gas.csv` | 10842 | `1aacb305422116856481195c3f7479726b88fe4e0ffd0b6158fae9fced9afb54` |
| `Paper2/Micro-CT data/Data/P4_LOF.csv` | 4123 | `36fad923af94d445dd6366e8e935663680fe123a542ad61656d78fed33bb1fec` |
| `Paper2/Micro-CT data/Data/P4_gas.csv` | 11212 | `2c089ce8c2b8e73a5e828a09e288797ba233c5800d14834b147537fceb09b6e3` |
| `Paper2/Micro-CT data/P1/ReadMe_ScanDetails.txt` | 821 | `4a2b3d571c2dd76e070b4d345c304c0e1b86221b44fcbdac4ea15db2c2e2ac24` |
| `Paper2/Micro-CT data/P2/ReadMe_ScanDetails.txt` | 821 | `1a86f5c07c2b8d28d4f1fab3ba76bbea80ac1ae17157f0e5e41b1537b622e3da` |
| `Paper2/Micro-CT data/P4/ReadMe_ScanDetails.txt` | 821 | `c5d179d02a9f5825b6bc8e9275b4fd4885b99cfd6b62a84281ac77c9ef894136` |
| `Paper2/Tensile data/elong.txt` | 691 | `6c2b953368e702ae63124c41ef90966af681569c9d0670339f0a491ff68b3fdc` |
| `Paper2/Tensile data/modulus.txt` | 843 | `eba4382d7778d982588e147685981b8f76eda2778928d924f7bc3c1bf640db9e` |
| `Paper2/Tensile data/specimen list for mini tensile testing of LM IN718.xlsx` | 26348 | `08e13213107f232dbbe508293535d40931e5fe4904423cd5d73b72ae57497c2a` |
| `Paper2/Tensile data/strainultimate.txt` | 695 | `227b1cfb13fb44be9f1e724341a4766a1c6b8df3dc57e907a281edd8b0c06fe9` |
| `Paper2/Tensile data/tensile.txt` | 618 | `daafc273d52883831d6700187a24b6ef3997ba7b8ab30d4cccc1d76c4d269379` |
| `Paper2/Tensile data/ultimate.txt` | 665 | `88d0a9d6e32d5c40dbf557e2fd3b5411a9055c63050c12885845612c16896ede` |
| `Paper2/Tensile data/ANOVA/stress.TXT` | 3492 | `54272db0bf3d92970024d6c0336b35de4e157359c586c6c56257483d9da1c0c1` |
| `Paper2/Tensile data/ANOVA/total_e.TXT` | 3505 | `c8802248b3685500f1255d654c0bf266260b6cea157cf4f64dc11dad76a5adb3` |
| `Paper2/Tensile data/ANOVA/uniform_e.TXT` | 3506 | `db60873af5644240e432685b5d70de751450de222339e7c09520d3fdf1d89542` |
| `Paper2/Tensile data/ANOVA/uts.TXT` | 3347 | `87e63086f57b1eb005f56b6ea43b7f27a06ac57d14dbf12e80023daa0d700d31` |
| `Paper2/Tensile data/ANOVA/youngs.TXT` | 3518 | `e33bbdf6ba8f7dd0feb54ff19d351bc66572877a82032e7fda74c6a3c6899a29` |
| `Paper3/paper3_readme.txt` | 4894 | `d2d47326cfe863d330adfb92d7c804b755dc8e41af00371aded40e790f412b05` |
| `Paper3/Meso-Load-Strain/IN718_P1-meso1_load-strain.xlsx` | 45678 | `772180e73c4d6fe3a247e805cc37e98a5eb94eb107e7aa08a5be989f57223a2a` |
| `Paper3/Meso-Load-Strain/IN718_P1-meso1_load-strain_actual_areas.csv` | 6238 | `6a0430a05b0367f9b36d951baf86d7d47362b2ea89d6dbdc99e5220b7a8618b8` |
| `Paper3/Meso-Load-Strain/IN718_P1-meso2_load_strain.xlsx` | 61819 | `7f6aaec7867dd300919fc0e393267de3d1c34d069bf932c78a2712e4fbc39930` |
| `Paper3/Meso-Load-Strain/IN718_P1-meso2_with_interruption.csv` | 4608 | `e062ce8090fc83fba1982475fec5ba3bc663c6f58df488075fc1159e519ce0c2` |
| `Paper3/Meso-Load-Strain/IN718_P1-meso2_with_interruption_actual_areas.csv` | 5626 | `eee3195eb7355b91f7e56f077df5955fa1fc28937f059cb04fda3d256040278a` |
| `Paper3/Meso-Load-Strain/IN718_P1-meso3_load-strain.xlsx` | 88360 | `c67a29fbf5fa33fdea5c1dcf84494169ea2a1f0cb9f1d64b9b4174d0bf034265` |
| `Paper3/Meso-Load-Strain/IN718_P1-meso3_load-strain_actual_areas.csv` | 7881 | `de14be33498b7bc8e3b8577e72d7bbc3de761bfba978c1c9800eb3803cdf1fc4` |
| `Paper3/Meso-Load-Strain/IN718_P2-meso1_load-strain.xlsx` | 62813 | `c2349bfccee0b436689d7e153a0bbc884ac7c86428162265ba6d0b03bc868a40` |
| `Paper3/Meso-Load-Strain/IN718_P2-meso1_load-strain_actual_areas.csv` | 9785 | `1d778b659df174e22a9d02e1d29af7be7c5bc8752476577ce857ae99699fb3ba` |
| `Paper3/Meso-Load-Strain/IN718_P2-meso2_load-strain.xlsx` | 57456 | `6d2770b7d5ca5d0a4d616a7cd3f31e170d5ae55077e1065d7ef7110ba99e6b99` |
| `Paper3/Meso-Load-Strain/IN718_P2-meso2_load-strain_actual_areas.csv` | 8937 | `4ac8362033c442f90bf8976ac15a9b29d568b7cbff64035033e35190bc01137e` |
| `Paper3/Meso-Load-Strain/IN718_P2-meso3_load-strain.xlsx` | 62988 | `6c9ab7c933f0ba4c9a12c52186177d4add07ec17b2cae2caaba23bc62592bc1e` |
| `Paper3/Meso-Load-Strain/IN718_P2-meso3_load-strain_actual_areas.csv` | 11369 | `401b211c8142c7620d90235e7be244e242d2cfb77d1ff57b92aa46af150e9378` |
| `Paper3/Meso-Load-Strain/Summary-combined-stress-strain-curves.csv` | 19583 | `89c1f0641cc543c152eee58e7cffe345449ee5b09f7b749d2d14cf183666007e` |
| `Paper3/Meso-Load-Strain/Summary-combined-stress-strain-curves.xlsx` | 106122 | `1a8a45021b4b52440a4f9ed5d418231d5e70bb9bce59b42c8755b499e2b07d65` |
| `Paper3/Micro-CT data/P1/Sample1/Post/Side1-no-scratch-side/ReadMe_ScanDetails_s1_CT_post.txt` | 813 | `a20616f10b6b21123849c783e7fbff184755731d45caab85aff854c157d2d13c` |
| `Paper3/Micro-CT data/P1/Sample1/Post/Side2-scratch-side/ReadMe_ScanDetails_s2_CT_post.txt` | 817 | `1730a0ff262bc25df97596b8fcd681353c167cb4624bfd2e8dc221cfbfc25eea` |
| `Paper3/Micro-CT data/P1/Sample1/Pre/ReadMe_ScanDetails_p1_meso1_CT_pre.txt` | 816 | `281afe722415b544a75479d1cf966f7a4aca7623beca808cd2e4934cf460d62f` |
| `Paper3/Micro-CT data/P1/Sample2/Post/Side1-scratch-side/ReadMe_ScanDetails_P1_m2_side1_CT_post.txt` | 813 | `3ea89696feb2fc2c86513ef1562462a38ebae4bd9b170dd826bc264d755ef721` |
| `Paper3/Micro-CT data/P1/Sample2/Post/Side2-no-scratch-side/ReadMe_ScanDetails_P1_m2_side2_CT_post.txt` | 812 | `a52e79c44553a506cd91fcfdc0a80f44437c76c8bb09efa566f90b8baa6e01af` |
| `Paper3/Micro-CT data/P1/Sample2/Pre/ReadMe_ScanDetails_p1_meso2_CT_pre.txt` | 815 | `24f8d0740b18dc07406fba353cfb34a8040d570574f1af57e05ada766a2ed8e7` |
| `Paper3/Micro-CT data/P1/Sample3/Post/Side1-scratch-side/ReadMe_ScanDetails_s1_CT_post.txt` | 813 | `ca96462618d9ccb9deeffa818fc92d9a847471f092ae54ba6467f5fa902df43d` |
| `Paper3/Micro-CT data/P1/Sample3/Post/Side2-no-scratch-side/ReadMe_ScanDetails_s2_CT_post.txt` | 813 | `b04a7355a58d0b2cdebd1f6d476ee3856bb2b5588866ced71bbe25b00d603bcc` |
| `Paper3/Micro-CT data/P1/Sample3/Pre/ReadMe_ScanDetails_p1_meso3_CT_pre.txt` | 815 | `997a3051bb06a1ed3310cbe7642683d31a390bc0365a69b7547e7dd637cabc45` |
| `Paper3/Micro-CT data/P2/Sample1/Post/Side1-scratch-side/ReadMe_ScanDetails_P2_m1_side1_CT_post.txt` | 813 | `54dacb461ecbec8d8d638b2e866faac3d4eaf4853cae7ad5fff9777f552239c9` |
| `Paper3/Micro-CT data/P2/Sample1/Post/Side2-no-scratch-side/ReadMe_ScanDetails_P2_m1_side2_CT_post.txt` | 808 | `8a43420e247f04544aaebad169dfa3cdc6da406c0d448fb1ebd49872d3095c4c` |
| `Paper3/Micro-CT data/P2/Sample1/Pre/ReadMe_ScanDetails_P2_m1_pre.txt` | 805 | `6959fac19ea953faf52d838cf2e34ce7040ed8204c36a60cdcb96e4d5b24a6b5` |
| `Paper3/Micro-CT data/P2/Sample2/Post/Side1-no-scratch-side/ReadMe_ScanDetails_P2_m2_side1_CT_post.txt` | 813 | `90c3f5f380a16ae85acad58dfe1830108ac792b58140c5b8845644e16444d967` |
| `Paper3/Micro-CT data/P2/Sample2/Post/Side2-scratch-side/ReadMe_ScanDetails_P2_m2_side2_CT_post.txt` | 806 | `086c39c242ce00345814dbd67afe00f16915b15c7b36bb6a5db789e6d3d1596b` |
| `Paper3/Micro-CT data/P2/Sample2/Pre/ReadMe_ScanDetails_P2_m2_pre.txt` | 807 | `775574ad62130cb295e2c331b17ee3a401e36d9741d47edc8bbbaa7c768f142d` |
| `Paper3/Micro-CT data/P2/Sample3/Post/Side1/ReadMe_ScanDetails_P2_m3_side1_CT_post.txt` | 813 | `28e64358d46a889af2cd1c3ef30a665866a2269cf8f8d9e2e031e09a2b8ac00c` |
| `Paper3/Micro-CT data/P2/Sample3/Post/Side2/ReadMe_ScanDetails_P2_m3_side1_CT_post.txt` | 814 | `25a2d0ed57a8a06aa0f03744fecd082d881087dedacebec401ea93eff196ef7f` |
| `Paper3/Micro-CT data/P2/Sample3/Pre/ReadMe_ScanDetails_P2_m3_pre.txt` | 806 | `d8ca8c279eb9250e9a82fd3e3341056c3f27b1c6a42449555a23a9950c1101d9` |
| `Paper3/Micro-CT data/P4/Sample1/pre/ReadMe_ScanDetails_P4_m1_pre.txt` | 806 | `647cbed73fa7d61e4796eeab03e9a33146b4ff3a0a44510eb476cffeb33dfdc6` |
| `Paper3/Micro-CT data/P4/Sample2/pre/ReadMe_ScanDetails_P4_m2_pre.txt` | 804 | `e5235208674a24b94fd51eb63bc6ad6d552e2a8bac92ecbcb420557793b07707` |
| `Paper3/Micro-CT data/P4/Sample3/pre/ReadMe_ScanDetails_P4_m3_pre.txt` | 806 | `25cf7dd0828c1359d9800080a7a4580b9845799f01733ceda7160031f0856f2b` |
| `Paper3/Simulated Dataset/pores/fracture.csv` | 1392 | `fec33030ee0f3797041805a83a1207306c5ef9464315b8c9f18f27ef4e6146eb` |
| `Paper3/Simulated Dataset/pores/plot_vdfs_simulated_data.m` | 1766 | `debe664025aae3fb623db037131bf43692a2fea3e34476e8330b52dd9b1359a9` |
| `Paper3/Simulated Dataset/pores/vdf_values.csv` | 603120 | `1494cd0057044915b7516629f0a110270a983378411bd68298a4474a3d7ecdc7` |
| `Paper3/Simulated Dataset/tensile/plot_all_stress_strains.m` | 333 | `5ec72f99cfa85dfdf50b1cdd1d534cebac93df6afe43475f53d32cfb589f647d` |
