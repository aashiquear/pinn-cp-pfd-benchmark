# MeltpoolNet (Akbari et al. 2022): melt-pool mode / defect classification and geometry compilation

- Paper: P. Akbari, F. Ogoke, N.-Y. Kao, K. Meidani, C.-Y. Yeh, W. Lee, A. Barati Farimani, "MeltpoolNet: Melt pool characteristic prediction in Metal
  Additive Manufacturing using machine learning", Additive Manufacturing 55 (2022) 102817, DOI 10.1016/j.addma.2022.102817
- Date accessed: 2026-11-25
- Source A: https://github.com/BaratiLab/MeltpoolNet (branch main, commit 6d68e2acaad8d074134a7f3d843278649774d661, 2022-07-13), files
  `https://raw.githubusercontent.com/BaratiLab/MeltpoolNet/main/{Data/meltpoolnet_classification.csv,Data/meltpoolnet_regression.csv,README.md}`.
  Licence: NONE stated (no LICENSE file; the README only asks for citation). Reuse terms are therefore unclear.
- Source B: Hugging Face dataset `ppak10/Additive-Manufacturing-Benchmark`, revision c98eefe7bb5315023b38f2a7f79bb0f25fdba16c, files
  `https://huggingface.co/datasets/ppak10/Additive-Manufacturing-Benchmark/resolve/<revision>/{source/meltpoolnet/meltpoolclassification.csv,source/meltpoolnet/meltpoolgeometry.csv,source/meltpoolnet/citation.bib,README.md}`.
  Licence: MIT (dataset card). This is a third-party re-packaging and a SUBSET of the GitHub tables.

## Content (checked by opening the csv files)

Literature compilation: many alloys, many machines (machine not recorded per row). Units are not written in the GitHub files; ranges indicate W, mm/s and um.

- A `meltpoolnet_classification.csv`: 1249 rows, 50 columns (5 unnamed trailing columns with stray notes). Material (21; SS316L 305, Ti-6Al-4V 228, SS17-4PH 228, IN718 169,
  IN738LC 47, IN625 24, ...), Process (PBF), Sub-process (SLM 1221, EBM 28), Power 40-2700, Velocity 0.6-10000, Hatch spacing (232 filled), layer thickness (762 filled),
  beam D (1097 filled), material properties, derived ratios (contain `#DIV/0!`), composition wt%, paper ID.
  Label `meltpool shape`: keyhole 439, desirable 364, LOF 305, balling 133, spatter formation 8. `relative density` filled in 28 rows and `porosity` in 4 rows only.
- A `meltpoolnet_regression.csv`: 1796 rows, 54 columns. Same parameters plus depth (1443 filled), width (1148), length (322) of the melt pool; PBF 1688 and DED 108 rows;
  SS316L 639, Ti-6Al-4V 458, IN718 147, IN625 97, IN738LC 50, Hastelloy X 41; `meltpool shape` filled in 786 rows; `relative density` 95 rows; `porosity` 16 rows.
- B `meltpoolclassification.csv`: 340 rows, 39 columns, LPBF only; SS316L 100, IN718 70 (written as `Inconel 718` 36 + `IN718` 34), Ni-5Nb 40, WE43 35, ...;
  label: desirable 112, LOF 96, keyhole 78, balling 43, empty 11.
- B `meltpoolgeometry.csv`: 464 rows, 16 columns with units in the header (power_w, velocity_mm_s, beam_diameter_micron, layer_height_micron, hatch_spacing_micron,
  measurement method, depth/width/length in micron, energy density) and a DOI per row; Ti-6Al-4V 164, SS316L 94, IN718 47, Hastelloy X 41.

## Caveats

Labels are the melt-pool mode assigned by the compilers from literature (criteria noted in a comment column, e.g. `d/t>1.5 LOF`), not measured porosity.
Tables are untidy (stray columns, `#DIV/0!`, mixed material spellings). Not downloaded: `MeltpoolNet.ipynb` (51 MB) and the parquet files of the Hugging Face dataset.

## Files (downloaded 2026-11-25)

| file | bytes | sha256 |
|---|---:|---|
| `github_BaratiLab_MeltpoolNet/README.md` | 1781 | `77f379a5ceaaea07b3fa556cd04ed9cba60c89d441f24378f670ab18fca42ccb` |
| `github_BaratiLab_MeltpoolNet/meltpoolnet_classification.csv` | 239371 | `5c2a2a719af790acd5f948f7f17c22c806b41f48283ca31e7cdf439c63a2bf80` |
| `github_BaratiLab_MeltpoolNet/meltpoolnet_regression.csv` | 362459 | `197b111bcb3c531d38c8dd7ea0d4e1f78c8b6b03e0fd900a5505c2b18414a4f6` |
| `huggingface_ppak10_Additive-Manufacturing-Benchmark/README.md` | 6109 | `4c915d775250b481535c8d572eee6f5bd393301342d223bed29357b8074fb2de` |
| `huggingface_ppak10_Additive-Manufacturing-Benchmark/citation.bib` | 1918 | `7bf68f188caa1ef1fd72c324184ed26970072af792b36aeebaa5feff8202b55e` |
| `huggingface_ppak10_Additive-Manufacturing-Benchmark/meltpoolclassification.csv` | 47276 | `5376751cb1e92dad63c7dce859afbec25318f6bcbfadf2125737fac136d1752b` |
| `huggingface_ppak10_Additive-Manufacturing-Benchmark/meltpoolgeometry.csv` | 118475 | `2abe96982d4c2abb04ce912b3502bd3f7a9f80ac4a7eddaf04b6869663a13ccd` |
