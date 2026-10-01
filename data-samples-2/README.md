# Sample dataset

A small example dataset that exercises the whole **NPLab** workflow without
needing real experimental data. It is the folder pre-filled in the app sidebar,
so `streamlit run app.py` works out of the box.

Every file is a `.txt` emission scan in the instrument's native export format:
a short text header (`Labels`, `Type`, `Start`, `Stop`, `Step`, …) followed by
`wavelength,counts` rows. All spectra span **495–525 nm** in 1 nm steps, which
covers the 508 nm fluorescence peak used by default.

## Contents

| Prefix | Files | Tab | What it covers |
|--------|------:|-----|----------------|
| `BLK`  | 30 | Blank | 10 blank series (`BLK01`–`BLK10`) × 3 replicas — enough to pass the ≥ 10-file LOD/LOQ requirement |
| `CCN`  | 18 | Calibration | 6 concentrations (10, 20, 25, 50, 100, 150 mg/L) × 3 replicas, matrix `M00` |
| `DRF`  | 24 | Drift | Type-1 (`N00`) and Type-2 (`N50`), 4 session time-points × 3 replicas each |
| `CTRL` |  8 | Control | Type-1 only (`N00`, matrix `M05`): adsorbents `BTO` and `MXB` at 2 g/L, conditions `C01` and `C02`, 2 replicas |
| `AKN`  | 28 | Kinetics | N = 150 mg/L on `MXB` 1 g/L, 8 contact times (10–420 min), replicas per time-point |
| `AIN`  | 20 | Isotherms | 6 initial concentrations (10, 25, 50, 100, 150, 200 mg/L) on `MXB` 1 g/L, replicas per concentration |

## Suggested run-through

1. **Calibration** — load the `CCN*` files, keep the peak wavelength at 508 nm,
   pick degree 1, and click **Set calibration**. Every other tab unlocks.
2. **Blank** — load the `BLK*` files to obtain LOD and LOQ.
3. **Drift** — the `DRFN00*` files give the Type-1 trend (with the LOD/LOQ lines
   carried over from step 2); `DRFN50*` gives the Type-2 recovery.
4. **Control** — the `CTRL*` files show the Type-1 bar chart comparing the two
   adsorbents across conditions `C01` and `C02`.
5. **Kinetics** — the `AKN*` files fit pseudo-1st and pseudo-2nd order models.
6. **Isotherms** — the `AIN*` files fit Langmuir and Freundlich.
7. **Download Results** — export the TXT report and the per-plot CSVs.

## Not covered here

- **Type-2 control** files (`N > 0`) and the reserved **`C000`** Co-correction
  sample, so the Co correction stays greyed out in the Kinetics and Isotherms
  tabs with this dataset.
- The **`M07`** odd/even replica-parity mode, which is triggered by a sub-folder
  named `M07` in the data path.
- The `.asc` / `.akn` headerless formats — all samples here are `.txt`.

See the [root README](../README.md) for the full file naming convention.
