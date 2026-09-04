"""
tabs/help.py — Help tab.

Self-contained user guide: workflow, models, troubleshooting and the full
file naming convention that the parsers in utils.py implement.
"""

import streamlit as st


def render() -> None:

    st.subheader("📖 Help & Documentation")

    st.markdown("""
This application processes fluorescence emission spectra to characterise
nanoplastic adsorption. It guides you through instrument quality control,
calibration, and then the actual adsorption experiments.

---

### 📥 Loading your data

Data are loaded from the **sidebar**, in one of two modes:

- **📁 Folder** — type a path to a folder on the machine running the app. A
  selector lets you drill down into its sub-folders (two levels deep).
- **🗜️ ZIP file** — upload a single ZIP containing all your spectral files.
  Sub-folders inside the ZIP are flattened automatically.

Independently of the sidebar, **every analysis tab also has its own
`📁 Folder` / `⬆️ Upload` switch**, so you can drag files straight into a
single tab without changing the global data source.

Three spectrum formats are recognised — `.txt` (text header, then
`wavelength,counts` rows), `.asc` and `.akn` (headerless, four comma-separated
integer fields per line, where the comma doubles as the decimal separator).

---

### 🗺️ Tab overview

| Tab | Purpose | File prefix |
|-----|---------|-------------|
| **Visualization** | Inspect any spectral file visually | any |
| **Calibration** | Build a counts-vs-concentration curve | `CCN`, `HCCN`, `LCCN` |
| **Blank** | Compute LOD / LOQ from blank spectra | `BLK` |
| **Drift** | Monitor laser intensity drift across a session | `DRF` |
| **Control** | Evaluate the effect of experimental factors | `CTRL` |
| **Kinetics** | Fit adsorption kinetics models | `AKN` |
| **Isotherms** | Fit Langmuir / Freundlich isotherm models | `AIN` |
| **Download Results** | Export the report (TXT) and all CSVs | — |
| **Help** | This page — workflow, models and naming convention | — |

---

### 🔬 Recommended workflow

#### 1 — Visualization tab
1. Point the app to your data folder (sidebar) **or** upload files directly.
2. Select any spectral files to preview their emission spectra.
3. Two plots are shown per selection — counts vs wavelength and
   log₁₀(counts) vs wavelength.
4. Use this tab to perform a quick sanity-check before committing to analysis.

#### 2 — Calibration tab
1. Select (or upload) calibration files (`CCN*`, `HCCN*`, `LCCN*`).
2. Verify that concentrations were parsed correctly from filenames in the
   summary table.
3. Set the **peak wavelength** with the slider (default **508 nm**, when the
   loaded spectra cover it). This value is reused by every other tab, so change
   it here if your fluorophore peaks elsewhere.
4. Inspect the **emission spectra** — all curves should be plausible.
5. Replicas of the same (N, M) group are averaged; choose a **polynomial
   degree** (up to *number of averaged points − 1*) and inspect the fit
   (R² and MAE shown).
6. Click **✅ Set calibration** to lock the model.
   > ⚠️ The Blank, Drift, Control, Kinetics and Isotherms tabs all use this
   > calibration to convert counts to concentration, and stay blocked until
   > it is set.

#### 3 — Blank tab
1. Select (or upload) blank files (`BLK*`).
2. **Requirements:** a linear (degree = 1) calibration must be active, and at
   least **10 BLK files** must be loaded.
3. The app checks for **negative counts** and applies an offset where needed
   (the offset per file is listed in the summary table; offset-corrected points
   appear as orange diamonds in the plot).
4. LOD and LOQ are computed as:  
   `LOD = 3 × sB / b`  and  `LOQ = 10 × sB / b`  
   where `sB` is the sample standard deviation of blank max-counts and `b`
   is the calibration slope.
5. The results are stored and reused as reference lines in the Drift plots.

#### 4 — Drift tab
1. Select (or upload) DRF files (`DRF*`).
2. Set the **wavelength window** for Type-1 median extraction (default 495–525 nm).
3. **Type-1 (N = 0 — no nanoplastics):**
   - Treated as blanks: negative-count offset applied if needed.
   - The **median** of counts within the wavelength window is extracted per
     file (more robust than mean for noisy blank spectra).
   - Medians are averaged across replicas for each session time-point
     (t1, t2, t3 …).
   - Two plots: signal (cps) vs t, and equivalent concentration (mg/L) vs t.
     Both include a dashed linear-trend line (slope and R² reported).
   - LOD / LOQ reference lines from the Blank tab are overlaid if available.
4. **Type-2 (N > 0 — with nanoplastics):**
   - Count at the calibration peak wavelength extracted per file and averaged
     across replicas per (N, t) group.
   - Two plots: signal (cps) vs t and concentration (mg/L) vs t, per N-series,
     with dotted nominal-concentration reference lines and a deviation table.
5. A linear drift trend slope close to zero indicates stable laser performance.

#### 5 — Control tab
1. Select (or upload) CTRL files (`CTRL*`).
2. Each **C-code** identifies one specific experimental condition (e.g. C01 =
   4 000 rpm centrifugation, C02 = 10 000 rpm). Add optional text descriptions
   per code in the expandable label widget — these appear in chart axes.
3. Set the **wavelength window** used for the Type-1 median (default 495–525 nm).
4. **Type-1 (N = 0 — adsorbent present, no nanoplastics):**
   - Blank-like treatment: offset correction, then median within the wavelength
     window.
   - Mean of medians per (adsorbent series × C-code) shown as a **grouped bar
     chart** (signal in cps).
5. **Type-2 (N > 0 — nanoplastics present, no adsorbent):**
   - Count at the calibration peak wavelength averaged across replicas per
     (N × C-code) group.
   - Converted to concentration (mg/L) using the calibration curve.
   - Grouped bar chart of recovered concentration vs C-code, with dotted
     nominal-N reference lines and a deviation-from-nominal table.
6. When both types are present a **cross-type summary table** is shown.
7. **Co correction (reserved code `C000`)** — see the dedicated section below.

#### 6 — Kinetics tab
1. Select (or upload) kinetics files (`AKN*`).
2. Optionally enable the **Co correction** (see below) before reading the results.
3. The app:
   - Reads counts at the calibration peak wavelength from each spectrum.
   - Inverts the calibration curve to obtain Cₜ.
   - Computes q(t) = (C₀ − Cₜ) / C_ads, where C₀ and C_ads are parsed from
     the filenames.
   - Averages replicas (R1, R2, …). If the selected data **folder path**
     contains a sub-folder named `M07`, replicates are instead split into
     **odd** (R1, R3, …) and **even** (R2, R4, …) sub-series.
4. Pseudo-1st and pseudo-2nd order models are fitted per series.
5. Linearised diagnostic plots are shown alongside the main q(t) chart.

#### 7 — Isotherms tab
1. Select (or upload) isotherm files (`AIN*`).
2. Optionally enable the **Co correction** (see below) before reading the results.
3. The app:
   - Reads counts at the calibration peak wavelength from each spectrum.
   - Uses the calibration curve to compute **Cₑ** (equilibrium concentration).
   - Computes qₑ = (N − Cₑ) / C_ads and **%R = (N − Cₑ) / N × 100** (percentage
     removal), where N and C_ads are extracted from the filename.
   - Groups and averages replicates by adsorbent, adsorbent concentration, and pH
     (when multiple pH values are present, grouping includes pH automatically).
4. Fits **Langmuir** and **Freundlich** models for each series. Additional models
   (Sips, Tóth, Halsey, Harkins-Jura, Janovics) can be selected on demand.
5. Displays nonlinear fitted curves (qₑ vs Cₑ) and linearised diagnostic plots.

#### 8 — Download Results tab
1. Review the status tiles (one per tab) to confirm which analyses are complete.
2. Click **⬇️ Download results as TXT** for a plain-text summary of all
   computed parameters.
3. Download individual **CSV files** for every plot and table produced across
   all tabs.

---

### 🧪 Co correction (optional)

Agitation and centrifugation can change the *actual* initial nanoplastic
concentration with respect to the nominal value written in the filename. To
quantify that, measure a control sample with the **reserved condition code
`C000`** (e.g. `CTRLN100M10pH4ADS0C000R1.txt`) and process it in the
**Control** tab. The measured-vs-nominal table it produces then unlocks an
**"Apply Co correction"** checkbox in the **Kinetics** and **Isotherms** tabs.

When enabled, the measured value replaces the nominal N in the qₑ, %R and q(t)
calculations only — the calibration curve itself is never modified. Matching
rules, in order:

1. Exact (pH, N) match → its measured value is used directly.
2. Same pH but no exact N → N is scaled by the mean measured/nominal ratio at
   that pH.
3. No rows for that pH → the whole table is used, ignoring pH.

Whether the correction was applied is recorded in the exported TXT report.

---

### 📐 Kinetic models

**Pseudo-1st-order**:
$$q(t) = q_e \\left(1 - e^{-k_1 t}\\right)$$

**Pseudo-2nd-order**:
$$q(t) = \\frac{k_2\\,q_e^2\\,t}{1 + k_2\\,q_e\\,t}$$

| Parameter | Units | Description |
|-----------|-------|-------------|
| qₑ | mg/g | Equilibrium adsorption capacity |
| k₁ | min⁻¹ | Pseudo-1st-order rate constant |
| k₂ | g·mg⁻¹·min⁻¹ | Pseudo-2nd-order rate constant |

---

### 📐 Isotherm models

**Langmuir**:
$$q_e = \\frac{q_{max} K_L C_e}{1 + K_L C_e}$$

Separation factor:
$$R_L = \\frac{1}{1 + K_L C_0}$$

**Freundlich**:
$$q_e = K_F C_e^{1/n}$$

| Parameter | Units | Description |
|-----------|-------|-------------|
| q_max | mg/g | Maximum adsorption capacity |
| K_L | L/mg | Langmuir constant |
| R_L | — | Separation factor (adsorption favorability) |
| K_F | (mg/g)(L/mg)^(1/n) | Freundlich constant |
| n | — | Adsorption intensity |
| %R | % | Percentage removal = (N − Cₑ) / N × 100 |

**Interpretation:**
- Langmuir:
  - 0 < R_L < 1 → favorable
  - R_L = 1 → linear
  - R_L > 1 → unfavorable
  - R_L → 0 → irreversible

- Freundlich:
  - 1/n < 1 → favorable
  - 1/n = 1 → linear
  - 1/n > 1 → unfavorable

---

### ❓ Troubleshooting

| Symptom | Likely cause |
|---------|-------------|
| Files not listed in the folder selector | Wrong folder path in the sidebar, or filenames lacking the expected prefix |
| Concentration parsed as "–" | Filename does not match the expected pattern (see naming convention below) |
| "Fit did not converge" | Too few time points, or data are noisy / inconsistent |
| Counts at the peak wavelength returns None | Spectrum does not cover the selected wavelength |
| Polynomial degree slider missing | Only two averaged calibration points — a linear fit is forced |
| LOD/LOQ section is disabled | Calibration degree > 1, or fewer than 10 BLK files loaded |
| Negative counts in blanks / DRF / CTRL | Normal for some instruments — the app applies an automatic offset and reports it |
| Any analysis tab blocked | Calibration not set yet — go to the Calibration tab first |
| "Apply Co correction" greyed out | No `C000` control file has been processed in the Control tab |

---

### 🗂 File naming convention

This is the complete, authoritative reference — it describes exactly what the
parsers in `utils.py` accept.

**General rules**

- Use **alphanumeric characters only**, in one continuous sequence. No spaces,
  hyphens, underscores, dots, slashes, commas or parentheses.
- Blocks appear in a **fixed order** and are recognised by their letter prefix,
  not by a fixed width — so `N5`, `N50` and `N150` are all valid.
- Every file ends with the **replica** block `R{n}` (`R1`, `R2`, …). Use `R1`
  even when there is only one measurement.
- **Decimals use `p` instead of a dot**: `10p5` → 10.5, `N0p5` → 0.5 mg/L.
  There is no leading-zero convention — write `N0p5`, *not* `N050`.
- An optional **`pH{value}`** token goes **immediately after the matrix code**
  (`M05pH4` → matrix M05, pH 4). It is only needed when more than one pH is
  measured in the same session.
- Matching is **case-insensitive**.
- The **matrix** block `M{XX}` is required for calibration, drift and control
  files, and optional for kinetics and isotherm files.
- A parameter that is constant across the whole session can be left out of the
  filename entirely.

> **Note on `CC` / `HCC` / `LCC`:** the app accepts all three prefixes but does
> **not** treat them differently — calibration points are grouped by
> (concentration, matrix) alone. Keep High and Low curves in separate folders,
> or give them distinct matrix codes, if you do not want them averaged together.

#### Calibration files — `[H|L]CCN{conc}M{XX}R{n}.txt`
| Token | Meaning | Example |
|-------|---------|---------|
| `CC` / `HCC` / `LCC` | Prefix (optional H = high range, L = low range) | — |
| `N{conc}` | Analyte concentration (mg/L) | `N5`, `N50` |
| `M{XX}` | Matrix / run code | `M00` |
| `R{n}` | Replica number | `R1` |

Examples: `CCN5M00R1.txt` → 5 mg/L, replica 1;  `HCCN50M00R1.txt` → 50 mg/L (high range)

#### Blank files — `BLK…R{n}.txt`
Any filename starting with `BLK` is treated as a blank.  
Example: `BLKM00R1.txt`, `BLK_run1_R3.asc`

#### Drift files — `DRFN{N}M{M}t{t}R{n}.txt`
| Token | Meaning | Example |
|-------|---------|---------|
| `N{N}` | NP concentration (mg/L); 0 = Type-1 (no NPs) | `N00`, `N50` |
| `M{M}` | Matrix code | `M00` |
| `t{t}` | Session time-point index | `t1`, `t2`, `t3` |
| `R{n}` | Replica number | `R1` |

Examples: `DRFN00M00t1R1.txt` (Type-1),  `DRFN50M00t2R1.txt` (Type-2)

#### Control files — `CTRLN{N}M{M}{AdsCode}{AdsConc}C{code}R{n}.txt`
| Token | Meaning | Example |
|-------|---------|---------|
| `N{N}` | NP concentration; 0 = Type-1 (no NPs) | `N00`, `N50` |
| `M{M}` | Matrix code | `M00` |
| `{AdsCode}` | Adsorbent identifier (letters only) | `BTO`, `MXB` |
| `{AdsConc}` | Adsorbent concentration (g/L); 0 = no adsorbent | `3`, `0` |
| `C{code}` | Condition code (user-defined; `C000` is reserved) | `C01`, `C02` |
| `R{n}` | Replica number | `R1` |

Examples: `CTRLN00M00BTO3C01R1.txt` (Type-1),  `CTRLN50M00ADS0C02R1.txt` (Type-2),
`CTRLN100M10pH4ADS0C000R1.txt` (Co correction)

#### Kinetics files — `AKN{N}M{M}{AdsCode}{AdsConc}t{time}R{n}.txt`
| Token | Meaning | Example |
|-------|---------|---------|
| `N{N}` | Initial analyte concentration (mg/L) | `N50`, `N150` |
| `M{M}` | Matrix code — **optional** | `M00` |
| `{AdsCode}` | Adsorbent material code (letters only) | `BTO`, `MXB` |
| `{AdsConc}` | Adsorbent concentration (g/L) | `2`, `1` |
| `t{time}` | Contact time (see note below) | `t05`, `t120`, `t24p0` |
| `R{n}` | Replica number | `R1` |

Examples: `AKN50BTO2t05R1.txt` (no matrix token),
`AKN150M00MXB1t420R1.txt` → N=150 mg/L, matrix M00, MXB 1 g/L, t=420 min, replica 1

> **Note on time tokens:** plain integers are minutes (`t120` = 120 min);
> p-decimal values are hours (`t24p0` = 24 h = 1440 min, `t0p25` = 15 min).

#### Isotherm files — `AIN{N}M{M}{AdsCode}{AdsConc}R{n}.txt`
| Token | Meaning | Example |
|-------|---------|---------|
| `N{N}` | Initial concentration (mg/L) | `N100`, `N200` |
| `M{M}` | Matrix code — **optional** | `M00` |
| `{AdsCode}` | Adsorbent material code (letters only) | `MXB`, `BTO` |
| `{AdsConc}` | Adsorbent concentration (g/L) | `1` |
| `R{n}` | Replica number | `R1` |

Examples: `AIN100MXA1R2.txt` → N = 100 mg/L, MXA 1 g/L, replica 2;
`AIN50M05pH4BTO10R1.txt` → N = 50 mg/L, matrix M05, pH 4, BTO 10 g/L, replica 1

---
""")
