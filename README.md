# NPLab — NanoPlastic Adsorption Analysis

A [Streamlit](https://streamlit.io) application for processing fluorescence
emission spectra of nanoplastic adsorption experiments. It takes raw spectra
straight from the spectrofluorometer and walks through the whole analysis:

1. **Instrument quality control** — blank, drift, and control experiments.
2. **Calibration curve** — counts at the fluorescence peak vs. known concentration.
3. **Adsorption kinetics** — pseudo-1st and pseudo-2nd order model fitting from
   time-resolved spectra.
4. **Adsorption isotherms** — Langmuir, Freundlich, and additional models.
5. **Export** — a plain-text report plus one CSV per plot and table.

All experimental metadata (concentration, matrix, pH, adsorbent, contact time,
replica) is read directly from the **filenames**, following the
[naming convention](#-file-naming-convention) described below. The same
reference is available in the app's **Help** tab.

---

## 🚀 Quick start

```bash
git clone https://github.com/FS-UPC/NPLab.git
cd NPLab

python -m venv .venv && source .venv/bin/activate   # optional but recommended
pip install -r requirements.txt

streamlit run app.py
```

The app opens at <http://localhost:8501>.

Requires **Python 3.10+** (the code uses `X | None` type syntax).

### Authentication

`app.py` gates the interface behind a simple login backed by Streamlit's
secrets file. Create `.streamlit/secrets.toml` before the first run — it is
deliberately **not** tracked by git:

```toml
[credentials.usernames.alice]
name     = "Alice Example"
password = "choose-a-password"

[cookie]
name        = "nplab_auth"
key         = "some-random-string"
expiry_days = 30
```

> ⚠️ Passwords are compared in plain text. This is a lightweight barrier for a
> lab-internal deployment, not a security boundary — do not expose the app to
> an untrusted network without putting a proper authentication layer in front
> of it.

### Sample data

The [`data-samples/`](data-samples) folder ships a small anonymised dataset
covering every file type, so the whole workflow can be exercised out of the
box. It is the default folder in the sidebar.

---

## 📥 Loading data

Data are selected from the **sidebar** in one of two modes:

- **📁 Folder** — a path on the machine running the app, with a selector to
  drill down two levels of sub-folders.
- **🗜️ ZIP file** — a single upload containing all spectral files; sub-folders
  inside the archive are flattened automatically.

Each analysis tab additionally offers its own `📁 Folder` / `⬆️ Upload` switch,
so files can be dropped into one tab without changing the global data source.

Three spectrum formats are recognised:

| Extension | Layout |
|-----------|--------|
| `.txt` | Text header, then `wavelength,counts` rows |
| `.asc` / `.akn` | Headerless; four comma-separated integer fields per line, the comma doubling as decimal separator |

---

## 🗺️ Tab overview

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
| **Help** | In-app guide, models and naming convention | — |

---

## 🔬 Workflow

### 1 — Visualization tab
1. Point the app to your data folder (sidebar) or upload a ZIP file.
2. Select any spectral files to preview their emission spectra.
3. Counts vs wavelength and log₁₀(counts) vs wavelength are plotted side by side.
4. Use this tab for a quick sanity-check before committing to analysis.

### 2 — Calibration tab
1. Select (or upload) calibration files (`CCN*`, `HCCN*`, `LCCN*`).
2. Check the file summary table; verify concentrations were parsed correctly.
3. Set the **peak wavelength** with the slider — default **508 nm** when the
   loaded spectra cover it. Every other tab reads counts at this wavelength, so
   change it here if the fluorophore peaks elsewhere.
4. Inspect the **emission spectra** — all curves should be plausible.
5. Replicas of the same (N, M) group are averaged. Choose a **polynomial
   degree** (at most *number of averaged points − 1*) and inspect the
   calibration curve fit (R² and MAE shown).
6. Click **✅ Set calibration** to lock the model.  
   ⚠️ The Blank, Drift, Control, Kinetics and Isotherms tabs stay blocked until
   this step is completed.

### 3 — Blank tab
1. Select (or upload) blank files (`BLK*`).
2. **Requirements:** a linear (degree = 1) calibration must be active, and at
   least **10 BLK files** must be loaded.
3. Negative counts are corrected with a per-file offset; the offset is listed
   in the summary table and corrected points appear as orange diamonds.
4. LOD and LOQ are computed as:  
   `LOD = 3 × sB / b`  and  `LOQ = 10 × sB / b`  
   where `sB` is the sample standard deviation of blank max-counts and `b` is the calibration slope.
5. The results are stored and reused as reference lines in the Drift plots.

### 4 — Drift tab
1. Select (or upload) DRF files (`DRF*`). Set the **wavelength window** for Type-1 median extraction (default 495–525 nm).
2. **Type-1 (N = 0 — no nanoplastics):** blank-like treatment; median of counts within the wavelength window extracted per file and averaged across replicas per session time-point; signal (cps) and equivalent concentration (mg/L) vs. t with a dashed linear-trend line; LOD/LOQ reference lines overlaid if available.
3. **Type-2 (N > 0 — with nanoplastics):** counts at the calibration peak wavelength extracted and averaged across replicas per (N, t) group; signal and concentration vs. t per N-series with nominal reference lines and a deviation table.
4. A linear drift trend slope close to zero indicates stable laser performance.

### 5 — Control tab
1. Select (or upload) CTRL files (`CTRL*`). Set the **wavelength window** used for the Type-1 median (default 495–525 nm).
2. Each **C-code** identifies one specific experimental condition (e.g. C01 = 4 000 rpm centrifugation, C02 = 10 000 rpm). Add optional text labels per code.
3. **Type-1 (N = 0 — adsorbent present, no nanoplastics):** blank-like treatment; grouped bar chart of mean median signal (cps) per C-code and adsorbent series.
4. **Type-2 (N > 0 — nanoplastics present, no adsorbent):** recovered concentration per C-code with nominal reference lines and a deviation table.
5. When both types are present a **cross-type summary table** is shown.
6. The reserved code **`C000`** feeds the optional Co correction — see below.

### 6 — Kinetics tab
1. Select (or upload) kinetics files (`AKN*`).
2. Optionally enable the **Co correction** (see below).
3. The app:
   - Reads counts at the calibration peak wavelength from each spectrum.
   - Inverts the calibration curve to obtain Cₜ.
   - Computes q(t) = (C₀ − Cₜ) / C_ads, where C₀ and C_ads are parsed from the filenames.
   - Averages replicas (R1, R2, …). If the selected data **folder path** contains a sub-folder named `M07`, replicates are instead grouped into **odd** (R1, R3, …) and **even** (R2, R4, …) sub-series.
4. Pseudo-1st and pseudo-2nd order models are fitted per series.
5. Linearised diagnostic plots are shown alongside the main q(t) chart.

### 7 — Isotherms tab
1. Select (or upload) isotherm files (`AIN*`).
2. Optionally enable the **Co correction** (see below).
3. The app:
   - Reads counts at the calibration peak wavelength from each spectrum.
   - Uses the calibration curve to compute **Cₑ** (equilibrium concentration).
   - Computes qₑ = (N − Cₑ) / C_ads and **%R = (N − Cₑ) / N × 100** (percentage removal), where N and C_ads are extracted from the filename.
   - Groups and averages replicates by adsorbent, adsorbent concentration, and (when multiple values are present) pH.
4. Fits **Langmuir** and **Freundlich** models for each series. Additional models (Sips, Tóth, Halsey, Harkins-Jura, Janovics) are available on demand.
5. Displays nonlinear fitted curves (qₑ vs Cₑ) and linearised diagnostic plots.

### 8 — Download Results tab
1. Review the status tiles (one per tab) to confirm which analyses are complete.
2. Click **⬇️ Download results as TXT** for a plain-text summary of all computed parameters.
3. Download individual **CSV files** for every plot and table produced across all tabs.

---

## 🧪 Co correction (optional)

Agitation and centrifugation can shift the *actual* initial nanoplastic
concentration away from the nominal value encoded in the filename. Measuring a
control sample with the **reserved condition code `C000`** (e.g.
`CTRLN100M10pH4ADS0C000R1.txt`) and processing it in the **Control** tab
produces a measured-vs-nominal table, which in turn unlocks an **"Apply Co
correction"** checkbox in the **Kinetics** and **Isotherms** tabs.

When enabled, the measured value replaces the nominal N in the qₑ, %R and q(t)
calculations only — the calibration curve is never modified. Values are matched
in this order:

1. Exact (pH, N) match → its measured value is used directly.
2. Same pH but no exact N → N is scaled by the mean measured/nominal ratio at that pH.
3. No rows for that pH → the whole table is used, ignoring pH.

Whether the correction was applied is recorded in the exported TXT report.

---

## 📐 Kinetic models

**Pseudo-1st-order**:
$$q(t) = q_e \left(1 - e^{-k_1 t}\right)$$

**Pseudo-2nd-order**:
$$q(t) = \frac{k_2\,q_e^2\,t}{1 + k_2\,q_e\,t}$$

| Parameter | Units | Description |
|-----------|-------|-------------|
| qₑ | mg/g | Equilibrium adsorption capacity |
| k₁ | min⁻¹ | Pseudo-1st-order rate constant |
| k₂ | g·mg⁻¹·min⁻¹ | Pseudo-2nd-order rate constant |

---

## 📐 Isotherm models

**Langmuir**:
$$q_e = \frac{q_{max} K_L C_e}{1 + K_L C_e}$$

Separation factor:
$$R_L = \frac{1}{1 + K_L C_0}$$

**Freundlich**:
$$q_e = K_F C_e^{1/n}$$

| Parameter | Units | Description |
|----------|------|-------------|
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

Sips, Tóth, Halsey, Harkins-Jura and Janovics models can additionally be
enabled from the Isotherms tab.

---

## ❓ Troubleshooting

| Symptom | Likely cause |
|---------|-------------|
| Files not listed in the folder selector | Wrong folder path in the sidebar, or filenames lacking the expected prefix |
| Concentration parsed as "–" | Filename does not match the expected pattern (see naming convention below) |
| "Fit did not converge" | Too few time points, or data are noisy / inconsistent |
| Counts at the peak wavelength returns None | Spectrum does not cover the selected wavelength |
| Polynomial degree slider missing | Only two averaged calibration points — a linear fit is forced |
| LOD/LOQ section is disabled | Calibration degree > 1, or fewer than 10 BLK files loaded |
| Negative counts in blanks / DRF / CTRL | Normal for some instruments — the app applies an automatic offset |
| Any analysis tab blocked | Calibration not set yet — go to the Calibration tab first |
| "Apply Co correction" greyed out | No `C000` control file processed in the Control tab |

---

## 🗂 File naming convention

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

### Calibration files — `[H|L]CCN{conc}M{XX}R{n}.txt`
| Token | Meaning | Example |
|-------|---------|---------|
| `CC` / `HCC` / `LCC` | Prefix (optional H = high range, L = low range) | — |
| `N{conc}` | Analyte concentration (mg/L) | `N5`, `N50` |
| `M{XX}` | Matrix / run code | `M00` |
| `R{n}` | Replica number | `R1` |

Examples: `CCN5M00R1.txt` → 5 mg/L,  `HCCN50M00R1.txt` → 50 mg/L (high-range)

### Blank files — `BLK…R{n}.txt`
Any filename starting with `BLK` is treated as a blank.  
Example: `BLKM00R1.txt`, `BLKR3.txt`

### Drift files — `DRFN{N}M{M}t{t}R{n}.txt`
| Token | Meaning | Example |
|-------|---------|---------|
| `N{N}` | NP concentration (mg/L); 0 = Type-1 (no NPs) | `N00`, `N50` |
| `M{M}` | Matrix code | `M00` |
| `t{t}` | Session time-point index | `t1`, `t2`, `t3` |
| `R{n}` | Replica number | `R1` |

Examples: `DRFN00M00t1R1.txt` (Type-1),  `DRFN50M00t2R1.txt` (Type-2)

### Control files — `CTRLN{N}M{M}{AdsCode}{AdsConc}C{code}R{n}.txt`
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

### Kinetics files — `AKN{N}M{M}{AdsCode}{AdsConc}t{time}R{n}.txt`
| Token | Meaning | Example |
|-------|---------|---------|
| `{N}` | Initial analyte concentration (mg/L) | `50`, `150` |
| `M{M}` | Matrix code — **optional** | `M00` |
| `{AdsCode}` | Adsorbent material code (letters only) | `BTO`, `MXB` |
| `{AdsConc}` | Adsorbent concentration (g/L) | `2`, `1` |
| `t{time}` | Contact time in minutes (or hours with p-decimal) | `t05`, `t120`, `t24p0` |
| `R{n}` | Replica number | `R1` |

Examples: `AKN50BTO2t05R1.txt` → N=50 mg/L, AdsCode=BTO, AdsConc=2 g/L, t=5 min, R=1;
`AKN150M00MXB1t420R1.txt` → N=150 mg/L, matrix M00, MXB 1 g/L, t=420 min, R=1

> **Note on time tokens:** plain integers are minutes (`t120` = 120 min);
> p-decimal values are hours (`t24p0` = 24 h = 1440 min, `t0p25` = 15 min).

### Isotherm files — `AIN{N}M{M}{AdsCode}{AdsConc}R{n}.txt`
| Token | Meaning | Example |
|-------|---------|---------|
| `{N}` | Initial concentration (mg/L) | `100`, `200` |
| `M{M}` | Matrix code — **optional** | `M00` |
| `{AdsCode}` | Adsorbent material code (letters only) | `BTO`, `MXB` |
| `{AdsConc}` | Adsorbent concentration (g/L) | `1` |
| `R{n}` | Replica number | `R1` |

Examples: `AIN100BTO1R2.txt` → N=100 mg/L, BTO 1 g/L, replica 2;
`AIN50M05pH4BTO10R1.txt` → N=50 mg/L, matrix M05, pH 4, BTO 10 g/L, replica 1

---

## 🗃️ Repository layout

```
app.py                  Entry point — login, sidebar data source, tab layout
utils.py                Shared helpers: spectrum I/O, filename parsers,
                        calibration inversion, kinetic models, colour scales
tabs/visualization.py   Visualization tab
tabs/calibration.py     Calibration tab
tabs/blank.py           Blank / LOD / LOQ tab
tabs/drift.py           Laser drift monitoring tab
tabs/control.py         Control experiments tab (incl. Co correction)
tabs/kinetics.py        Adsorption kinetics tab
tabs/isotherms.py       Adsorption isotherms tab
tabs/download.py        TXT report and CSV export tab
tabs/help.py            In-app help tab
data-samples/           Example dataset covering every file type
```

---

## 📄 Licence

Released under the [MIT License](LICENSE).

## 📚 Citation

If you use this software in academic work, please cite the accompanying
publication.

<!-- TODO: add the journal reference / DOI here once published. -->
