"""
utils.py — shared helpers for parsing, file I/O, kinetic models and calibration inversion.
"""

import os
import re                          # module-level regex helpers
import re as _re                   # alias used inside list_files to avoid shadowing
import numpy as np
import pandas as pd
import streamlit as st
from scipy.optimize import curve_fit
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.linear_model import LinearRegression


# ── Recognised file extensions ────────────────────────────────────────────────
# MUST be a tuple (not a set): str.endswith() only accepts a str or a tuple.

_SPECTRAL_EXTENSIONS = (".txt", ".asc", ".akn")


# ═════════════════════════════════════════════════════════════════════════════
# ASC format detection and reading
# ═════════════════════════════════════════════════════════════════════════════

def _is_asc_format(text):
    """
    Return True if *text* looks like an ASC file.
    Criterion: the first non-empty line has exactly 4 comma-separated fields,
    all convertible to int.
    TXT files carry a text header, so the conversion fails for them.
    """
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) != 4:
            return False
        try:
            for p in parts:
                int(p)
            return True     # first non-empty line → 4 integers → ASC
        except ValueError:
            return False    # contains text → standard TXT format
    return False


def _parse_asc_line(line):
    """
    Parse one ASC line of 4 comma-separated integer fields:
        wl_int , wl_dec ,   counts_int , counts_dec

    The comma acts simultaneously as column separator and as decimal point.
    The real value is reconstructed as:
        value = int_part + dec_part / 10**len(dec_str)

    Examples:
        "500,0,   67,849"  ->  (500.0,   67.849)
        "501,0,   67,839"  ->  (501.0,   67.839)
        "508,25,  68,004"  ->  (508.25,  68.004)   <- leading zeros preserved

    Returns (wavelength, counts) as floats, or None if the line is invalid.
    """
    parts = [p.strip() for p in line.split(",")]
    if len(parts) != 4:
        return None
    try:
        wl_i, wl_d, ct_i, ct_d = parts
        wavelength = int(wl_i) + int(wl_d) / (10 ** max(len(wl_d), 1))
        counts     = 1000 * (int(ct_i) + int(ct_d) / (10 ** max(len(ct_d), 1)))
        return wavelength, counts
    except (ValueError, ZeroDivisionError):
        return None


def read_asc_content(text):
    """
    Convert the full text of an ASC file into a DataFrame with columns
    ['wavelength', 'counts'].  Lines that cannot be parsed are ignored.
    """
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        result = _parse_asc_line(line)
        if result is not None:
            rows.append({"wavelength": result[0], "counts": result[1]})
    return pd.DataFrame(rows, columns=["wavelength", "counts"])


# ═════════════════════════════════════════════════════════════════════════════
# File I/O
# ═════════════════════════════════════════════════════════════════════════════

def read_spectral_file(path):
    """
    Read a spectral file in TXT format (with a text header) or ASC format
    (headerless, 4 comma-separated integer fields per line).
    Returns a DataFrame with columns ['wavelength', 'counts'].
    """
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return pd.DataFrame(columns=["wavelength", "counts"])

    if _is_asc_format(text):
        return read_asc_content(text)

    # TXT format with header: look for "number,number,..." lines
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(",")
        if len(parts) < 2:
            continue
        try:
            rows.append({"wavelength": float(parts[0]), "counts": float(parts[1])})
        except ValueError:
            continue    # header / comment line -> ignore
    return pd.DataFrame(rows, columns=["wavelength", "counts"])


def read_spectral_bytes(uploaded_file):
    """
    Variant for files supplied by st.file_uploader.
    Accepts TXT format (with header) or ASC (4 integer fields per line).
    Returns a DataFrame with columns ['wavelength', 'counts'].
    """
    try:
        raw  = uploaded_file.read()
        text = raw.decode("utf-8", errors="replace")
        if hasattr(uploaded_file, "seek"):
            uploaded_file.seek(0)
    except Exception:
        return pd.DataFrame(columns=["wavelength", "counts"])

    if _is_asc_format(text):
        return read_asc_content(text)

    # TXT format with header
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(",")
        if len(parts) < 2:
            continue
        try:
            rows.append({"wavelength": float(parts[0]), "counts": float(parts[1])})
        except ValueError:
            continue
    return pd.DataFrame(rows, columns=["wavelength", "counts"])


def list_files(folder, pattern=None):
    """
    Return the spectral files (.txt, .asc, .akn) inside *folder* whose name
    matches *pattern* (optional regular expression, case-insensitive).
    """
    if not os.path.isdir(folder):
        return []
    files = [
        f for f in os.listdir(folder)
        if f.lower().endswith(_SPECTRAL_EXTENSIONS)
    ]
    if pattern:
        files = [f for f in files if _re.search(pattern, f, _re.IGNORECASE)]
    return sorted(files)


# ═════════════════════════════════════════════════════════════════════════════
# Filename parsers
# ═════════════════════════════════════════════════════════════════════════════

def parse_concentration_ccn(filename):
    """CCN5M00.txt → 5.0,  CCN10p5M00.txt → 10.5."""
    m = re.search(r"CCN(\d+(?:p\d+)?)M", os.path.basename(filename), re.IGNORECASE)
    return _p_to_float(m.group(1)) if m else None


def parse_time_akn(filename):
    """
    AKN50BTO2t120R1.txt → 120.0 min
    AKN50BTO2t24p0R1.txt → 1440.0 min  (24 h expressed as p-decimal)
    """
    m = re.search(r"t(\d+(?:p\d+)?)R", os.path.basename(filename), re.IGNORECASE)
    return _p_time_to_minutes(m.group(1)) if m else None


def count_at_wavelength(df, target=508.0):
    """Return counts at the wavelength nearest to *target*."""
    if df.empty:
        return None
    idx = (df["wavelength"] - target).abs().idxmin()
    return float(df.loc[idx, "counts"])


# ── pH prefix parser ──────────────────────────────────────────────────────────

_PH_RE = re.compile(r"M\d+pH(\d+(?:p\d+)?)", re.IGNORECASE)


def parse_ph_prefix(fname):
    """
    Extract the pH from the optional field placed after the matrix code.
    'AIN10M05pH4BTO10R1.txt' -> 4.0,  'AIN10M05BTO10R1.txt' -> None.
    """
    m = _PH_RE.search(os.path.basename(fname))
    return float(m.group(1).replace("p", ".")) if m else None


# ── Unified filename parser ───────────────────────────────────────────────────
#
# Filename conventions:
#   Calibration : [H|L]CCN<N>M<M>R<R>.txt      e.g. HCCN50M10R1.txt
#   Kinetics    : AKN<N><AdsCode><AdsConc>t<t>R<R>.txt  e.g. AKN50BTO2t05R1.txt
#   Isotherm    : AIN<N><Ads><Cads>R<R>.txt      e.g. AIN100AC10R3.txt
#
# All three accept an optional leading pH prefix (parsed by parse_ph_prefix).


# ── p-decimal helpers ─────────────────────────────────────────────────────────

def _p_to_float(s: str) -> float:
    """Convert a p-decimal token to float: '10p5' → 10.5,  '10' → 10.0."""
    return float(s.replace("p", "."))


def _p_time_to_minutes(s: str) -> float:
    """
    Convert a time token to minutes.
    Tokens containing 'p' represent hours  → multiply by 60.
      e.g. '24p0' → 24.0 h → 1440 min,  '0p25' → 0.25 h → 15 min.
    Plain integer tokens represent minutes directly.
      e.g. '120' → 120 min,  '5' → 5 min.
    """
    val = _p_to_float(s)
    return val * 60.0 if "p" in s.lower() else val


_CCN_RE = re.compile(
    r"(?:H|L)?CC"               # optional H or L prefix, then CC
    r"N(\d+(?:p\d+)?)"          # N  - analyte concentration  (mg/L)
    r"M(\d+(?:p\d+)?)"          # M  - matrix code
    r"(?:pH\d+(?:p\d+)?)?"      # optional pH token  (non-capturing; parsed by parse_ph_prefix)
    r"R(\d+)",                  # R  - replica number
    re.IGNORECASE,
)

_AKN_RE = re.compile(
    r"AKN(\d+(?:p\d+)?)"        # N       - analyte concentration (mg/L)
    r"(?:M\d+(?:p\d+)?)?"       # optional matrix code (non-capturing)
    r"(?:pH\d+(?:p\d+)?)?"      # optional pH token  (non-capturing; parsed by parse_ph_prefix)
    r"([A-Za-z]+)"              # AdsCode - adsorbent identifier (letters only)
    r"(\d+(?:p\d+)?)"           # AdsConc - adsorbent concentration (g/L)
    r"t(\d+(?:p\d+)?)"          # t       - contact time (min or h via p-decimal)
    r"R(\d+)",                  # R       - replica number
    re.IGNORECASE,
)

_AIN_RE = re.compile(
    r"AIN(\d+(?:p\d+)?)"        # N    - initial analyte concentration (mg/L)
    r"(?:M\d+(?:p\d+)?)?"       # optional matrix code (non-capturing)
    r"(?:pH\d+(?:p\d+)?)?"      # optional pH token  (non-capturing; parsed by parse_ph_prefix)
    r"([A-Za-z]+)"              # Ads  - adsorbent code (letters only)
    r"(\d+(?:p\d+)?)"           # Cads - adsorbent concentration (g/L)
    r"R(\d+)",                  # R    - replica number
    re.IGNORECASE,
)

_DRF_RE = re.compile(
    r"DRF"                       # DRF  - drift measurement prefix
    r"N(\d+(?:p\d+)?)"          # N    - nanoplastic concentration (0 = blank-like)
    r"M(\d+(?:p\d+)?)"          # M    - matrix code
    r"(?:pH\d+(?:p\d+)?)?"      # optional pH token  (non-capturing; parsed by parse_ph_prefix)
    r"t(\d+)"                   # t    - session time-point index (1, 2, 3 …)
    r"R(\d+)",                  # R    - replica number
    re.IGNORECASE,
)

_CTRL_RE = re.compile(
    r"CTRL"                      # CTRL    - control measurement prefix
    r"N(\d+(?:p\d+)?)"          # N       - nanoplastic conc (0=Type-1, >0=Type-2)
    r"M(\d+(?:p\d+)?)"          # M       - matrix code
    r"(?:pH\d+(?:p\d+)?)?"      # optional pH token  (non-capturing; parsed by parse_ph_prefix)
    r"([A-Za-z]+)"              # AdsCode - adsorbent identifier (letters only)
    r"(\d+(?:p\d+)?)"           # AdsConc - adsorbent conc; 0 = no adsorbent
    r"C(\d+)"                   # Ccode   - control type / operational condition code
    r"R(\d+)",                  # R       - replica number
    re.IGNORECASE,
)


def parse_filename(fname):
    """
    Detect the file type from *fname* and return a dict of parsed parameters,
    or None if the name matches none of the known conventions.

    The returned dict always contains:
      "type"  - one of "calibration", "kinetics", "isotherm"
      "pH"    - float or None  (from optional Ph<X> prefix)

    Plus type-specific fields — same keys as the original per-tab parsers.
    Time ("t") is always returned in **minutes**.
    """
    base = os.path.basename(fname)
    ph   = parse_ph_prefix(fname)

    m = _CCN_RE.search(base)
    if m:
        return {
            "type": "calibration",
            "pH":   ph,
            "N":    _p_to_float(m.group(1)),
            "M":    _p_to_float(m.group(2)),
            "R":    int(m.group(3)),
        }

    m = _AKN_RE.search(base)
    if m:
        return {
            "type":    "kinetics",
            "pH":      ph,
            "N":       _p_to_float(m.group(1)),
            "AdsCode": m.group(2).upper(),
            "AdsConc": _p_to_float(m.group(3)),
            "t":       _p_time_to_minutes(m.group(4)),   # always minutes
            "R":       int(m.group(5)),
        }

    m = _AIN_RE.search(base)
    if m:
        return {
            "type": "isotherm",
            "pH":   ph,
            "N":    _p_to_float(m.group(1)),
            "Ads":  m.group(2).upper(),
            "Cads": _p_to_float(m.group(3)),
            "R":    int(m.group(4)),
        }

    m = _DRF_RE.search(base)
    if m:
        return {
            "type": "drift",
            "pH":   ph,
            "N":    _p_to_float(m.group(1)),
            "M":    _p_to_float(m.group(2)),
            "t":    int(m.group(3)),
            "R":    int(m.group(4)),
        }

    m = _CTRL_RE.search(base)
    if m:
        return {
            "type":    "control",
            "pH":      ph,
            "N":       _p_to_float(m.group(1)),
            "M":       _p_to_float(m.group(2)),
            "AdsCode": m.group(3).upper(),
            "AdsConc": _p_to_float(m.group(4)),
            "Ccode":   m.group(5),          # kept as string to preserve leading zeros
            "R":       int(m.group(6)),
        }

    return None   # filename not recognised


# ═════════════════════════════════════════════════════════════════════════════
# UI helper
# ═════════════════════════════════════════════════════════════════════════════

def folder_picker(label, default, key):
    """
    Sidebar folder picker: text input + live subfolder selector.
    Returns the chosen path.
    """
    path = st.sidebar.text_input(label, value=default, key=f"{key}_input")

    if os.path.isdir(path):
        try:
            subdirs = sorted(
                d for d in os.listdir(path)
                if os.path.isdir(os.path.join(path, d)) and not d.startswith(".")
            )
        except PermissionError:
            subdirs = []

        if subdirs:
            choice = st.sidebar.selectbox(
                f"Subfolders of `{os.path.basename(path) or path}`",
                options=["— stay here —"] + subdirs,
                key=f"{key}_sub",
            )
            if choice != "— stay here —":
                path = os.path.join(path, choice)
                try:
                    subdirs2 = sorted(
                        d for d in os.listdir(path)
                        if os.path.isdir(os.path.join(path, d)) and not d.startswith(".")
                    )
                except PermissionError:
                    subdirs2 = []
                if subdirs2:
                    choice2 = st.sidebar.selectbox(
                        f"Subfolders of `{os.path.basename(path)}`",
                        options=["— stay here —"] + subdirs2,
                        key=f"{key}_sub2",
                    )
                    if choice2 != "— stay here —":
                        path = os.path.join(path, choice2)
        else:
            st.sidebar.caption("No subfolders found.")
    else:
        st.sidebar.caption("Path does not exist.")

    if path != st.session_state.get(f"{key}_input", path):
        st.sidebar.caption(f"{path}")

    return path


# ═════════════════════════════════════════════════════════════════════════════
# Calibration helpers
# ═════════════════════════════════════════════════════════════════════════════

def invert_calibration(cal_model, counts_value, conc_range, n_points=10000):
    """
    Numerically invert  counts = f(conc)  by scanning a dense grid and finding
    the concentration that gives the closest counts prediction.

    The scan grid is extended well beyond the training conc_range in both
    directions so that cps values that fall outside the calibration range are
    not silently clamped to the nearest boundary.  The lower bound is always
    clamped to 0 (concentration cannot be negative).

    conc_range: (min_conc, max_conc) as a plain tuple — used only to set the
    scale of the extrapolation, not as hard limits.
    """
    span      = conc_range[1] - conc_range[0]
    c_lo      = max(0.0, conc_range[0] - span)   # 1× span below min, ≥ 0
    c_hi      = conc_range[1] + 2 * span          # 2× span above max
    c_scan    = np.linspace(c_lo, c_hi, n_points).reshape(-1, 1)
    counts_pred = cal_model.predict(c_scan).flatten()
    idx         = int(np.argmin(np.abs(counts_pred - counts_value)))
    return float(c_scan[idx, 0])


# ═════════════════════════════════════════════════════════════════════════════
# Co correction (optional) — nominal-concentration fix from CTRL ...C000... files
# ═════════════════════════════════════════════════════════════════════════════

def apply_co_correction(co_df, ph, n_nominal):
    """
    Resolve the corrected nominal concentration to use for *n_nominal* at *ph*,
    using the Co-correction table built in the Control tab from special
    CTRL...C000... files (see tabs/control.py).

    co_df columns: "pH", "N_nominal", "N_measured".

    Resolution order:
      1. Exact (pH, N_nominal) match — return its measured value directly.
      2. Same pH, no exact N match   — scale n_nominal by the average
         measured/nominal ratio found at that pH (the "same percentage
         variation" rule used for isotherm series).
      3. No rows for that pH         — ignore pH and use the whole table
         (covers single-pH experiments where filenames carry no pH token).
      4. Table empty                 — None (caller keeps the nominal value).
    """
    if co_df is None or co_df.empty:
        return None

    sub = co_df[co_df["pH"] == ph] if ph is not None else co_df[co_df["pH"].isna()]
    if sub.empty:
        sub = co_df

    exact = sub[np.isclose(sub["N_nominal"].astype(float), float(n_nominal))]
    if not exact.empty:
        return float(exact["N_measured"].iloc[0])

    ratio = float((sub["N_measured"] / sub["N_nominal"]).mean())
    return float(n_nominal) * ratio


# ═════════════════════════════════════════════════════════════════════════════
# Kinetic models
# ═════════════════════════════════════════════════════════════════════════════

def model_first_order(t, qe, k1):
    return qe * (1.0 - np.exp(-k1 * t))


def model_second_order(t, qe, k2):
    denom = 1.0 + k2 * qe * t
    return (k2 * qe**2 * t) / np.where(denom == 0, 1e-12, denom)


def fit_kinetic(model_fn, t_arr, q_arr, p0, bounds):
    """Wrapper around curve_fit; returns (popt, r2, mae) or (None, None, None)."""
    try:
        popt, _ = curve_fit(model_fn, t_arr, q_arr,
                            p0=p0, maxfev=20_000, bounds=bounds)
        q_pred = model_fn(t_arr, *popt)
        return (popt,
                float(r2_score(q_arr, q_pred)),
                float(mean_absolute_error(q_arr, q_pred)))
    except Exception:
        return None, None, None


def linear_first_order(t_arr, q_arr, qe):
    """
    Linear form of pseudo-1st-order:  ln(qe - q) = ln(qe) - k1*t
    Uses a qe slightly above max(q_arr) when qe leaves no room.
    """
    qe_safe = max(qe, q_arr.max() * 1.001) * 1.001
    y   = np.log(qe_safe - q_arr)
    reg = LinearRegression().fit(t_arr.reshape(-1, 1), y)
    return float(reg.coef_[0]), float(reg.intercept_), qe_safe


def linear_second_order(t_arr, q_arr):
    """Linear form of pseudo-2nd-order:  t/q = 1/(k2*qe^2) + t/qe"""
    mask = q_arr > 1e-12
    if mask.sum() < 2:
        return None, None
    y   = t_arr[mask] / q_arr[mask]
    reg = LinearRegression().fit(t_arr[mask].reshape(-1, 1), y)
    return float(reg.coef_[0]), float(reg.intercept_)


# ═════════════════════════════════════════════════════════════════════════════
# Colour scales
# ═════════════════════════════════════════════════════════════════════════════

def gray_color(i, n):
    """Darkest for the highest concentration (index n-1), lightest for lowest (0)."""
    v = int(40 + (210 - 40) * (1 - i / max(n - 1, 1)))
    return f"rgb({v},{v},{v})"


def blue_color(i, n):
    """Darkest blue for highest concentration, lightest blue for lowest."""
    t = i / max(n - 1, 1)
    r = int(180 - 100 * t)
    g = int(190 - 100 * t)
    b = 255
    return f"rgb({r},{g},{b})"
