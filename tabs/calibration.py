"""
tabs/calibration.py — Calibration tab.

Filename convention: [H|L]CCN{N}M{M}R{R}.txt
  e.g.  CCN100M10R1.txt  → N=100, M=10, R=1
        HCCN1M100R2.txt  → N=1,   M=100, R=2
        LCCN5M50R4.txt   → N=5,   M=50,  R=4

Grouping:
  Files are grouped by (N, M).  If multiple replicas (R) exist for the same
  (N, M) pair their counts @ selected wavelength are averaged before building
  the calibration curve.  The polynomial fit (counts → concentration) is
  performed on those averaged points.

Renders:
  1. File selector (multiselect) with per-file summary table.
  2. Wavelength selector slider (default 508 nm if data permits).
  3. Averaged summary table (one row per (N, M) group).
  4. Emission spectra plot (log counts vs wavelength, gray-scale, peak-wl line).
  5. Calibration curve (avg counts @ peak_wl vs concentration N) + polynomial fit.
  6. Full polynomial equation with all coefficients.
  7. "Set calibration" button to lock the model into session_state.
"""

import os
import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.metrics import r2_score, mean_absolute_error

from utils import (
    list_files,
    read_spectral_file,
    read_spectral_bytes,
    parse_filename,
    count_at_wavelength,
    gray_color,
)


# ── Polynomial equation builder ────────────────────────────────────────────────

def _build_poly_equation(model, poly_deg: int) -> str:
    """
    Extract coefficients from the sklearn pipeline and format a pretty string.
    Pipeline must have been built with PolynomialFeatures(include_bias=False).
    """
    lr    = model.named_steps["linearregression"]
    a0    = lr.intercept_
    coefs = lr.coef_          # [a1, a2, …, an]

    sup_map = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

    parts = [f"{a0:.6g}"]
    for i, c in enumerate(coefs, start=1):
        sign  = "+" if c >= 0 else "−"
        power = str(i).translate(sup_map) if i > 1 else ""
        parts.append(f" {sign} {abs(c):.6g}·x{power}")
    return "y = " + "".join(parts)


# ── Main render ────────────────────────────────────────────────────────────────

def render(data_dir: str) -> None:
    """Main entry-point; call from app.py inside the Calibration tab context."""

    # ── 1. File selector ──────────────────────────────────────────────────────
    st.subheader("📂 Input Files — Calibration")

    source = st.radio("File source:", ["📁 Folder", "⬆️ Upload"],
                      horizontal=True, key="cal_source")

    # Build a unified list: {"name": str, "df": DataFrame, "meta": dict|None}
    file_entries = []

    if source == "📁 Folder":
        all_files = list_files(data_dir, r"(CCN|HCCN|LCCN)\d")
        if not all_files:
            st.warning(f"No CC* nor HCC* nor LCC* files found in **{data_dir}**. "
                       "Check the calibration folder path in the sidebar.")
            st.session_state["selected_cal_files"] = []
            return
        saved_cal = st.session_state.get("selected_cal_files", all_files)
        selected = st.multiselect(
            "Choose Calibration (CC, LC, HCC) files:",
            options=all_files,
            default=[f for f in saved_cal if f in all_files] or all_files,
            key="cal_file_selector",
        )
        st.session_state["selected_cal_files"] = selected
        if not selected:
            st.info("Select at least one file to proceed.")
            return
        for fname in selected:
            df = read_spectral_file(os.path.join(data_dir, fname))
            file_entries.append({"name": fname, "df": df,
                                  "meta": parse_filename(fname)})
    else:
        uploaded = st.file_uploader(
            "Upload CCN calibration files (.txt)",
            type="txt", accept_multiple_files=True, key="cal_uploader",
        )
        if not uploaded:
            st.info("Upload at least one CCN calibration file.")
            return
        st.session_state["selected_cal_files"] = [f.name for f in uploaded]
        for f in uploaded:
            df = read_spectral_bytes(f)
            file_entries.append({"name": f.name, "df": df,
                                  "meta": parse_filename(f.name)})

    # ── 2. Wavelength selector ────────────────────────────────────────────────
    # Determine the wavelength range covered by the loaded files.
    _DEFAULT_WL = 508.0
    all_wl: list[float] = []
    for e in file_entries:
        if not e["df"].empty and "wavelength" in e["df"].columns:
            all_wl.extend(e["df"]["wavelength"].dropna().tolist())

    if all_wl:
        wl_min = int(np.floor(min(all_wl)))
        wl_max = int(np.ceil(max(all_wl)))
        # Keep the previously chosen value if it is still in range; otherwise
        # fall back to 508 nm (or the midpoint if 508 is outside the range).
        prev_wl = st.session_state.get("calibration_wavelength", _DEFAULT_WL)
        if wl_min <= prev_wl <= wl_max:
            default_wl = int(round(prev_wl))
        elif wl_min <= _DEFAULT_WL <= wl_max:
            default_wl = int(_DEFAULT_WL)
        else:
            default_wl = (wl_min + wl_max) // 2
    else:
        wl_min, wl_max, default_wl = 400, 800, int(_DEFAULT_WL)

    target_wl: float = float(
        st.slider(
            "Peak wavelength (nm)",
            min_value=wl_min,
            max_value=wl_max,
            value=default_wl,
            step=1,
            key="cal_wavelength_slider",
            help=(
                "Wavelength used to extract counts for the calibration curve. "
                f"Default: {int(_DEFAULT_WL)} nm (if covered by the data)."
            ),
        )
    )
    # Persist so other tabs can read it without re-importing this module.
    st.session_state["calibration_wavelength"] = target_wl

    # Convenience label reused in column names and annotations.
    wl_label = f"Counts @ {target_wl:.0f} nm"

    # ── 3. Per-file summary table ─────────────────────────────────────────────
    rows = []
    for e in file_entries:
        meta = e["meta"]
        df   = e["df"]
        c_peak = count_at_wavelength(df, target=target_wl)
        rows.append({
            "Filename":                  e["name"],
            "Concentration N (mg/L)":    meta["N"] if meta else "–",
            "Matrix M":                  meta["M"] if meta else "–",
            "Replica R":                 meta["R"] if meta else "–",
            "λ range (nm)":              (f"{df['wavelength'].min():.0f}–"
                                          f"{df['wavelength'].max():.0f}")
                                          if not df.empty else "–",
            "Points":                    len(df),
            wl_label:                    f"{c_peak:.1f}" if c_peak is not None else "–",
        })
    st.dataframe(pd.DataFrame(rows), width='stretch')
    st.divider()

    # ── 4. Parse files & group by (N, M), averaging replicas ─────────────────
    raw_records = []
    skipped     = []
    for e in file_entries:
        meta = e["meta"]
        if meta is None:
            skipped.append(
                f"{e['name']} (cannot parse filename — expected [H|L]CCN<N>M<M>R<R>.txt)"
            )
            continue
        c_peak = count_at_wavelength(e["df"], target=target_wl)
        if c_peak is None:
            skipped.append(f"{e['name']} (no data near {target_wl:.0f} nm)")
            continue
        raw_records.append({
            "Filename":  e["name"],
            "pH":        meta["pH"],
            "N (mg/L)":  meta["N"],
            "M":         meta["M"],
            "R":         meta["R"],
            wl_label:    c_peak,
            "df":        e["df"],   # keep spectrum for plotting
        })

    for s in skipped:
        st.warning(f"Skipped: {s}")

    if not raw_records:
        st.warning("None of the selected files could be parsed. "
                   "Expected filenames like CCN100M10R1.txt, HCCN5M50R2.txt, …")
        return

    # ── 5. Averaged table (one row per (N, M) group) ──────────────────────────

    # Detect if there are multiple pH values
    ph_values = {r["pH"] for r in raw_records if r["pH"] is not None}
    multi_ph  = len(ph_values) > 1

    if multi_ph:
        st.info(f"🧪 Detected **{len(ph_values)} pH values**: {sorted(ph_values)}. "
                "Data will also be grouped by pH.")

    raw_df = pd.DataFrame([
        {k: v for k, v in r.items() if k != "df"} for r in raw_records
    ])

    grp_cols = (["pH", "N (mg/L)", "M"] if multi_ph else ["N (mg/L)", "M"])

    avg_df = (
        raw_df.groupby(grp_cols, as_index=False)
        .agg(
            Counts_mean=(wl_label, "mean"),
            Counts_std =(wl_label, "std"),
            Reps       =("R",      "count"),
        )
        .rename(columns={"Counts_mean": f"Avg {wl_label}"})
        .sort_values(grp_cols)
        .reset_index(drop=True)
    )

    with st.expander("📋 Averaged data (mean per N × M group)", expanded=True):
        st.dataframe(avg_df.round(3), width='stretch')

    # ── 6. Emission spectra ───────────────────────────────────────────────────
    st.subheader("Calibration Spectra")
    sorted_records = sorted(raw_records, key=lambda r: r["N (mg/L)"])
    n_spec = len(sorted_records)

    fig1 = go.Figure()
    for i, rec in enumerate(sorted_records):
        df = rec["df"]

        # Filter out zero-count points only for plotting
        df_plot = df[df["counts"] != 0]
        if df_plot.empty:
            continue

        log_c = np.log10(df_plot["counts"])
        label = f"N={rec['N (mg/L)']:.4g} M={rec['M']:.4g} R={rec['R']}"

        fig1.add_trace(go.Scatter(
            x=df_plot["wavelength"],
            y=log_c,
            mode="lines",
            name=label,
            line=dict(color=gray_color(i, n_spec), width=2),
        ))

    fig1.add_vline(
        x=target_wl,
        line=dict(color="red", width=2, dash="dash"),
        annotation_text=f"{target_wl:.0f} nm",
        annotation_font_color="red",
        annotation_position="top right",
    )
    fig1.update_layout(
        xaxis_title="Wavelength (nm)", yaxis_title="log₁₀(Counts)",
        legend_title="File", template="plotly_white", height=460,
    )
    st.plotly_chart(fig1, width='stretch')

    # ── 7. Calibration curve + polynomial fit ─────────────────────────────────
    st.subheader("Calibration Curve")

    avg_counts_col = f"Avg {wl_label}"
    concs_arr  = avg_df["N (mg/L)"].values.astype(float)
    counts_arr = avg_df[avg_counts_col].values.astype(float)
    counts_std = avg_df["Counts_std"].fillna(0).values.astype(float)

    n_points = len(avg_df)
    max_deg  = n_points - 1

    if max_deg < 1:
        st.warning("Only 1 averaged data point available for calibration. "
                   "At least 2 distinct (N, M) groups are needed.")
        return
    elif max_deg == 1:
        poly_deg = 1
        st.info("ℹ️ Only 2 averaged data points — using a linear fit (degree 1).")
    else:
        poly_deg = st.selectbox(
            "Polynomial degree (n)",
            options=list(range(1, max_deg + 1)),
            index=0,
            key="cal_poly_deg",
        )

    # Fit model
    cal_model = make_pipeline(
        PolynomialFeatures(poly_deg, include_bias=False),
        LinearRegression(fit_intercept=True),
    )
    cal_model.fit(concs_arr.reshape(-1, 1), counts_arr)

    x_fit  = np.linspace(concs_arr.min(), concs_arr.max(), 400).reshape(-1, 1)
    y_fit  = cal_model.predict(x_fit)
    y_pred = cal_model.predict(concs_arr.reshape(-1, 1))
    r2c  = float(r2_score(counts_arr, y_pred))         if n_points > 1 else float("nan")
    maec = float(mean_absolute_error(counts_arr, y_pred)) if n_points > 1 else 0.0

    eq_str = _build_poly_equation(cal_model, poly_deg)

    has_err = bool((counts_std > 0).any())
    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(
        x=concs_arr, y=counts_arr,
        error_y=dict(type="data", array=counts_std, visible=has_err),
        mode="markers", name="Averaged data",
        marker=dict(color="steelblue", size=9),
    ))
    fig2.add_trace(go.Scatter(
        x=x_fit.flatten(), y=y_fit, mode="lines",
        name=f"Poly deg={poly_deg}  |  R²={r2c:.4f}  |  MAE={maec:.2f}",
        line=dict(color="crimson", width=2.5),
    ))
    fig2.update_layout(
        xaxis_title="Concentration N (mg/L)",
        yaxis_title=f"Avg {wl_label}",
        template="plotly_white", height=400,
    )
    st.plotly_chart(fig2, width='stretch')

    st.markdown(
        f"**Calibration equation:**  `{eq_str}` "
        f"(where `x` = concentration (mg/L) and `y` = {wl_label})"
    )

    col1, col2, col3 = st.columns(3)
    col1.metric("Polynomial degree", poly_deg)
    col2.metric("R²",  f"{r2c:.5f}" if not np.isnan(r2c) else "N/A (1 point)")
    col3.metric("MAE", f"{maec:.4f}")

    # ── 8. Set calibration ────────────────────────────────────────────────────
    st.divider()
    if st.button("✅  Set calibration", type="primary", width='stretch'):
        st.session_state["calibration_model"]        = cal_model
        st.session_state["calibration_poly_deg"]     = poly_deg
        st.session_state["calibration_r2"]           = r2c
        st.session_state["calibration_mae"]          = maec
        st.session_state["calibration_equation"]     = eq_str
        st.session_state["calibration_conc_range"]   = (float(concs_arr.min()),
                                                         float(concs_arr.max()))
        st.session_state["calibration_source_dir"]   = data_dir
        st.session_state["calibration_wavelength"]   = target_wl   # already set above
        # Calibration curve data (for CSV export)
        cal_curve_df = pd.DataFrame({
            "Concentration_N_mg_per_L": concs_arr,
            f"Avg_Counts_{target_wl:.0f}nm": counts_arr,
            "Counts_std":               counts_std,
            "Fitted_Counts":            cal_model.predict(concs_arr.reshape(-1, 1)),
        })
        st.session_state["cal_curve_df"] = cal_curve_df
        st.success(f"✅ Calibration locked — degree {poly_deg} | "
                   f"R²={r2c:.5f} | MAE={maec:.4f} | λ={target_wl:.0f} nm")

    # Status banner
    if "calibration_model" in st.session_state:
        locked_deg = st.session_state["calibration_poly_deg"]
        locked_r2  = st.session_state["calibration_r2"]
        locked_eq  = st.session_state.get("calibration_equation", "")
        locked_wl  = st.session_state.get("calibration_wavelength", _DEFAULT_WL)
        st.info(
            f"🔒 **Active calibration** — degree **{locked_deg}** | "
            f"R²={locked_r2:.5f} | λ={locked_wl:.0f} nm\n\n`{locked_eq}`"
        )

    # LOD / LOQ has moved to the dedicated "Blank" tab (tabs/blank.py).
