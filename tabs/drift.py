"""
tabs/drift.py — Laser Drift Monitoring tab.

Filename convention:  DRFN<N>M<M>t<t>R<R>.ext
  e.g.  DRFN00M00t1R1.txt  → Type-1 (N=0, no NPs), session t=1, replica R=1
        DRFN50M00t2R1.txt  → Type-2 (N=50 mg/L),   session t=2, replica R=1

Purpose:
  Monitors possible drift of the excitation laser intensity over the course of
  a measurement session.  DRF files are taken at least at the beginning (t=1)
  and end (t=last) of each session, and optionally in the middle (t=2, …).

Type-1 (N = 0) — blank-like files:
  - Apply negative-count offset if needed (same rule as BLK in calibration).
  - Extract the MEDIAN of counts inside a configurable wavelength window
    (default 495–525 nm) — more robust than mean for noisy blank spectra.
  - Average those medians across replicas for each session time-point t.
  - Plot "mean of medians" vs t  →  signal (cps) + linear drift trend.
  - Convert to concentration (mg/L) via calibration curve + plot vs t.

Type-2 (N > 0) — files with nanoplastics:
  - Apply negative-count offset if needed.
  - Extract count at 508 nm (fluorescence peak).
  - Average across replicas for each (N, M, t) group.
  - Plot avg counts vs t per N-series  →  signal (cps) + linear drift trend.
  - Convert to concentration (mg/L) + plot vs t; compare with nominal N.
"""

import os
import re as _re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from utils import (
    list_files,
    read_spectral_file,
    read_spectral_bytes,
    invert_calibration,
    count_at_wavelength,
)


# ── Filename parser (local, self-contained) ────────────────────────────────────

def _parse_drf(fname: str) -> dict | None:
    """
    Parse a DRF filename.  Returns a dict with keys N, M, t, R or None.
    Accepts p-decimal notation: DRFN10p5M00t2R1 → N=10.5.
    """
    base = os.path.basename(fname)
    m = _re.search(
        r"DRF"
        r"N(\d+(?:p\d+)?)"
        r"M(\d+(?:p\d+)?)"
        r"(?:pH\d+(?:p\d+)?)?"
        r"t(\d+)"
        r"R(\d+)",
        base, _re.IGNORECASE,
    )
    if not m:
        return None

    def _p(s: str) -> float:
        return float(s.replace("p", "."))

    return {
        "N": _p(m.group(1)),
        "M": _p(m.group(2)),
        "t": int(m.group(3)),
        "R": int(m.group(4)),
    }


# ── Signal helpers ─────────────────────────────────────────────────────────────

def _apply_offset(df: pd.DataFrame, fname: str) -> tuple[pd.DataFrame, float]:
    """
    If the counts column has negatives, shift all values up by |min|.
    Prints an st.info message and returns (corrected_df, offset).
    Returns (df, 0.0) when no offset is needed.
    """
    df = df.copy()
    if df.empty:
        return df, 0.0
    min_c = float(df["counts"].min())
    if min_c < 0:
        offset = abs(min_c)
        df["counts"] = df["counts"] + offset
        return df, offset
    return df, 0.0


def _median_in_window(df: pd.DataFrame, wl_lo: float, wl_hi: float) -> float | None:
    """Return the median of counts in [wl_lo, wl_hi] nm, or None if empty."""
    mask = (df["wavelength"] >= wl_lo) & (df["wavelength"] <= wl_hi)
    sub  = df.loc[mask, "counts"]
    return float(sub.median()) if not sub.empty else None


def _linear_trend(x, y) -> tuple[float, float, float]:
    """Fit y = slope·x + intercept and return (slope, intercept, R²)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or np.all(x == x[0]):
        return float("nan"), float("nan"), float("nan")
    coeffs            = np.polyfit(x, y, 1)
    slope, intercept  = float(coeffs[0]), float(coeffs[1])
    y_pred            = np.polyval(coeffs, x)
    ss_res            = float(np.sum((y - y_pred) ** 2))
    ss_tot            = float(np.sum((y - y.mean()) ** 2))
    r2                = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return slope, intercept, r2


def _add_trend_trace(fig, x_vals, slope, intercept, name, color):
    """Overlay a dashed linear-trend line on a Plotly figure."""
    if np.isnan(slope):
        return
    x_arr = np.asarray(x_vals, dtype=float)
    y_fit = slope * x_arr + intercept
    fig.add_trace(go.Scatter(
        x=x_arr.tolist(), y=y_fit.tolist(),
        mode="lines",
        name=name,
        line=dict(color=color, dash="dash", width=1.5),
    ))


def _xaxis_tick_cfg(t_vals):
    """Return xaxis tick kwargs for integer time-point labels."""
    tv = sorted(set(int(t) for t in t_vals))
    return dict(tickvals=tv, ticktext=[str(t) for t in tv])


def _drift_warning(conc_first: float, conc_last: float, label: str = "") -> None:
    """Show a colored Streamlit alert for the relative concentration drift first→last."""
    if conc_first == 0:
        st.info(
            f"{'**' + label + '**: ' if label else ''}"
            "Cannot compute relative drift — first concentration is zero."
        )
        return
    rel_pct = (conc_last - conc_first) / abs(conc_first) * 100
    abs_pct = abs(rel_pct)
    sign = "+" if rel_pct >= 0 else ""
    prefix = f"**{label}**: " if label else ""
    msg = f"{prefix}Concentration drift (last vs first): **{sign}{rel_pct:.2f}%**"
    if abs_pct < 5.0:
        st.success(msg)
    elif abs_pct < 10.0:
        st.warning(msg)
    else:
        st.error(msg)


# ── Main render ────────────────────────────────────────────────────────────────

def render(data_dir: str) -> None:
    """Entry-point called from app.py inside the Drift tab context."""

    # ── Guard ──────────────────────────────────────────────────────────────────
    if "calibration_model" not in st.session_state:
        st.warning("⚠️ No calibration set yet.  "
                   "Go to the **Calibration** tab, choose a polynomial degree "
                   "and click **Set calibration**.")
        return

    cal     = st.session_state["calibration_model"]
    cal_deg = st.session_state["calibration_poly_deg"]
    cal_r2  = st.session_state["calibration_r2"]
    cal_eq  = st.session_state.get("calibration_equation", "")
    crange  = st.session_state["calibration_conc_range"]
    cal_dir = st.session_state.get("calibration_source_dir", "–")
    _peak_wl = st.session_state.get("calibration_wavelength", 508.0) 
    st.success(
    f"🔒 Using calibration: polynomial degree {cal_deg} | "
    f"R²={cal_r2:.5f}\n\n"
    f"`{cal_eq}`\n\n"
    f"📁 Calibration folder: `{cal_dir}`"
    )

    st.divider()

    st.subheader("📂 Input Files — Drift")

    st.markdown(
        #"Analyses **DRF** files that are measured at the beginning, (optionally) "
        #"the middle, and the end of a measurement session to detect possible drift "
        #"of the excitation laser intensity.\n\n"
        "| Type | Prefix example | Treatment |\n"
        "|------|----------------|-----------|\n"
        "| **Type 1** | `DRFN00M00t1R1` | No NPs — blank-like; median in λ window |\n"
        f"| **Type 2** | `DRFN50M00t2R1` | With NPs — peak at {_peak_wl:.0f} nm |"
    )

    # ── 1. File source ─────────────────────────────────────────────────────────
    source = st.radio("File source", ["📁 Folder", "⬆️ Upload"],
                      horizontal=True, key="drf_source")

    file_entries: list[dict] = []

    if source == "📁 Folder":
        drf_files = list_files(data_dir, r"DRF")
        if not drf_files:
            st.warning(f"No DRF* files found in **{data_dir}**.")
            return
        saved = st.session_state.get("_drf_sel", drf_files)
        selected = st.multiselect(
            "Choose DRF files:",
            options=drf_files,
            default=[f for f in saved if f in drf_files] or drf_files,
            key="_drf_selector",
        )
        st.session_state["_drf_sel"] = selected
        if not selected:
            st.info("Select at least one DRF file.")
            return
        for fname in selected:
            df = read_spectral_file(os.path.join(data_dir, fname))
            file_entries.append({"name": fname, "df": df})
    else:
        uploaded = st.file_uploader(
            "Upload DRF spectral files",
            type=["txt", "asc", "akn"],
            accept_multiple_files=True,
            key="_drf_uploader",
        )
        if not uploaded:
            st.info("Upload at least one DRF file.")
            return
        st.session_state["_drf_sel"] = [f.name for f in uploaded]
        for f in uploaded:
            df = read_spectral_bytes(f)
            file_entries.append({"name": f.name, "df": df})

    # ── 2. Parse filenames & build per-file summary ───────────────────────────
    parsed: list[dict] = []
    skipped: list[str] = []
    for e in file_entries:
        meta = _parse_drf(e["name"])
        if meta is None:
            skipped.append(e["name"])
        else:
            parsed.append({**e, "meta": meta})

    for s in skipped:
        st.warning(
            f"Skipped — cannot parse filename: **{s}**  "
            "(expected `DRFN<N>M<M>t<t>R<R>.ext`)"
        )
    if not parsed:
        return

    summary_rows = [
        {
            "Filename":  e["name"],
            "N (mg/L)":  e["meta"]["N"],
            "M":         e["meta"]["M"],
            "t":         e["meta"]["t"],
            "R":         e["meta"]["R"],
            "Type":      "Type 1 (no NPs)" if e["meta"]["N"] == 0
                         else f"Type 2  N={e['meta']['N']} mg/L",
        }
        for e in parsed
    ]
    with st.expander("📋 Loaded DRF files", expanded=True):
        st.dataframe(pd.DataFrame(summary_rows), width='stretch')

    # ── 3. Wavelength window for Type-1 median ────────────────────────────────
    st.divider()
    st.markdown("**Type-1 wavelength window** (for median extraction):")
    _wc1, _wc2 = st.columns(2)
    wl_lo = _wc1.number_input("λ min (nm)", value=495.0, step=1.0, key="drf_wl_lo")
    wl_hi = _wc2.number_input("λ max (nm)", value=525.0, step=1.0, key="drf_wl_hi")

    # ── Calibration availability ──────────────────────────────────────────────
    cal_ready = "calibration_model" in st.session_state
    if not cal_ready:
        st.info(
            "ℹ️ Set a calibration in the **Calibration** tab to also view "
            "drift expressed in **concentration units (mg/L)**."
        )

    type1 = [e for e in parsed if e["meta"]["N"] == 0]
    type2 = [e for e in parsed if e["meta"]["N"] != 0]

    # ════════════════════════════════════════════════════════════════════════════
    # TYPE 1 — blank-like (N = 0)
    # ════════════════════════════════════════════════════════════════════════════
    if type1:
        st.divider()
        st.subheader("Type 1 — DRF without nanoplastics (N = 0)")

        # ── Per-file: offset correction + median extraction ────────────────
        t1_records: list[dict] = []
        for e in type1:
            df_e, off = _apply_offset(e["df"], e["name"])
            med = _median_in_window(df_e, wl_lo, wl_hi)
            if med is None:
                st.warning(
                    f"No data in {wl_lo}–{wl_hi} nm for **{e['name']}** — skipped."
                )
                continue
            t1_records.append({
                "Filename":          e["name"],
                "t":                 e["meta"]["t"],
                "R":                 e["meta"]["R"],
                "M":                 e["meta"]["M"],
                "Median (cps)":      med,
                "Offset applied":    off,
            })

        if not t1_records:
            st.warning("No usable Type-1 data.")
        else:
            t1_file_df = pd.DataFrame(t1_records)
            with st.expander("📋 Type-1  per-file medians", expanded=False):
                st.dataframe(t1_file_df.round(4), width='stretch')

            # ── Mean of medians per time-point ─────────────────────────────
            t1_avg = (
                t1_file_df.groupby("t", as_index=False)
                .agg(
                    Mean_of_medians=("Median (cps)", "mean"),
                    Std_of_medians =("Median (cps)", "std"),
                    n_replicas     =("R",            "count"),
                )
                .sort_values("t")
                .reset_index(drop=True)
            )
            t1_avg["Std_of_medians"] = t1_avg["Std_of_medians"].fillna(0)

            with st.expander("📋 Type-1  mean of medians per time-point", expanded=True):
                st.dataframe(t1_avg.round(4), width='stretch')

            t_v   = t1_avg["t"].values.astype(float)
            sig_v = t1_avg["Mean_of_medians"].values.astype(float)
            std_v = t1_avg["Std_of_medians"].values.astype(float)

            slope_s1, int_s1, r2_s1 = _linear_trend(t_v, sig_v)

            # ── Plot 1a: signal (cps) ──────────────────────────────────────
            st.markdown("#### Type 1 — Signal drift (cps)")
            fig1a = go.Figure()
            fig1a.add_trace(go.Scatter(
                x=t_v.tolist(), y=sig_v.tolist(),
                error_y=dict(type="data", array=std_v.tolist(), visible=True),
                mode="markers+lines",
                name="Mean of medians",
                marker=dict(size=9, color="steelblue"),
                line=dict(color="steelblue"),
            ))
            _add_trend_trace(
                fig1a, t_v, slope_s1, int_s1,
                f"Linear trend  (slope = {slope_s1:.4f} cps/step,  R² = {r2_s1:.4f})",
                "crimson",
            )
            fig1a.update_layout(
                xaxis=dict(title="Session time-point (t)", **_xaxis_tick_cfg(t_v)),
                yaxis_title=f"Mean of medians — {wl_lo:.0f}–{wl_hi:.0f} nm  (cps)",
                template="plotly_white", height=380,
                legend=dict(orientation="h", yanchor="bottom", y=1.02,
                            xanchor="left", x=0),
            )
            st.plotly_chart(fig1a, width='stretch')
            st.caption(
                f"Linear drift trend: **slope = {slope_s1:.5f} cps / step** | "
                f"intercept = {int_s1:.4f} | R² = {r2_s1:.4f}"
            )

            # ── Plot 1b: concentration (mg/L) ──────────────────────────────
            if cal_ready:
                st.markdown("#### Type 1 — Equivalent concentration drift (mg/L)")
                cal_model  = st.session_state["calibration_model"]
                conc_range = st.session_state["calibration_conc_range"]

                # Type-1 signals are blank-like (near-zero concentration).
                # The calibration conc_range starts at the lowest calibration
                # standard (e.g. 5 mg/L), so invert_calibration would clamp
                # every blank value to that minimum.
                # Extend the lower bound to 0 so the inversion can resolve
                # near-zero equivalent concentrations correctly.
                t1_conc_range = (0.0, conc_range[1] * 1.5)       
                conc_v  = np.array([
                    invert_calibration(cal_model, v, t1_conc_range) for v in sig_v
                ])
                conc_hi = np.array([
                    invert_calibration(cal_model, v + s, t1_conc_range)
                    for v, s in zip(sig_v, std_v)
                ])
                conc_err = np.abs(conc_hi - conc_v)

                slope_c1, int_c1, r2_c1 = _linear_trend(t_v, conc_v)

                fig1b = go.Figure()
                fig1b.add_trace(go.Scatter(
                    x=t_v.tolist(), y=conc_v.tolist(),
                    error_y=dict(type="data", array=conc_err.tolist(), visible=True),
                    mode="markers+lines",
                    name="Equivalent concentration",
                    marker=dict(size=9, color="teal"),
                    line=dict(color="teal"),
                ))
                _add_trend_trace(
                    fig1b, t_v, slope_c1, int_c1,
                    f"Linear trend  (slope = {slope_c1:.4f} mg/L/step,  "
                    f"R² = {r2_c1:.4f})",
                    "darkorange",
                )

                # LOD / LOQ reference lines if computed
                lod_data = st.session_state.get("lod_loq", {})
                if lod_data:
                    fig1b.add_hline(
                        y=lod_data["LOD"],
                        line=dict(color="gray", dash="dot", width=1.5),
                        annotation_text="LOD", annotation_position="right",
                    )
                    fig1b.add_hline(
                        y=lod_data["LOQ"],
                        line=dict(color="lightcoral", dash="dot", width=1.5),
                        annotation_text="LOQ", annotation_position="right",
                    )

                fig1b.update_layout(
                    xaxis=dict(title="Session time-point (t)", **_xaxis_tick_cfg(t_v)),
                    yaxis_title="Concentration (mg/L)",
                    template="plotly_white", height=380,
                    legend=dict(orientation="h", yanchor="bottom", y=1.02,
                                xanchor="left", x=0),
                )
                st.plotly_chart(fig1b, width='stretch')
                st.caption(
                    f"Linear drift trend: **slope = {slope_c1:.5f} mg/L / step** | "
                    f"intercept = {int_c1:.4f} | R² = {r2_c1:.4f}"
                )
                if len(conc_v) >= 2:
                    _drift_warning(float(conc_v[0]), float(conc_v[-1]), label="Type 1")

                t1_avg_out             = t1_avg.copy()
                t1_avg_out["Conc (mg/L)"] = conc_v
                t1_avg_out["t"]        = t1_avg_out["t"].astype(int)
            else:
                slope_c1 = int_c1 = r2_c1 = float("nan")
                conc_v   = np.full(len(sig_v), float("nan"))
                t1_avg_out = t1_avg.copy()
                t1_avg_out["t"] = t1_avg_out["t"].astype(int)

            # ── Persist Type-1 results ─────────────────────────────────────
            st.session_state.setdefault("drift_results", {})
            st.session_state["drift_results"].update({
                "type1_per_file":   t1_file_df,
                "type1_avg":        t1_avg_out,
                "type1_trend_sig":  {
                    "slope":     slope_s1,
                    "intercept": int_s1,
                    "r2":        r2_s1,
                    "unit":      "cps/step",
                },
                "type1_trend_conc": {
                    "slope":     slope_c1,
                    "intercept": int_c1,
                    "r2":        r2_c1,
                    "unit":      "mg/L/step",
                },
                "wl_window": (wl_lo, wl_hi),
            })

    # ════════════════════════════════════════════════════════════════════════════
    # TYPE 2 — with nanoplastics (N > 0)
    # ════════════════════════════════════════════════════════════════════════════
    if type2:
        st.divider()
        st.subheader("Type 2 — DRF with nanoplastics (N > 0)")

        # ── Per-file: offset correction + 508 nm extraction ───────────────
        t2_records: list[dict] = []
        for e in type2:
            df_e, off = _apply_offset(e["df"], e["name"])
            c508 = count_at_wavelength(df_e, target=_peak_wl)
            if c508 is None:
                st.warning(f"No data near {_peak_wl:.0f} nm for **{e['name']}** — skipped.")
                continue
            t2_records.append({
                "Filename":          e["name"],
                "N (mg/L)":          e["meta"]["N"],
                "M":                 e["meta"]["M"],
                "t":                 e["meta"]["t"],
                "R":                 e["meta"]["R"],
                f"Counts @ {_peak_wl:.0f} nm":   c508,
                "Offset applied":    off,
            })

        if not t2_records:
            st.warning("No usable Type-2 data.")
        else:
            t2_file_df = pd.DataFrame(t2_records)
            with st.expander(f"📋 Type-2  per-file counts @ {_peak_wl:.0f} nm", expanded=False):
                st.dataframe(t2_file_df.round(4), width='stretch')

            # ── Average across replicas per (N, M, t) ─────────────────────
            t2_avg = (
                t2_file_df.groupby(["N (mg/L)", "M", "t"], as_index=False)
                .agg(
                    Avg_counts=(f"Counts @ {_peak_wl:.0f} nm", "mean"),
                    Std_counts=(f"Counts @ {_peak_wl:.0f} nm", "std"),
                    n_replicas=("R",               "count"),
                )
                .sort_values(["N (mg/L)", "t"])
                .reset_index(drop=True)
            )
            t2_avg["Std_counts"] = t2_avg["Std_counts"].fillna(0)

            with st.expander("📋 Type-2  averaged counts per (N, t) group", expanded=True):
                st.dataframe(t2_avg.round(4), width='stretch')

            n_groups = sorted(t2_avg["N (mg/L)"].unique())

            # ── Plot 2a: signal (cps) per N-series ────────────────────────
            st.markdown("#### Type 2 — Signal drift (cps)")
            fig2a = go.Figure()
            t2_trend_sig: dict = {}

            for i, N_val in enumerate(n_groups):
                sub   = t2_avg[t2_avg["N (mg/L)"] == N_val].sort_values("t")
                t_v   = sub["t"].values.astype(float)
                sig_v = sub["Avg_counts"].values.astype(float)
                std_v = sub["Std_counts"].values.astype(float)
                col   = "steelblue"

                fig2a.add_trace(go.Scatter(
                    x=t_v.tolist(), y=sig_v.tolist(),
                    error_y=dict(type="data", array=std_v.tolist(), visible=True),
                    mode="markers+lines",
                    name=f"N = {N_val} mg/L",
                    marker=dict(size=9, color=col),
                    line=dict(color=col),
                ))

                slope_s, int_s, r2_s = _linear_trend(t_v, sig_v)
                t2_trend_sig[N_val]  = {
                    "slope":     slope_s,
                    "intercept": int_s,
                    "r2":        r2_s,
                    "unit":      "cps/step",
                }
                _add_trend_trace(
                    fig2a, t_v, slope_s, int_s,
                    f"Trend N={N_val}  (slope={slope_s:.4f},  R²={r2_s:.4f})",
                    col,
                )

            fig2a.update_layout(
                xaxis=dict(title="Session time-point (t)",
                           **_xaxis_tick_cfg(t2_avg["t"].values)),
                yaxis_title=f"Avg Counts @ {_peak_wl:.0f} nm (cps)",
                template="plotly_white", height=400,
                legend=dict(orientation="h", yanchor="bottom", y=1.02,
                            xanchor="left", x=0),
            )
            st.plotly_chart(fig2a, width='stretch')

            # Trend summary caption
            for N_val, tr in t2_trend_sig.items():
                st.caption(
                    f"N = {N_val} mg/L — slope = **{tr['slope']:.5f} cps/step** | "
                    f"intercept = {tr['intercept']:.4f} | R² = {tr['r2']:.4f}"
                )

            # ── Plot 2b: concentration (mg/L) per N-series ────────────────
            t2_trend_conc: dict   = {}
            t2_deviation_rows: list[dict] = []
            t2_drift_fl: dict = {}

            if cal_ready:
                st.markdown("#### Type 2 — Equivalent concentration drift (mg/L)")
                cal_model  = st.session_state["calibration_model"]
                conc_range = st.session_state["calibration_conc_range"]
                t2_conc_range = (0.0, conc_range[1] * 1.5)
                fig2b = go.Figure()

                for i, N_val in enumerate(n_groups):
                    sub   = t2_avg[t2_avg["N (mg/L)"] == N_val].sort_values("t")
                    t_v   = sub["t"].values.astype(float)
                    sig_v = sub["Avg_counts"].values.astype(float)
                    std_v = sub["Std_counts"].values.astype(float)
                    col   = "steelblue"

                    conc_v  = np.array([
                        invert_calibration(cal_model, v, t2_conc_range) for v in sig_v
                    ])
                    conc_hi = np.array([
                        invert_calibration(cal_model, v + s, t2_conc_range)
                        for v, s in zip(sig_v, std_v)
                    ])
                    conc_err = np.abs(conc_hi - conc_v)

                    slope_c, int_c, r2_c = _linear_trend(t_v, conc_v)
                    t2_trend_conc[N_val] = {
                        "slope":     slope_c,
                        "intercept": int_c,
                        "r2":        r2_c,
                        "unit":      "mg/L/step",
                    }
                    if len(conc_v) >= 2:
                        t2_drift_fl[N_val] = (float(conc_v[0]), float(conc_v[-1]))

                    fig2b.add_trace(go.Scatter(
                        x=t_v.tolist(), y=conc_v.tolist(),
                        error_y=dict(type="data", array=conc_err.tolist(), visible=True),
                        mode="markers+lines",
                        name=f"N = {N_val} mg/L",
                        marker=dict(size=9, color=col),
                        line=dict(color=col),
                    ))
                    _add_trend_trace(
                        fig2b, t_v, slope_c, int_c,
                        f"Trend N={N_val}  (slope={slope_c:.4f},  R²={r2_c:.4f})",
                        col,
                    )

                    # Nominal concentration reference line
                    fig2b.add_hline(
                        y=N_val,
                        line=dict(color=col, dash="dot", width=1),
                        annotation_text=f"Nominal {N_val} mg/L",
                        annotation_position="right",
                    )

                    # Collect deviation rows
                    for t_pt, cv, sv in zip(t_v, conc_v, sig_v):
                        dev = abs(cv - N_val)
                        pct = 100.0 * dev / N_val if N_val != 0 else float("nan")
                        t2_deviation_rows.append({
                            "N_nominal (mg/L)":  N_val,
                            "t":                 int(t_pt),
                            "Avg_counts (cps)":  round(float(sv), 4),
                            "Conc_calc (mg/L)":  round(float(cv), 5),
                            "Dev_from_nom (mg/L)": round(float(dev), 5),
                            "Dev_%":             round(float(pct), 3),
                        })

                fig2b.update_layout(
                    xaxis=dict(title="Session time-point (t)",
                               **_xaxis_tick_cfg(t2_avg["t"].values)),
                    yaxis_title="Concentration (mg/L)",
                    template="plotly_white", height=400,
                    legend=dict(orientation="h", yanchor="bottom", y=1.02,
                                xanchor="left", x=0),
                )
                st.plotly_chart(fig2b, width='stretch')

                for N_val, tr in t2_trend_conc.items():
                    st.caption(
                        f"N = {N_val} mg/L — slope = **{tr['slope']:.5f} mg/L/step** | "
                        f"intercept = {tr['intercept']:.4f} | R² = {tr['r2']:.4f}"
                    )
                for N_val, (c0, c1) in t2_drift_fl.items():
                    _drift_warning(c0, c1, label=f"N = {N_val} mg/L")

                # Deviation from nominal table
                if t2_deviation_rows:
                    dev_df = pd.DataFrame(t2_deviation_rows)
                    with st.expander(
                        "📋 Type-2  deviation of recovered concentration from nominal",
                        expanded=True,
                    ):
                        st.dataframe(dev_df, width='stretch')
                    st.session_state.setdefault("drift_results", {})
                    st.session_state["drift_results"]["type2_deviation"] = dev_df

            # ── Persist Type-2 results ─────────────────────────────────────
            st.session_state.setdefault("drift_results", {})
            st.session_state["drift_results"].update({
                "type2_per_file":   t2_file_df,
                "type2_avg":        t2_avg,
                "type2_trend_sig":  t2_trend_sig,
                "type2_trend_conc": t2_trend_conc,
            })

    # ── Completion notice ──────────────────────────────────────────────────────
    if "drift_results" in st.session_state and st.session_state["drift_results"]:
        st.divider()
        st.success(
            "✅ Drift analysis complete. "
            "Results and CSV downloads are available in the **Download Results** tab."
        )
