"""
tabs/control.py — Control experiments tab.

Filename convention:  CTRL N<N> M<M> <AdsCode><AdsConc> C<Ccode> R<R> .ext
  e.g.  CTRLN00M00BTO3C01R1.txt  → Type-1 (N=0, adsorbent BTO 3 g/L), condition C01
        CTRLN50M00ADS0C02R1.txt  → Type-2 (N=50 mg/L, no adsorbent),  condition C02

Purpose:
  Evaluates the effect of different experimental factors (centrifugation speed, pH,
  temperature, adsorbent concentration, etc.) on the fluorescence signal.
  Each Cxx code represents one specific operational condition defined by the user.

Type-1 (N = 0) — adsorbent present, no nanoplastics:
  - Blank-like spectra: apply negative-count offset if needed.
  - Extract the MEDIAN of counts in a configurable wavelength window (default 495–525 nm).
  - Average medians across replicas for each (AdsCode, AdsConc, Ccode) group.
  - Bar chart: mean of medians per Ccode, grouped by adsorbent series.

Type-2 (N > 0) — nanoplastics present, no adsorbent:
  - Count at 508 nm → average across replicas per (N, AdsCode, AdsConc, Ccode).
  - Convert to concentration (mg/L) via calibration curve.
  - Bar chart: recovered concentration per Ccode, grouped by N series.
  - Only concentration in mg/L is reported (not raw cps) for Type-2.

Cross-type view:
  When both types are present the section ends with a side-by-side summary table.

Co correction (optional) — CTRL ...C000... :
  Reserved condition code C000 marks a special control sample used to check
  whether agitation/centrifugation altered the nominal nanoplastic
  concentration. Its measured concentration is stored in
  st.session_state["co_correction_df"] so the Kinetics and Isotherms tabs can
  optionally use it in place of the nominal N (qe, %R, q(t) only — not
  applied to calibration curves).
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
    parse_ph_prefix,
    blue_color,
    gray_color,
)


# ── Filename parser ────────────────────────────────────────────────────────────

def _parse_ctrl(fname: str) -> dict | None:
    """
    Parse a CTRL filename. Returns a dict or None.
    CTRL N<N> M<M> <AdsCode><AdsConc> C<Ccode> R<R>
    Accepts p-decimal notation for N, M, AdsConc.
    Ccode is kept as a zero-padded string (e.g. '01', '12').
    """
    base = os.path.basename(fname)
    m = _re.search(
        r"CTRL"
        r"N(\d+(?:p\d+)?)"
        r"M(\d+(?:p\d+)?)"
        r"(?:pH\d+(?:p\d+)?)?"
        r"([A-Za-z]+)"
        r"(\d+(?:p\d+)?)"
        r"C(\d+)"
        r"R(\d+)",
        base, _re.IGNORECASE,
    )
    if not m:
        return None

    def _p(s: str) -> float:
        return float(s.replace("p", "."))

    return {
        "N":       _p(m.group(1)),
        "M":       _p(m.group(2)),
        "AdsCode": m.group(3).upper(),
        "AdsConc": _p(m.group(4)),
        "Ccode":   m.group(5),       # kept as string — preserves leading zeros
        "R":       int(m.group(6)),
        "pH":      parse_ph_prefix(base),
    }


# ── Signal helpers (shared with drift pattern) ────────────────────────────────

def _apply_offset(df: pd.DataFrame, fname: str) -> tuple[pd.DataFrame, float]:
    """Shift counts up by |min| if any negative values exist. Prints st.info."""
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
    """Median of counts in [wl_lo, wl_hi] nm, or None if window is empty."""
    mask = (df["wavelength"] >= wl_lo) & (df["wavelength"] <= wl_hi)
    sub  = df.loc[mask, "counts"]
    return float(sub.median()) if not sub.empty else None


def _series_label(ads_code: str, ads_conc: float) -> str:
    """Human-readable adsorbent series label."""
    if ads_conc == 0:
        return f"{ads_code} (no adsorbent)"
    return f"{ads_code} {ads_conc:g} g/L"


def _x_label(ccode: str, ccode_labels: dict) -> str:
    """
    Return the display label for a C-code.
    Uses the user-provided description when available, otherwise falls back to
    the plain 'Cxx' label.  Passed explicitly so there is no closure dependency.
    """
    default = _ccode_label(ccode)
    user    = ccode_labels.get(ccode, default)
    return f"{default}: {user}" if user != default else default


def _ccode_label(ccode: str) -> str:
    """Display label for a C-code, e.g. '01' → 'C01'."""
    return f"C{ccode}"


# ── Plotly bar-chart builder ───────────────────────────────────────────────────

def _bar_chart(
    grouped_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    err_col: str,
    group_col: str,
    y_label: str,
    title: str,
    n_colors: int,
    color_fn=blue_color,
    text_col: str | None = None,
) -> go.Figure:
    """
    Grouped bar chart. Each unique value of *group_col* becomes a separate trace.
    x_col    : categorical axis (C-codes)
    y_col    : bar height
    err_col  : std — shown as error bars
    group_col: the column that defines bar groups (series)
    text_col : optional column whose values are shown as annotations above each bar
    """
    fig = go.Figure()
    groups = sorted(grouped_df[group_col].unique())

    for i, grp in enumerate(groups):
        sub = grouped_df[grouped_df[group_col] == grp].sort_values(x_col)
        col = color_fn(i, max(n_colors, 2))
        bar_text = sub[text_col].tolist() if text_col and text_col in sub.columns else None
        fig.add_trace(go.Bar(
            x=sub[x_col].tolist(),
            y=sub[y_col].tolist(),
            error_y=dict(type="data", array=sub[err_col].tolist(), visible=True),
            name=str(grp),
            marker_color=col,
            text=bar_text,
            textposition="outside",
            textfont=dict(size=12),
            cliponaxis=False,
        ))

    fig.update_layout(
        barmode="group",
        title=title,
        xaxis_title="Control condition (C-code)",
        yaxis_title=y_label,
        template="plotly_white",
        height=460,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        uniformtext=dict(mode="hide", minsize=8),
    )
    return fig


# ── Main render ────────────────────────────────────────────────────────────────

def render(data_dir: str) -> None:
    """Entry-point called from app.py inside the Control tab context."""

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

    st.subheader("📂 Input Files — Control")

    st.markdown(
        #"Evaluates the effect of different experimental factors on the fluorescence "
        #"signal. Each **C-code** identifies one specific operational condition "
        #"(e.g. centrifugation speed, pH, temperature, adsorbent concentration).\n\n"
        "| Type | Filename example | Treatment |\n"
        "|------|------------------|-----------|\n"
        "| **Type 1** | `CTRLN00M00BTO3C01R1` | N=0, adsorbent present — "
        "blank-like; median in λ window |\n"
        "| **Type 2** | `CTRLN50M00ADS0C02R1` | N>0, no adsorbent — "
        f"peak at {_peak_wl:.0f} nm → concentration (mg/L) |\n"
        "| **Co correction** | `CTRLN100M10pH4ADS0C000R1` | Reserved code "
        "`C000` — optional nominal-concentration check |"
    )

    # ── 1. File source ─────────────────────────────────────────────────────────
    source = st.radio("File source:", ["📁 Folder", "⬆️ Upload"],
                      horizontal=True, key="ctrl_source")

    file_entries: list[dict] = []

    if source == "📁 Folder":
        ctrl_files = list_files(data_dir, r"CTRL")
        if not ctrl_files:
            st.warning(f"No CTRL* files found in **{data_dir}**.")
            return
        saved = st.session_state.get("_ctrl_sel", ctrl_files)
        selected = st.multiselect(
            "Choose Control (CTRL) files:",
            options=ctrl_files,
            default=[f for f in saved if f in ctrl_files] or ctrl_files,
            key="_ctrl_selector",
        )
        st.session_state["_ctrl_sel"] = selected
        if not selected:
            st.info("Select at least one CTRL file.")
            return
        for fname in selected:
            df = read_spectral_file(os.path.join(data_dir, fname))
            file_entries.append({"name": fname, "df": df})
    else:
        uploaded = st.file_uploader(
            "Upload CTRL spectral files",
            type=["txt", "asc", "akn"],
            accept_multiple_files=True,
            key="_ctrl_uploader",
        )
        if not uploaded:
            st.info("Upload at least one CTRL file.")
            return
        st.session_state["_ctrl_sel"] = [f.name for f in uploaded]
        for f in uploaded:
            df = read_spectral_bytes(f)
            file_entries.append({"name": f.name, "df": df})

    # ── 2. Parse filenames ─────────────────────────────────────────────────────
    parsed: list[dict] = []
    skipped: list[str] = []
    for e in file_entries:
        meta = _parse_ctrl(e["name"])
        if meta is None:
            skipped.append(e["name"])
        else:
            parsed.append({**e, "meta": meta})

    for s in skipped:
        st.warning(
            f"Skipped — cannot parse: **{s}**  "
            "(expected `CTRLN<N>M<M><AdsCode><AdsConc>C<Ccode>R<R>.ext`)"
        )
    if not parsed:
        return

    # ── 3. Per-file summary table ──────────────────────────────────────────────
    summary_rows = []
    for e in parsed:
        m = e["meta"]
        summary_rows.append({
            "Filename":  e["name"],
            "N (mg/L)":  m["N"],
            "M":         m["M"],
            "AdsCode":   m["AdsCode"],
            "AdsConc":   m["AdsConc"],
            "C-code":    _ccode_label(m["Ccode"]),
            "R":         m["R"],
            "Type":      "Co correction (C000)" if m["Ccode"] == "000"
                         else "Type 1 (N=0)" if m["N"] == 0
                         else f"Type 2 (N={m['N']} mg/L)",
        })
    with st.expander("📋 Loaded CTRL files", expanded=True):
        st.dataframe(pd.DataFrame(summary_rows), width='stretch')

    # ── 4. Wavelength window for Type-1 ───────────────────────────────────────
    st.divider()
    st.markdown("**Type-1 wavelength window** (for median extraction):")
    _wc1, _wc2 = st.columns(2)
    wl_lo = _wc1.number_input("λ min (nm)", value=495.0, step=1.0, key="ctrl_wl_lo")
    wl_hi = _wc2.number_input("λ max (nm)", value=525.0, step=1.0, key="ctrl_wl_hi")

    # ── Optional C-code annotations ───────────────────────────────────────────
    all_ccodes = sorted({e["meta"]["Ccode"] for e in parsed})
    ccode_labels: dict[str, str] = {cc: _ccode_label(cc) for cc in all_ccodes}  # defaults

    with st.expander("🏷️  Optional: add a description for each C-code", expanded=False):
        st.markdown(
            "Provide a short label for each C-code to appear in charts. "
            "Leave blank to use the default `Cxx` label."
        )
        for cc in all_ccodes:
            label = st.text_input(
                f"C{cc} description", value="", key=f"ctrl_clabel_{cc}"
            )
            if label.strip():
                ccode_labels[cc] = label.strip()

    # ── Calibration availability ───────────────────────────────────────────────
    cal_ready = "calibration_model" in st.session_state
    if not cal_ready:
        st.info(
            "ℹ️ Set a calibration in the **Calibration** tab to also view "
            "Type-2 results as **concentration (mg/L)**."
        )

    co_entries = [e for e in parsed if e["meta"]["Ccode"] == "000"]
    type1 = [e for e in parsed if e["meta"]["N"] == 0 and e["meta"]["Ccode"] != "000"]
    type2 = [e for e in parsed if e["meta"]["N"] != 0 and e["meta"]["Ccode"] != "000"]

    # ══════════════════════════════════════════════════════════════════════════
    # TYPE 1 — blank-like (N = 0, adsorbent present)
    # ══════════════════════════════════════════════════════════════════════════
    t1_avg_df: pd.DataFrame | None = None

    if type1:
        st.divider()
        st.subheader("Type 1 — CTRL without nanoplastics (N = 0)")

        t1_records: list[dict] = []
        for e in type1:
            df_e, off = _apply_offset(e["df"], e["name"])
            med = _median_in_window(df_e, wl_lo, wl_hi)
            if med is None:
                st.warning(
                    f"No data in {wl_lo}–{wl_hi} nm for **{e['name']}** — skipped."
                )
                continue
            m = e["meta"]
            t1_records.append({
                "Filename":       e["name"],
                "N (mg/L)":       m["N"],
                "M":              m["M"],
                "AdsCode":        m["AdsCode"],
                "AdsConc":        m["AdsConc"],
                "Ads_series":     _series_label(m["AdsCode"], m["AdsConc"]),
                "Ccode":          m["Ccode"],
                "x_label":        _x_label(m["Ccode"], ccode_labels),
                "R":              m["R"],
                "Median (cps)":   med,
                "Offset applied": off,
            })

        if not t1_records:
            st.warning("No usable Type-1 data.")
        else:
            t1_file_df = pd.DataFrame(t1_records)
            with st.expander("📋 Type-1  per-file medians", expanded=False):
                st.dataframe(
                    t1_file_df.drop(columns=["x_label"]).round(4),
                    width='stretch',
                )

            # ── Average medians per (AdsCode, AdsConc, Ccode) ─────────────
            t1_avg = (
                t1_file_df
                .groupby(["Ads_series", "Ccode", "x_label"], as_index=False)
                .agg(
                    Mean_median=("Median (cps)", "mean"),
                    Std_median =("Median (cps)", "std"),
                    n_replicas =("R",            "count"),
                    N_val      =("N (mg/L)",     "first"),
                    M_val      =("M",            "first"),
                    AdsCode_val=("AdsCode",      "first"),
                    AdsConc_val=("AdsConc",      "first"),
                )
                .sort_values(["Ads_series", "Ccode"])
                .reset_index(drop=True)
            )
            t1_avg["Std_median"] = t1_avg["Std_median"].fillna(0)
            t1_avg["bar_label"] = t1_avg.apply(
                lambda r: f"N{r['N_val']:02.0f}M{r['M_val']:02.0f}{r['AdsCode_val']}{r['AdsConc_val']:g}",    
                #lambda r: f"N{r['N_val']:02.0f}M{r['M_val']:02.0f}ADS{r['AdsConc_val']:g}",
                axis=1,
            )
            t1_avg_df = t1_avg

            with st.expander(
                "📋 Type-1  mean of medians per condition (Ccode × Adsorbent)",
                expanded=True,
            ):
                st.dataframe(t1_avg.round(4), width='stretch')

            # ── Bar chart: signal (cps) ────────────────────────────────────
            st.markdown("#### Type 1 — Signal per control condition (cps)")
            n_series_t1 = len(t1_avg["Ads_series"].unique())
            fig_t1 = _bar_chart(
                grouped_df=t1_avg,
                x_col="x_label",
                y_col="Mean_median",
                err_col="Std_median",
                group_col="Ads_series",
                y_label=f"Mean of medians — {wl_lo:.0f}–{wl_hi:.0f} nm  (cps)",
                title="Type-1 CTRL: blank signal per condition",
                n_colors=n_series_t1,
                color_fn=gray_color,
                text_col="bar_label",
            )
            st.plotly_chart(fig_t1, width='stretch')

            st.session_state.setdefault("ctrl_results", {})
            st.session_state["ctrl_results"].update({
                "type1_per_file": t1_file_df.drop(columns=["x_label"]),
                "type1_avg":      t1_avg,
                "wl_window":      (wl_lo, wl_hi),
            })

    # ══════════════════════════════════════════════════════════════════════════
    # TYPE 2 — nanoplastics present (N > 0), no adsorbent
    # ══════════════════════════════════════════════════════════════════════════
    t2_avg_df: pd.DataFrame | None = None

    if type2 or co_entries:
        st.divider()
        st.subheader("Type 2 — CTRL with nanoplastics (N > 0)")

        t2_records: list[dict] = []
        for e in type2 + co_entries:
            df_e, off = _apply_offset(e["df"], e["name"])
            c508 = count_at_wavelength(df_e, target=_peak_wl)
            if c508 is None:
                st.warning(f"No data near {_peak_wl:.0f} nm for **{e['name']}** — skipped.")
                continue
            m = e["meta"]
            t2_records.append({
                "Filename":         e["name"],
                "N (mg/L)":         m["N"],
                "M":                m["M"],
                "AdsCode":          m["AdsCode"],
                "AdsConc":          m["AdsConc"],
                "Ads_series":       _series_label(m["AdsCode"], m["AdsConc"]),
                "N_series":         f"N={m['N']:g} mg/L",
                "Ccode":            m["Ccode"],
                "x_label":          _x_label(m["Ccode"], ccode_labels),
                "R":                m["R"],
                f"Counts @ {_peak_wl:.0f} nm":  c508,
                "Offset applied":   off,
            })

        if not t2_records:
            st.warning("No usable Type-2 data.")
        else:
            t2_file_df = pd.DataFrame(t2_records)
            with st.expander(f"📋 Type-2  per-file counts @ {_peak_wl:.0f} nm", expanded=False):
                st.dataframe(
                    t2_file_df.drop(columns=["x_label"]).round(4),
                    width='stretch',
                )

            # ── Average across replicas per (N, Ads, Ccode) ───────────────
            t2_avg = (
                t2_file_df
                .groupby(["N_series", "Ads_series", "Ccode", "x_label"], as_index=False)
                .agg(
                    Avg_counts=(f"Counts @ {_peak_wl:.0f} nm", "mean"),
                    Std_counts=(f"Counts @ {_peak_wl:.0f} nm", "std"),
                    n_replicas=("R",               "count"),
                    N_val     =("N (mg/L)",        "first"),
                    M_val     =("M",               "first"),
                    AdsCode_val=("AdsCode",        "first"),
                    AdsConc_val=("AdsConc",        "first"),
                )
                .sort_values(["N_series", "Ads_series", "Ccode"])
                .reset_index(drop=True)
            )
            t2_avg["Std_counts"] = t2_avg["Std_counts"].fillna(0)
            t2_avg["bar_label"] = t2_avg.apply(
                lambda r: f"N{r['N_val']:02.0f}M{r['M_val']:02.0f}{r['AdsCode_val']}{r['AdsConc_val']:g}",    
                #lambda r: f"N{r['N_val']:02.0f}M{r['M_val']:02.0f}ADS{r['AdsConc_val']:g}",
                axis=1,
            )
            t2_avg_df = t2_avg

            with st.expander(
                "📋 Type-2  averaged counts per (N × Condition) group",
                expanded=True,
            ):
                st.dataframe(t2_avg.round(4), width='stretch')

            # ── Convert to concentration ───────────────────────────────────
            if cal_ready:
                cal_model  = st.session_state["calibration_model"]
                conc_range = st.session_state["calibration_conc_range"]

                conc_vals = np.array([
                    invert_calibration(cal_model, v, conc_range)
                    for v in t2_avg["Avg_counts"].values
                ])
                conc_hi = np.array([
                    invert_calibration(cal_model, v + s, conc_range)
                    for v, s in zip(
                        t2_avg["Avg_counts"].values,
                        t2_avg["Std_counts"].values,
                    )
                ])
                conc_err = np.abs(conc_hi - conc_vals)

                t2_avg = t2_avg.copy()
                t2_avg["Conc_calc (mg/L)"]    = conc_vals
                t2_avg["Conc_err (mg/L)"]     = conc_err
                t2_avg_df = t2_avg

                # ── Bar chart: concentration (mg/L) ────────────────────────
                st.markdown("#### Type 2 — Recovered concentration per control condition")

                # Build a temp df with error column for the chart helper
                _t2_chart = t2_avg.copy()
                _t2_chart["_err"] = conc_err
                n_series_t2 = len(t2_avg["N_series"].unique())

                fig_t2 = _bar_chart(
                    grouped_df=_t2_chart,
                    x_col="x_label",
                    y_col="Conc_calc (mg/L)",
                    err_col="_err",
                    group_col="N_series",
                    y_label="Recovered concentration (mg/L)",
                    title="Type-2 CTRL: recovered concentration per condition",
                    n_colors=n_series_t2,
                    color_fn=blue_color,
                    text_col="bar_label",
                )

                # Add nominal N horizontal reference lines
                for N_val_str in t2_avg["N_series"].unique():
                    try:
                        N_val = float(N_val_str.split("=")[1].split(" ")[0])
                        fig_t2.add_hline(
                            y=N_val,
                            line=dict(color="gray", dash="dot", width=1),
                            annotation_text=f"Nominal {N_val_str}",
                            annotation_position="right",
                        )
                    except (IndexError, ValueError):
                        pass

                st.plotly_chart(fig_t2, width='stretch')

                # ── Deviation from nominal ────────────────────────────────
                n_nominal_map = {
                    row["N_series"]: float(row["N_series"].split("=")[1].split(" ")[0])
                    for _, row in t2_avg.iterrows()
                    if "=" in row["N_series"]
                }
                dev_rows = []
                for _, row in t2_avg.iterrows():
                    N_nom = n_nominal_map.get(row["N_series"], float("nan"))
                    cv    = float(row["Conc_calc (mg/L)"])
                    dev   = abs(cv - N_nom)
                    pct   = 100.0 * dev / N_nom if N_nom and N_nom != 0 else float("nan")
                    dev_rows.append({
                        "N_nominal (mg/L)":  N_nom,
                        "Ads_series":        row["Ads_series"],
                        "C-code":            _ccode_label(row["Ccode"]),
                        "Conc_calc (mg/L)":  round(cv, 5),
                        "Dev_from_nom (mg/L)": round(dev, 5),
                        "Dev_%":             round(pct, 3),
                    })
                dev_df = pd.DataFrame(dev_rows)
                with st.expander(
                    "📋 Type-2  deviation from nominal concentration",
                    expanded=True,
                ):
                    st.dataframe(dev_df, width='stretch')

                st.session_state.setdefault("ctrl_results", {})
                st.session_state["ctrl_results"]["type2_deviation"] = dev_df
            else:
                st.info(
                    "Set a calibration to convert Type-2 counts to concentration (mg/L)."
                )

            st.session_state.setdefault("ctrl_results", {})
            st.session_state["ctrl_results"].update({
                "type2_per_file": t2_file_df.drop(columns=["x_label"]),
                "type2_avg":      t2_avg_df,
            })

    # ══════════════════════════════════════════════════════════════════════════
    # Co CORRECTION — special control C000 (optional nominal-concentration fix)
    # ══════════════════════════════════════════════════════════════════════════
    if co_entries:
        st.divider()
        st.subheader("🧪 Co correction — initial concentration check (C000)")
        st.markdown(
            "Special control (reserved code **`C000`**) used to check whether "
            "agitation/centrifugation altered the nominal nanoplastic "
            "concentration. The measured concentration found here can "
            "**optionally** replace the nominal N in the **Kinetics** and "
            "**Isotherms** tabs — for the qe, %R and q(t) calculations only. "
            "*Not applied to calibration curves.*"
        )

        if not cal_ready:
            st.info("Set a calibration to compute the corrected concentration.")
        else:
            co_records = []
            for e in co_entries:
                df_e, off = _apply_offset(e["df"], e["name"])
                c508 = count_at_wavelength(df_e, target=_peak_wl)
                if c508 is None:
                    st.warning(f"No data near {_peak_wl:.0f} nm for **{e['name']}** — skipped.")
                    continue
                m = e["meta"]
                co_records.append({
                    "Filename":  e["name"],
                    "pH":        m["pH"],
                    "N_nominal": m["N"],
                    "R":         m["R"],
                    f"Counts @ {_peak_wl:.0f} nm": c508,
                    "Offset applied": off,
                })

            if not co_records:
                st.warning("No usable Co-correction data.")
            else:
                co_file_df = pd.DataFrame(co_records)
                with st.expander("📋 Co correction — per-file counts", expanded=False):
                    st.dataframe(co_file_df.round(4), width='stretch')

                co_avg = (
                    co_file_df
                    .groupby(["pH", "N_nominal"], as_index=False, dropna=False)
                    .agg(
                        Avg_counts=(f"Counts @ {_peak_wl:.0f} nm", "mean"),
                        Std_counts=(f"Counts @ {_peak_wl:.0f} nm", "std"),
                        n_replicas=("R", "count"),
                    )
                    .sort_values(["N_nominal"])
                    .reset_index(drop=True)
                )
                co_avg["Std_counts"] = co_avg["Std_counts"].fillna(0)
                co_avg["N_measured"] = [
                    invert_calibration(cal, v, crange) for v in co_avg["Avg_counts"]
                ]
                co_avg["Ratio (measured/nominal)"] = co_avg["N_measured"] / co_avg["N_nominal"]
                co_avg["Deviation (%)"] = (co_avg["Ratio (measured/nominal)"] - 1.0) * 100.0

                with st.expander(
                    "📋 Co correction — measured vs. nominal concentration",
                    expanded=True,
                ):
                    st.dataframe(co_avg.round(4), width='stretch')

                st.session_state["co_correction_df"] = co_avg[
                    ["pH", "N_nominal", "N_measured"]
                ].copy()

                st.session_state.setdefault("ctrl_results", {})
                st.session_state["ctrl_results"]["co_correction"] = co_avg

                st.success(
                    "✅ Co correction available. Enable it in the **Kinetics** "
                    "and **Isotherms** tabs (optional checkbox) to apply it to "
                    "qe, %R and q(t)."
                )

    # ══════════════════════════════════════════════════════════════════════════
    # CROSS-TYPE SUMMARY (both types present)
    # ══════════════════════════════════════════════════════════════════════════
    if type1 and type2 and t1_avg_df is not None and t2_avg_df is not None:
        st.divider()
        st.subheader("📊 Cross-type summary")
        st.markdown(
            "Both Type-1 and Type-2 files were loaded. "
            "Below is a combined summary of all conditions. "
            "Units differ between types: Type-1 in **cps** (signal), "
            "Type-2 in **mg/L** (concentration, requires calibration)."
        )

        t1_summary = t1_avg_df[
            ["Ads_series", "Ccode", "Mean_median", "Std_median", "n_replicas"]
        ].copy()
        t1_summary.insert(0, "Type", "Type 1")
        t1_summary.rename(columns={
            "Mean_median": "Value",
            "Std_median":  "Std",
        }, inplace=True)
        t1_summary["Unit"] = "cps (median signal)"
        t1_summary["Ccode"] = t1_summary["Ccode"].apply(_ccode_label)

        if cal_ready and "Conc_calc (mg/L)" in t2_avg_df.columns:
            t2_summary = t2_avg_df[
                ["N_series", "Ccode", "Conc_calc (mg/L)", "Conc_err (mg/L)", "n_replicas"]
            ].copy()
            t2_summary.insert(0, "Type", "Type 2")
            t2_summary.rename(columns={
                "N_series":         "Ads_series",
                "Conc_calc (mg/L)": "Value",
                "Conc_err (mg/L)":  "Std",
            }, inplace=True)
            t2_summary["Unit"] = "mg/L (recovered conc.)"
            t2_summary["Ccode"] = t2_summary["Ccode"].apply(_ccode_label)
        else:
            t2_summary = t2_avg_df[
                ["N_series", "Ccode", "Avg_counts", "Std_counts", "n_replicas"]
            ].copy()
            t2_summary.insert(0, "Type", "Type 2")
            t2_summary.rename(columns={
                "N_series":   "Ads_series",
                "Avg_counts": "Value",
                "Std_counts": "Std",
            }, inplace=True)
            t2_summary["Unit"] = f"cps (counts @ {_peak_wl:.0f} nm)"
            t2_summary["Ccode"] = t2_summary["Ccode"].apply(_ccode_label)

        cross_df = pd.concat(
            [t1_summary, t2_summary], ignore_index=True
        ).round(4)
        with st.expander("📋 Cross-type condition summary", expanded=True):
            st.dataframe(cross_df, width='stretch')

        st.session_state.setdefault("ctrl_results", {})
        st.session_state["ctrl_results"]["cross_summary"] = cross_df

    # ── Completion notice ──────────────────────────────────────────────────────
    if st.session_state.get("ctrl_results"):
        st.divider()
        st.success(
            "✅ Control analysis complete. "
            "Results and CSV downloads are available in the **Download Results** tab."
        )
