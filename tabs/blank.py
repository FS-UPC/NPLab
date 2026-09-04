"""
tabs/blank.py — Blank files & LOD / LOQ tab.

Filename convention: BLK*.txt / BLK*.ASC  (e.g. BLKM00R1.txt)

Purpose:
  Loads blank (BLK) spectral files, applies a negative-count offset where
  necessary, then computes LOD and LOQ from the blank signal statistics and
  the slope of the active linear calibration curve.

Requirements:
  - A linear (degree = 1) calibration must have been set in the Calibration tab.
  - At least 10 valid BLK files are needed for a statistically reliable result.

Outputs (stored in st.session_state["lod_loq"]):
  n, y_mean, y_std, b, LOD, LOQ, blk_df
"""

import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from utils import (
    list_files,
    read_spectral_file,
    read_spectral_bytes,
)


def render(data_dir: str) -> None:
    """Entry-point; call from app.py inside the Blank tab context."""

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

    st.subheader("📂 Input Files — Blank")
    st.markdown(
        "**Requirements:**\n"
        "- A **linear (degree = 1)** calibration set in the **Calibration** tab.\n"
        "- At least **10 valid BLK files**."
    )

    # ── 1. Calibration gate ────────────────────────────────────────────────────
    _cal_ready  = "calibration_model" in st.session_state
    _cal_deg_ok = _cal_ready and (st.session_state["calibration_poly_deg"] == 1)

    if not _cal_ready:
        st.warning(
            "⚠️ No calibration has been set yet. "
            "Go to the **Calibration** tab, fit the curve, and press "
            "**Set calibration** before running the LOD/LOQ calculation."
        )
    elif not _cal_deg_ok:
        st.warning(
            f"⚠️ The active calibration uses polynomial degree "
            f"**{st.session_state['calibration_poly_deg']}**. "
            "LOD / LOQ requires a **linear fit (degree = 1)**. "
            "Return to Calibration, select degree 1 and press *Set calibration*."
        )
    else:
        locked_eq = st.session_state.get("calibration_equation", "N/A")
        st.success(
            f"✅ Linear calibration active — "
            f"R² = {st.session_state['calibration_r2']:.5f}  |  `{locked_eq}`"
        )

    st.divider()

    # ── 2. File source ─────────────────────────────────────────────────────────
    blk_source = st.radio(
        "File source:",
        ["📁 Folder", "⬆️ Upload files"],
        horizontal=True,
        key="blk_source",
    )

    blk_entries: list[dict] = []

    if blk_source == "📁 Folder":
        blk_files = list_files(data_dir, r"BLK")
        if not blk_files:
            st.warning(f"No BLK* files found in **{data_dir}**.")
        else:
            saved_blk = st.session_state.get("_blk_sel_files", blk_files)
            blk_sel = st.multiselect(
                "Choose Blank (BLK) files:",
                options=blk_files,
                default=[f for f in saved_blk if f in blk_files] or blk_files,
                key="_blk_selector",
            )
            st.session_state["_blk_sel_files"] = blk_sel
            for fname in blk_sel:
                _df = read_spectral_file(os.path.join(data_dir, fname))
                blk_entries.append({"name": fname, "df": _df})
    else:
        blk_uploaded = st.file_uploader(
            "Upload BLK files (blanks)",
            type=["txt", "asc", "akn"],
            accept_multiple_files=True,
            key="_blk_uploader",
        )
        if blk_uploaded:
            st.session_state["_blk_sel_files"] = [f.name for f in blk_uploaded]
            for f in blk_uploaded:
                _df = read_spectral_bytes(f)
                blk_entries.append({"name": f.name, "df": _df})

    if not blk_entries:
        st.info("Load BLK files to proceed.")
        return

    # ── 3. Per-file: offset correction + max-count extraction ─────────────────
    blk_counts: list[dict] = []
    blk_skipped: list[str] = []

    for e in blk_entries:
        fname = e["name"]
        df_e  = e["df"].copy()       # work on a copy — never mutate cache

        offset_applied = 0.0
        if not df_e.empty:
            _min_c = float(df_e["counts"].min())
            if _min_c < 0:
                offset_applied = abs(_min_c)
                df_e["counts"] = df_e["counts"] + offset_applied

        _c = float(df_e["counts"].max()) if not df_e.empty else None
        if _c is None:
            blk_skipped.append(fname)
        else:
            blk_counts.append({
                "Filename":       fname,
                "Counts":         _c,
                "offset_applied": offset_applied,
            })

    for s in blk_skipped:
        st.warning(f"Skipped (no usable data): {s}")

    if not blk_counts:
        st.warning("No valid blank data could be extracted.")
        return

    # ── 4. Summary table ───────────────────────────────────────────────────────
    blk_df = pd.DataFrame(blk_counts)
    _disp = blk_df[["Filename", "Counts", "offset_applied"]].rename(
        columns={"offset_applied": "Offset applied (counts)"}
    )
    with st.expander("📋 Blank counts (max peak per BLK file)", expanded=False):
        st.dataframe(_disp.round(4), width='stretch')

    n_blk = len(blk_counts)
    st.caption(f"**{n_blk}** valid blank file(s) loaded.")

    # ── 5. LOD / LOQ ──────────────────────────────────────────────────────────
    st.divider()
    st.subheader("📏 LOD / LOQ Results")

    if not _cal_deg_ok:
        st.info(
            "LOD/LOQ calculation will activate once a linear (degree = 1) "
            "calibration is set in the Calibration tab."
        )
        return

    if n_blk < 10:
        st.warning(
            f"⚠️ At least **10 BLK files** are required — currently {n_blk}. "
            "Load more files to enable the calculation."
        )
        return

    _y_vals = np.array([r["Counts"] for r in blk_counts])
    _y_mean = float(np.mean(_y_vals))
    _y_std  = float(np.std(_y_vals, ddof=1))

    _lr = st.session_state["calibration_model"].named_steps["linearregression"]
    _b  = float(_lr.coef_[0])

    _lod = 3.0  * _y_std / _b
    _loq = 10.0 * _y_std / _b

    # Metric tiles
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("n blanks",       n_blk)
    c2.metric("ȳ_B  (cps)",     f"{_y_mean:.2f}")
    c3.metric("sB   (cps)",     f"{_y_std:.4f}")
    c4.metric("b  (cps·L/mg)", f"{_b:.4f}")

    c5, c6 = st.columns(2)
    c5.metric("LOD  (mg/L)", f"{_lod:.5f}", help="LOD = 3·sB / b")
    c6.metric("LOQ  (mg/L)", f"{_loq:.5f}", help="LOQ = 10·sB / b")

    with st.expander("📐 Formulas and calculation details", expanded=False):
        st.latex(
            r"\bar{y}_B = \frac{1}{n}\sum_{i=1}^{n} y_{B,i}"
            rf"\;=\; {_y_mean:.4f} \text{{ cps}}"
        )
        st.latex(
            r"s_B = \sqrt{\frac{\sum(y_{B,i}-\bar{y}_B)^2}{n-1}}"
            rf"\;=\; {_y_std:.4f} \text{{ cps}}"
        )
        st.latex(
            rf"LOD = \frac{{3\,s_B}}{{b}} = "
            rf"\frac{{3 \times {_y_std:.4f}}}{{{_b:.4f}}}"
            rf"= {_lod:.5f} \text{{ mg/L}}"
        )
        st.latex(
            rf"LOQ = \frac{{10\,s_B}}{{b}} = "
            rf"\frac{{10 \times {_y_std:.4f}}}{{{_b:.4f}}}"
            rf"= {_loq:.5f} \text{{ mg/L}}"
        )

    # ── 6. Scatter plot ────────────────────────────────────────────────────────
    st.markdown("#### Blank signal per file")

    _x  = list(range(1, n_blk + 1))
    _y  = [r["Counts"]         for r in blk_counts]
    _nm = [r["Filename"]       for r in blk_counts]
    _of = [r["offset_applied"] for r in blk_counts]

    _x_norm = [x for x, o in zip(_x,  _of) if o == 0.0]
    _y_norm = [y for y, o in zip(_y,  _of) if o == 0.0]
    _n_norm = [n for n, o in zip(_nm, _of) if o == 0.0]

    _x_off  = [x for x, o in zip(_x,  _of) if o != 0.0]
    _y_off  = [y for y, o in zip(_y,  _of) if o != 0.0]
    _n_off  = [n for n, o in zip(_nm, _of) if o != 0.0]
    _o_off  = [o for o in _of if o != 0.0]

    _ht_off = [
        f"{n}<br>Counts: {y:.2f}<br>Offset: {o:.4f}"
        for n, y, o in zip(_n_off, _y_off, _o_off)
    ]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=_x, y=_y, mode="lines",
        name="_connector", showlegend=False,
        line=dict(color="lightgray", width=1.5),
    ))
    if _x_norm:
        fig.add_trace(go.Scatter(
            x=_x_norm, y=_y_norm, mode="markers",
            name="Blank signal (max counts)",
            marker=dict(size=8, color="steelblue"),
            text=_n_norm,
            hovertemplate="%{text}<br>Counts: %{y:.2f}<extra></extra>",
        ))
    if _x_off:
        fig.add_trace(go.Scatter(
            x=_x_off, y=_y_off, mode="markers",
            name="Blank signal — offset corrected",
            marker=dict(size=9, color="darkorange",
                        symbol="diamond", line=dict(color="black", width=1)),
            text=_ht_off,
            hovertemplate="%{text}<extra></extra>",
        ))

    fig.add_hline(y=_y_mean,
                  line=dict(color="gray",    dash="dash", width=1.5))
    fig.add_hline(y=_y_mean + 3  * _y_std,
                  line=dict(color="orange",  dash="dot",  width=1.5))
    fig.add_hline(y=_y_mean + 10 * _y_std,
                  line=dict(color="crimson", dash="dot",  width=1.5))

    # Dummy legend entries for the reference lines
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="lines",
        name=f"ȳ_B = {_y_mean:.2f}  (blank mean)",
        line=dict(color="gray", dash="dash", width=1.5),
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="lines",
        name=f"ȳ_B + 3·sB = {_y_mean + 3 * _y_std:.2f}  (LOD threshold)",
        line=dict(color="orange", dash="dot", width=1.5),
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="lines",
        name=f"ȳ_B + 10·sB = {_y_mean + 10 * _y_std:.2f}  (LOQ threshold)",
        line=dict(color="crimson", dash="dot", width=1.5),
    ))

    fig.update_layout(
        title="Blank signal — max counts per BLK file",
        xaxis_title="Blank file number",
        yaxis_title="Counts (max peak)",
        legend=dict(),
        template="plotly_white",
        height=420,
    )
    st.plotly_chart(fig, width='stretch')

    # ── 7. Persist results in session state ────────────────────────────────────
    st.session_state["lod_loq"] = {
        "n":      n_blk,
        "y_mean": _y_mean,
        "y_std":  _y_std,
        "b":      _b,
        "LOD":    _lod,
        "LOQ":    _loq,
        "blk_df": blk_df,
    }

    st.success(
        f"✅  LOD = **{_lod:.5f} mg/L**   |   LOQ = **{_loq:.5f} mg/L**"
    )
