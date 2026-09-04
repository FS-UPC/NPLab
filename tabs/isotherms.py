"""
tabs/isotherms.py — Isotherms tab.

Filename convention: AIN{N}{Ads}{y}R{n}.txt
  e.g.  AIN50BTO2R1.txt   →  N=50, Ads=BTO, Cads=2 g/L, R=1
        AIN100AC10R3.txt   →  N=100, Ads=AC, Cads=10 g/L, R=3

  N   = initial analyte concentration (mg/L) — numeric, any number of digits
  Ads = adsorbent code                        — letters only
  y   = adsorbent concentration (g/L)         — numeric, any number of digits
  R   = experiment replica number

Grouping:
  Files are grouped by (N, Ads, y).  If multiple replicas (R) exist for the
  same group, Cₑ and qₑ are averaged across replicas before model fitting.

Each file is an emission spectrum. The calibration curve is used to convert
counts @ 508 nm → Ce (equilibrium concentration). Then:
  qe = (N − Ce) / Cads   [mg/g]   where Cads = y (from filename)

The (Ce, qe) pairs are fitted with Langmuir and Freundlich isotherm models,
one independent fit per (N, Ads, y) group, shown in different colours.
"""

import os
import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.linear_model import LinearRegression

from utils import (
    list_files,
    read_spectral_file,
    parse_filename,
    read_spectral_bytes,
    count_at_wavelength,
    invert_calibration,
    apply_co_correction,
)

# ── Isotherm models ────────────────────────────────────────────────────────────

def langmuir(Ce, qmax, KL):
    return (qmax * KL * Ce) / (1.0 + KL * Ce)

def freundlich(Ce, KF, n):
    return KF * np.power(np.abs(Ce), 1.0 / n)

def _fit(model_fn, Ce, qe, p0, bounds):
    try:
        popt, _ = curve_fit(model_fn, Ce, qe, p0=p0,
                            maxfev=20_000, bounds=bounds)
        q_pred = model_fn(Ce, *popt)
        return popt, float(r2_score(qe, q_pred)), float(mean_absolute_error(qe, q_pred))
    except Exception:
        return None, None, None

# ── Colour palette ─────────────────────────────────────────────────────────────

_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]

def _color(idx: int) -> str:
    return _PALETTE[idx % len(_PALETTE)]

# ── Main render ────────────────────────────────────────────────────────────────

def render(isotherms_dir: str) -> None:

    # ── Guard: calibration must be set ────────────────────────────────────────
    if "calibration_model" not in st.session_state:
        st.warning("⚠️ No calibration set yet. "
                   "Go to the **Calibration** tab and click **Set calibration**.")
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

    # ── File selector ──────────────────────────────────────────────────────────
    st.subheader("📂 Input Files — Isotherms")

    source = st.radio("File source:", ["📁 Folder", "⬆️ Upload"],
                      horizontal=True, key="iso_source")

    iso_entries = []  # list of {"name": str, "df": DataFrame}

    if source == "📁 Folder":
        all_iso = list_files(isotherms_dir, r"AIN\d")
        if not all_iso:
            st.warning(f"No AI* files found in **{isotherms_dir}**.")
            return
        saved_iso = st.session_state.get("iso_selected_files", all_iso)
        iso_sel = st.multiselect(
            "Choose Isotherm (AI) files:",
            options=all_iso,
            default=[f for f in saved_iso if f in all_iso] or all_iso,
            key="iso_file_selector",
        )
        st.session_state["iso_selected_files"] = iso_sel
        if not iso_sel:
            st.info("Select at least one isotherm file.")
            return
        for fname in iso_sel:
            df = read_spectral_file(os.path.join(isotherms_dir, fname))
            iso_entries.append({"name": fname, "df": df})
    else:
        uploaded_iso = st.file_uploader(
            "Upload AIN isotherm files (.txt)",
            type="txt", accept_multiple_files=True, key="iso_uploader",
        )
        if not uploaded_iso:
            st.info("Upload at least one AIN isotherm file.")
            return
        st.session_state["iso_selected_files"] = [f.name for f in uploaded_iso]
        for f in uploaded_iso:
            df = read_spectral_bytes(f)
            iso_entries.append({"name": f.name, "df": df})

    # ── Optional Co correction (from Control CTRL...C000... samples) ──────────
    df_co = st.session_state.get("co_correction_df")
    co_available = df_co is not None and not df_co.empty
    apply_corr = False
    with st.expander("🧪 Optional: Co correction (initial concentration)", expanded=False):
        st.markdown(
            "If a special **CTRL ...C000...** control sample was processed in "
            "the **Control** tab, you can replace the nominal N used in the "
            "qe and %R calculations with its measured value. The concentration "
            "matching the control sample is replaced directly; the other "
            "nominal concentrations are scaled by the same percentage "
            "variation found in that control."
        )
        if co_available:
            apply_corr = st.checkbox(
                "Apply Co correction to N (nominal concentration)",
                value=False, key="iso_apply_co_corr",
            )
        else:
            st.caption(
                "No Co-correction data found. Process a CTRL...C000... file "
                "in the **Control** tab to enable this option."
            )

    # ── Parse & compute (Ce, qe) for every file ───────────────────────────────
    raw_rows = []
    skipped  = []
    for e in iso_entries:
        meta = parse_filename(e["name"])
        if meta is None:
            skipped.append(
                f"{e['name']} (cannot parse — expected AIN<N><AdsCode><Cads>R<R>.txt)"
            )
            continue
        c508 = count_at_wavelength(e["df"], target=_peak_wl)
        if c508 is None:
            skipped.append(f"{e['name']} (no data near {_peak_wl:.0f} nm)")
            continue
        Ce = invert_calibration(cal, c508, crange)

        N_used = meta["N"]
        if apply_corr and co_available:
            corrected = apply_co_correction(df_co, meta["pH"], meta["N"])
            if corrected is not None:
                N_used = corrected

        qe = (N_used - Ce) / meta["Cads"]
        raw_rows.append({
            "Filename":        e["name"],
            "pH":              meta["pH"],
            "N (mg/L)":        meta["N"],
            "Ads":             meta["Ads"],
            "Cads (g/L)":      meta["Cads"],
            "R":               meta["R"],
            f"Counts @ {_peak_wl:.0f} nm": round(c508, 2),
            "Cₑ (mg/L)":       round(Ce, 4),
            "N used (mg/L)":   round(N_used, 4),
            "%R":              round((N_used - Ce) / N_used * 100, 4) if N_used > 1e-12 else None,
            "qₑ (mg/g)":       round(qe, 6),
        })

    for s in skipped:
        st.warning(f"Skipped: {s}")

    if not raw_rows:
        st.error("No valid isotherm data could be processed.")
        return

    raw_df = (
        pd.DataFrame(raw_rows)
        .sort_values(["Ads", "Cads (g/L)", "N (mg/L)", "R"])
        .reset_index(drop=True)
    )

    with st.expander("📋 Raw per-file data", expanded=False):
        st.dataframe(raw_df, width='stretch')

    # ── Detect whether multiple pH values are present ──────────────────────────
    ph_values = set(raw_df["pH"].dropna().unique())
    multi_ph  = len(ph_values) > 1

    if multi_ph:
        st.info(f"🧪 Detected **{len(ph_values)} pH values**: "
                f"{sorted(ph_values)}. Data will also be grouped by pH.")

    # Grouping columns depend on whether multiple pH values are present
    grp_avg    = (["pH", "N (mg/L)", "Ads", "Cads (g/L)"] if multi_ph
                  else       ["N (mg/L)", "Ads", "Cads (g/L)"])
    grp_series = (["pH", "Ads", "Cads (g/L)"] if multi_ph
                  else       ["Ads", "Cads (g/L)"])
    sort_avg   = (["Ads", "Cads (g/L)", "pH", "N (mg/L)"] if multi_ph
                  else ["Ads", "Cads (g/L)", "N (mg/L)"])

    # ── Average replicates ─────────────────────────────────────────────────────
    avg_df = (
        raw_df.groupby(grp_avg, as_index=False)
        .agg(
            Ce_mean=("Cₑ (mg/L)",  "mean"),
            Ce_std =("Cₑ (mg/L)",  "std"),
            qe_mean=("qₑ (mg/g)",  "mean"),
            qe_std =("qₑ (mg/g)",  "std"),
            Reps   =("R",          "count"),
            R_mean =("%R",         "mean"),
            R_std  =("%R",         "std"),
        )
        .rename(columns={"Ce_mean": "Cₑ (mg/L)", "qe_mean": "qₑ (mg/g)", "R_mean": "%R", "R_std": "%R std"})
        .sort_values(sort_avg)
        .reset_index(drop=True)
    )

    lbl_avg = ("mean per pH × N × Ads × Cads group" if multi_ph
               else "mean per N × Ads × Cads group")
    with st.expander(f"📋 Averaged data ({lbl_avg})", expanded=True):
        col_order_iso = (
            ["pH", "N (mg/L)", "Ads", "Cads (g/L)", "Reps", "Cₑ (mg/L)", "Ce_std", "%R", "%R std", "qₑ (mg/g)", "qe_std"]
            if multi_ph else
            ["N (mg/L)", "Ads", "Cads (g/L)", "Reps", "Cₑ (mg/L)", "Ce_std", "%R", "%R std", "qₑ (mg/g)", "qe_std"]
        )
        st.dataframe(avg_df.round(5), width='stretch', column_order=col_order_iso)

    # ── group_keys: una serie = un (pH ×) Ads × Cads ──────────────────────────
    group_keys = (
        avg_df[grp_series]
        .drop_duplicates()
        .sort_values(grp_series)
        .reset_index(drop=True)
    ) 


    st.divider()
    st.subheader("Isotherm Model Fitting  —  per (Adsorbent × C_ads) series")

    mc1, mc2 = st.columns(2)
    with mc1:
        st.markdown("**Langmuir:**")
        st.latex(r"q_e = \frac{q_{\max}\,K_L\,C_e}{1 + K_L\,C_e}, \qquad R_L=\frac{1}{1+K_L C_0}")
    with mc2:
        st.markdown("**Freundlich:**")
        st.latex(r"q_e = K_F\,C_e^{1/n}")

    fig_iso   = go.Figure()   # all models on one plot
    fig_lin_L = go.Figure()   # Langmuir linear: Ce/qe vs Ce
    fig_lin_F = go.Figure()   # Freundlich linear: ln(qe) vs ln(Ce)

    all_series_results = []

    for idx, grow in group_keys.iterrows():
        ads_code = grow["Ads"]
        cads_val = grow["Cads (g/L)"]
        ph_val   = grow["pH"] if multi_ph else None   # ← añadir

        mask_series = (
            (avg_df["Ads"] == ads_code) &
            (avg_df["Cads (g/L)"] == cads_val)
        )
        if multi_ph:
            mask_series &= (avg_df["pH"] == ph_val)   # ← añadir filtro pH

        sub = avg_df[mask_series].sort_values("Cₑ (mg/L)")

        C0_ref = float(sub["N (mg/L)"].max())
        Ce    = sub["Cₑ (mg/L)"].values.astype(float)
        qe    = sub["qₑ (mg/g)"].values.astype(float)
        q_std = sub["qe_std"].fillna(0).values.astype(float)

        # keep only positive values for fitting / log transforms
        mask  = (Ce > 0) & (qe > 0)
        Ce_m, qe_m = Ce[mask], qe[mask]

        color = _color(idx)
        label = (f"pH={ph_val:.4g} | {ads_code}  C_ads={cads_val:.4g} g/L"
                 if multi_ph else
                 f"{ads_code}  C_ads={cads_val:.4g} g/L")        

        # ── Fits ──────────────────────────────────────────────────────────────
        qmax_guess = max(qe_m.max() * 1.2, 1e-6) if len(qe_m) else 25.0

        popt_L, r2_L, mae_L = _fit(
            langmuir, Ce_m, qe_m,
            p0=[qmax_guess, 1.0],
            bounds=([0, 0], [np.inf, np.inf]),
        )
        popt_F, r2_F, mae_F = _fit(
            freundlich, Ce_m, qe_m,
            p0=[10.0, 2.0],
            bounds=([0, 1e-6], [np.inf, np.inf]),
        )

        Ce_fit = np.linspace(max(Ce_m.min() * 0.5, 1e-3), Ce_m.max() * 1.1, 400) \
                 if len(Ce_m) else np.array([])

        # scatter data with error bars
        has_err = bool((q_std > 0).any())
        fig_iso.add_trace(go.Scatter(
            x=Ce, y=qe,
            error_y=dict(type="data", array=q_std, visible=has_err),
            mode="markers",
            name=f"{label} — Data",
            marker=dict(color=color, size=9, symbol="circle"),
        ))
        if popt_L is not None and len(Ce_fit):
            fig_iso.add_trace(go.Scatter(
                x=Ce_fit, y=langmuir(Ce_fit, *popt_L), mode="lines",
                name=f"{label} — Langmuir",
                line=dict(color=color, width=2),
            ))
        if popt_F is not None and len(Ce_fit):
            fig_iso.add_trace(go.Scatter(
                x=Ce_fit, y=freundlich(Ce_fit, *popt_F), mode="lines",
                name=f"{label} — Freundlich",
                line=dict(color=color, width=2, dash="dash"),
            ))

        # Langmuir linearised: Ce/qe vs Ce
        if len(Ce_m) >= 2:
            y_LL  = Ce_m / qe_m
            reg_L = LinearRegression().fit(Ce_m.reshape(-1, 1), y_LL)
            sl_L, ic_L = float(reg_L.coef_[0]), float(reg_L.intercept_)
            Ce_r  = np.linspace(Ce_m.min(), Ce_m.max(), 200)
            fig_lin_L.add_trace(go.Scatter(x=Ce_m, y=y_LL, mode="markers",
                                           name=label, marker=dict(color=color, size=8)))
            fig_lin_L.add_trace(go.Scatter(x=Ce_r, y=ic_L + sl_L * Ce_r, mode="lines",
                                           line=dict(color=color, width=2), showlegend=False))

        # Freundlich linearised: ln(qe) vs ln(Ce)
        if len(Ce_m) >= 2:
            log_Ce = np.log(Ce_m)
            log_qe = np.log(qe_m)
            reg_F  = LinearRegression().fit(log_Ce.reshape(-1, 1), log_qe)
            sl_F, ic_F = float(reg_F.coef_[0]), float(reg_F.intercept_)
            lCe_r  = np.linspace(log_Ce.min(), log_Ce.max(), 200)
            fig_lin_F.add_trace(go.Scatter(x=log_Ce, y=log_qe, mode="markers",
                                           name=label, marker=dict(color=color, size=8)))
            fig_lin_F.add_trace(go.Scatter(x=lCe_r, y=ic_F + sl_F * lCe_r, mode="lines",
                                           line=dict(color=color, width=2), showlegend=False))

        # ── Per-series parameter panel ─────────────────────────────────────────
        n_pts = len(Ce_m)
        n_reps_info = sub["Reps"].sum() if "Reps" in sub.columns else "?"
        with st.expander(
            f"📊  Series: **{label}**  ({n_pts} data points, {int(n_reps_info)} total replicas)",
            expanded=True,
        ):
            pc1, pc2 = st.columns(2)
            with pc1:
                st.markdown("**Langmuir:**")
                if popt_L is not None:
                    qmax_f, KL_f = popt_L
                    RL = 1.0 / (1.0 + KL_f * C0_ref)
                    st.markdown(f"- q_max = **{qmax_f:.5f}** mg/g")
                    st.markdown(f"- K_L   = **{KL_f:.6f}** L/mg")
                    st.markdown(f"- R²    = **{r2_L:.5f}**")
                    st.markdown(f"- MAE   = **{mae_L:.5f}** mg/g")
                    st.markdown(f"- **R_L = {RL:.5f}**")
                    if RL > 1:
                        st.error("R_L > 1 → unfavorable adsorption")
                    elif abs(RL - 1) < 1e-4:
                        st.warning("R_L ≈ 1 → linear")
                    elif RL > 0:
                        st.success("0 < R_L < 1 → ✅ favorable adsorption")
                    else:
                        st.info("R_L ≈ 0 → irreversible")
                else:
                    st.warning("Langmuir fit did not converge.")
            with pc2:
                st.markdown("**Freundlich:**")
                if popt_F is not None:
                    KF_f, n_f = popt_F
                    inv_n = 1.0 / n_f
                    st.markdown(f"- K_F = **{KF_f:.5f}**  (mg/g)·(L/mg)^(1/n)")
                    st.markdown(f"- n   = **{n_f:.5f}**")
                    st.markdown(f"- 1/n = **{inv_n:.5f}**")
                    st.markdown(f"- R²  = **{r2_F:.5f}**")
                    st.markdown(f"- MAE = **{mae_F:.5f}** mg/g")
                    if inv_n < 1:
                        st.success("1/n < 1 → ✅ favorable adsorption")
                    elif abs(inv_n - 1) < 1e-3:
                        st.warning("1/n ≈ 1 → linear")
                    else:
                        st.error("1/n > 1 → unfavorable adsorption")
                else:
                    st.warning("Freundlich fit did not converge.")

        # collect for download
        res = {
            "label":    label,
            "pH":       ph_val,
            "Ads":      ads_code,
            "Cads":     cads_val,
            "Ce":       Ce_m,
            "qe":       qe_m,
            "C0_ref":   C0_ref,
        }
        if popt_L is not None:
            qmax_f, KL_f = popt_L
            res["langmuir"] = {
                "qmax": qmax_f, "KL": KL_f,
                "RL": 1.0 / (1.0 + KL_f * C0_ref),
                "r2": r2_L, "mae": mae_L,
            }
        if popt_F is not None:
            KF_f, n_f = popt_F
            res["freundlich"] = {
                "KF": KF_f, "n": n_f, "inv_n": 1.0 / n_f,
                "r2": r2_F, "mae": mae_F,
            }
        all_series_results.append(res)

    # ── Combined plots ─────────────────────────────────────────────────────────
    st.divider()
    st.subheader("Isotherm curves (all series)")
    fig_iso.update_layout(
        xaxis_title="Cₑ (mg/L)", yaxis_title="qₑ (mg/g)",
        template="plotly_white", height=480,
    )
    st.plotly_chart(fig_iso, width='stretch')

    st.subheader("Linearised Forms")
    lc1, lc2 = st.columns(2)
    with lc1:
        fig_lin_L.update_layout(
            title="Langmuir: Cₑ/qₑ vs Cₑ",
            xaxis_title="Cₑ (mg/L)", yaxis_title="Cₑ / qₑ  (g/L)",
            template="plotly_white", height=360,
        )
        st.plotly_chart(fig_lin_L, width='stretch')
    with lc2:
        fig_lin_F.update_layout(
            title="Freundlich: ln(qₑ) vs ln(Cₑ)",
            xaxis_title="ln(Cₑ)", yaxis_title="ln(qₑ)",
            template="plotly_white", height=360,
        )
        st.plotly_chart(fig_lin_F, width='stretch')


    # ══════════════════════════════════════════════════════════════════════════
    # Extra isotherm models — user-selectable
    # ══════════════════════════════════════════════════════════════════════════

    # ── Local model & helper functions ────────────────────────────────────────

    def _sips(Ce, qm, KS, ne):
        Ce = np.asarray(Ce, float)
        return qm * (KS * Ce) ** ne / (1.0 + (KS * Ce) ** ne)

    def _toth(Ce, qm, KT, tp):
        Ce = np.asarray(Ce, float)
        return qm * KT * Ce / (1.0 + (KT * Ce) ** tp) ** (1.0 / tp)

    def _halsey(Ce, KH, nH):
        Ce = np.asarray(Ce, float)
        inner = np.log(1.0 / np.maximum(Ce, 1e-12))
        return np.where(inner > 1e-12,
                        (KH / np.maximum(inner, 1e-12)) ** (1.0 / nH),
                        np.nan)

    def _harkins_jura(Ce, A, B):
        Ce = np.asarray(Ce, float)
        inner = B / A - np.log10(np.maximum(Ce, 1e-12)) / A
        return np.where(inner > 1e-12,
                        1.0 / np.sqrt(np.maximum(inner, 1e-12)),
                        np.nan)

    def _janovics(Ce, qmax, KJ):
        Ce = np.asarray(Ce, float)
        return qmax * (1.0 - np.exp(-KJ * Ce))

    def _fit_extra_iso(fn, Ce, qe, p0, bounds):
        from scipy.optimize import curve_fit
        from sklearn.metrics import r2_score, mean_absolute_error as _mae
        try:
            # Pre-check: discard points where the model returns NaN at p0
            trial = np.asarray(fn(Ce, *p0), float)
            valid = np.isfinite(trial) & np.isfinite(qe)
            if valid.sum() < 2:
                return None, float("nan"), float("nan")
            Ce_v, qe_v = Ce[valid], qe[valid]
            popt, _ = curve_fit(fn, Ce_v, qe_v, p0=p0, bounds=bounds, maxfev=30_000)
            pred = fn(Ce_v, *popt)
            return popt, float(r2_score(qe_v, pred)), float(_mae(qe_v, pred))
        except Exception:
            return None, float("nan"), float("nan")

    # ── Model registry ────────────────────────────────────────────────────────
    _EXTRA_ISO_MODELS = {
        "Sips (Langmuir-Freundlich)": {
            "latex":        r"q_e = \frac{q_m(K_S C_e)^{n_e}}{1+(K_S C_e)^{n_e}}",
            "params_label": "qₘ, Kₛ, nₑ",
        },
        "Tóth": {
            "latex":        r"q_e = \frac{q_m K_T C_e}{\bigl[1+(K_T C_e)^t\bigr]^{1/t}}",
            "params_label": "qₘ, K_T, t",
        },
        "Halsey": {
            "latex":        r"q_e = \!\left(\frac{K_H}{\ln(1/C_e)}\right)^{\!1/n_H}",
            "params_label": "K_H, n_H",
            "note":         "Only valid for Cₑ < 1 mg/L; points out of range are ommited.",
        },
        "Harkins-Jura": {
            "latex":        r"\frac{1}{q_e^2}=\frac{B}{A}-\frac{1}{A}\log C_e",
            "params_label": "A, B",
        },
        "Janovics": {
            "latex":        r"q_e = q_{max}\!\left(1-e^{-K_J C_e}\right)",
            "params_label": "qₘₐₓ, K_J",
        },
    }

    # ── UI: heading + equation reference ──────────────────────────────────────
    st.divider()
    st.subheader("🔬 Additional Isotherms models")
    st.markdown(
        "**Langmuir** and **Freundlich** models are the default (parameters obtained in the previous section). "
        "Select here any additional model to fit:"
    )

    with st.expander("📋 Available additional models", expanded=False):
        for _mn, _mi in _EXTRA_ISO_MODELS.items():
            _ca, _cb = st.columns([1, 2])
            _ca.markdown(f"**{_mn}**")
            _ca.caption(f"parameters: {_mi['params_label']}")
            if "note" in _mi:
                _ca.caption(f"⚠️ {_mi['note']}")
            _cb.latex(_mi["latex"])

    selected_extra_iso = st.multiselect(
        "Additional models to fit:",
        options=list(_EXTRA_ISO_MODELS.keys()),
        default=st.session_state.get("_iso_extra_sel", []),
        key="_iso_extra_sel",
    )

    extra_iso_results = []   # persisted to session_state["iso_results"]["extra_models"]

    if selected_extra_iso:
        st.markdown("---")
        for _gi, _gr in group_keys.iterrows():
            _ads  = _gr["Ads"]
            _cads = _gr["Cads (g/L)"]
            _ph   = _gr.get("pH") if multi_ph else None

            _mask = (
                (avg_df["Ads"]        == _ads)  &
                (avg_df["Cads (g/L)"] == _cads)
            )
            if multi_ph:
                _mask &= (avg_df["pH"] == _ph)

            _sub  = avg_df[_mask].sort_values("Cₑ (mg/L)")
            _Ce   = _sub["Cₑ (mg/L)"].values.astype(float)
            _qe   = _sub["qₑ (mg/g)"].values.astype(float)
            _vmask = (_Ce > 0) & (_qe > 0)
            _Ce_m, _qe_m = _Ce[_vmask], _qe[_vmask]

            _label = (f"pH={_ph:.4g} | {_ads}  C_ads={_cads:.4g} g/L"
                      if multi_ph else f"{_ads}  C_ads={_cads:.4g} g/L")

            color_val   = _color(_gi)   # isotherms.py already defines _color()
            _qm_ref  = max(_qe_m.max() * 1.2, 1e-6) if len(_qe_m) else 25.0
            _Ce_plot = (np.linspace(max(_Ce_m.min() * 0.5, 1e-3),
                                    _Ce_m.max() * 1.1, 400)
                        if len(_Ce_m) else np.array([]))

            with st.expander(f"📊  Series: **{_label}**", expanded=True):
                _fig_ex = go.Figure()
                _fig_ex.add_trace(go.Scatter(
                    x=_Ce_m, y=_qe_m, mode="markers", name="Datos",
                    marker=dict(color=color_val, size=9, symbol="circle"),
                ))

                for _mn in selected_extra_iso:
                    _rec = {
                        "model": _mn,
                        "label": _label,
                        "Ads":   _ads,
                        "Cads":  _cads,
                    }

                    # ── Sips ──────────────────────────────────────────────────
                    if _mn == "Sips (Langmuir-Freundlich)":
                        _popt, _r2, _mae = _fit_extra_iso(
                            _sips, _Ce_m, _qe_m,
                            p0=[_qm_ref, 1.0, 1.0],
                            bounds=([0.0, 1e-9, 0.01], [np.inf, np.inf, 10.0]),
                        )
                        if _popt is not None:
                            _qm, _KS, _ne = _popt
                            _rec["params"] = {"qm": _qm, "KS": _KS, "ne": _ne}
                            _rec["r2"], _rec["mae"] = _r2, _mae
                            if len(_Ce_plot):
                                _fig_ex.add_trace(go.Scatter(
                                    x=_Ce_plot, y=_sips(_Ce_plot, _qm, _KS, _ne),
                                    mode="lines", name=_mn, line=dict(width=2),
                                ))
                            st.markdown(
                                f"**{_mn}** — qₘ=`{_qm:.5f}` mg/g | "
                                f"Kₛ=`{_KS:.6f}` L/mg | nₑ=`{_ne:.4f}` | "
                                f"R²=`{_r2:.5f}` | MAE=`{_mae:.5f}`"
                            )
                        else:
                            st.warning(f"**{_mn}**: fit did not converged.")
                            _rec["params"] = {}
                            _rec["r2"] = _rec["mae"] = float("nan")

                    # ── Tóth ──────────────────────────────────────────────────
                    elif _mn == "Tóth":
                        _popt, _r2, _mae = _fit_extra_iso(
                            _toth, _Ce_m, _qe_m,
                            p0=[_qm_ref, 1.0, 1.0],
                            bounds=([0.0, 1e-9, 0.01], [np.inf, np.inf, 20.0]),
                        )
                        if _popt is not None:
                            _qm, _KT, _tp = _popt
                            _rec["params"] = {"qm": _qm, "KT": _KT, "t": _tp}
                            _rec["r2"], _rec["mae"] = _r2, _mae
                            if len(_Ce_plot):
                                _fig_ex.add_trace(go.Scatter(
                                    x=_Ce_plot, y=_toth(_Ce_plot, _qm, _KT, _tp),
                                    mode="lines", name=_mn,
                                    line=dict(width=2, dash="dot"),
                                ))
                            st.markdown(
                                f"**{_mn}** — qₘ=`{_qm:.5f}` mg/g | "
                                f"K_T=`{_KT:.6f}` | t=`{_tp:.4f}` | "
                                f"R²=`{_r2:.5f}` | MAE=`{_mae:.5f}`"
                            )
                        else:
                            st.warning(f"**{_mn}**: fit did not converged.")
                            _rec["params"] = {}
                            _rec["r2"] = _rec["mae"] = float("nan")

                    # ── Halsey ─────────────────────────────────────────────────
                    elif _mn == "Halsey":
                        _h_mask = (_Ce_m < 1.0) & (_Ce_m > 0) & (_qe_m > 0)
                        if _h_mask.sum() < 2:
                            st.warning(
                                f"**{_mn}**: not enough points with Cₑ < 1 mg/L "
                                f"(found: {int(_h_mask.sum())})."
                            )
                            _rec["params"] = {}
                            _rec["r2"] = _rec["mae"] = float("nan")
                        else:
                            _popt, _r2, _mae = _fit_extra_iso(
                                _halsey, _Ce_m[_h_mask], _qe_m[_h_mask],
                                p0=[1.0, 2.0],
                                bounds=([1e-9, 0.01], [np.inf, 20.0]),
                            )
                            if _popt is not None:
                                _KH, _nH = _popt
                                _rec["params"] = {"KH": _KH, "nH": _nH}
                                _rec["r2"], _rec["mae"] = _r2, _mae
                                _Ce_h = np.linspace(
                                    1e-3, min(_Ce_m[_h_mask].max() * 1.1, 0.99), 200)
                                _fig_ex.add_trace(go.Scatter(
                                    x=_Ce_h, y=_halsey(_Ce_h, _KH, _nH),
                                    mode="lines", name=f"{_mn} (Cₑ<1)",
                                    line=dict(width=2, dash="dashdot"),
                                ))
                                st.markdown(
                                    f"**{_mn}** (fit for Cₑ<1 mg/L) — "
                                    f"K_H=`{_KH:.6f}` | n_H=`{_nH:.4f}` | "
                                    f"R²=`{_r2:.5f}` | MAE=`{_mae:.5f}`"
                                )
                            else:
                                st.warning(f"**{_mn}**: fit did not converged.")
                                _rec["params"] = {}
                                _rec["r2"] = _rec["mae"] = float("nan")

                    # ── Harkins-Jura ───────────────────────────────────────────
                    elif _mn == "Harkins-Jura":
                        _popt, _r2, _mae = _fit_extra_iso(
                            _harkins_jura, _Ce_m, _qe_m,
                            p0=[10.0, 1.0],
                            bounds=([1e-9, -np.inf], [np.inf, np.inf]),
                        )
                        if _popt is not None:
                            _A, _B = _popt
                            _rec["params"] = {"A": _A, "B": _B}
                            _rec["r2"], _rec["mae"] = _r2, _mae
                            if len(_Ce_plot):
                                _y_hj = _harkins_jura(_Ce_plot, _A, _B)
                                _vhj  = np.isfinite(_y_hj)
                                if _vhj.any():
                                    _fig_ex.add_trace(go.Scatter(
                                        x=_Ce_plot[_vhj], y=_y_hj[_vhj],
                                        mode="lines", name=_mn,
                                        line=dict(width=2, dash="longdash"),
                                    ))
                            st.markdown(
                                f"**{_mn}** — A=`{_A:.6f}` | B=`{_B:.6f}` | "
                                f"R²=`{_r2:.5f}` | MAE=`{_mae:.5f}`"
                            )
                        else:
                            st.warning(f"**{_mn}**: fit did not converged.")
                            _rec["params"] = {}
                            _rec["r2"] = _rec["mae"] = float("nan")

                    # ── Janovics ───────────────────────────────────────────────
                    elif _mn == "Janovics":
                        _popt, _r2, _mae = _fit_extra_iso(
                            _janovics, _Ce_m, _qe_m,
                            p0=[_qm_ref, 0.1],
                            bounds=([0.0, 1e-9], [np.inf, np.inf]),
                        )
                        if _popt is not None:
                            _qmax, _KJ = _popt
                            _rec["params"] = {"qmax": _qmax, "KJ": _KJ}
                            _rec["r2"], _rec["mae"] = _r2, _mae
                            if len(_Ce_plot):
                                _fig_ex.add_trace(go.Scatter(
                                    x=_Ce_plot, y=_janovics(_Ce_plot, _qmax, _KJ),
                                    mode="lines", name=_mn,
                                    line=dict(width=2),
                                ))
                            st.markdown(
                                f"**{_mn}** — qₘₐₓ=`{_qmax:.5f}` mg/g | "
                                f"K_J=`{_KJ:.6f}` L/mg | "
                                f"R²=`{_r2:.5f}` | MAE=`{_mae:.5f}`"
                            )
                        else:
                            st.warning(f"**{_mn}**: fit did not converged.")
                            _rec["params"] = {}
                            _rec["r2"] = _rec["mae"] = float("nan")

                    extra_iso_results.append(_rec)

                # ── Combined plot ──────────────────────────────────────────────
                _fig_ex.update_layout(
                    xaxis_title="Cₑ  (mg/L)", yaxis_title="qₑ  (mg/g)",
                    template="plotly_white", height=400,
                    legend=dict(orientation="h", yanchor="bottom", y=1.02,
                                font=dict(size=11)),
                )
                st.plotly_chart(_fig_ex, width='stretch')



    # ── Persist for Download tab ───────────────────────────────────────────────

    iso_lin_L_rows, iso_lin_F_rows = [], []
    for s in all_series_results:
        Ce_m = s["Ce"]
        qe_m = s["qe"]
        lbl  = s["label"]
        for c, q in zip(Ce_m, qe_m):
            iso_lin_L_rows.append({"Series": lbl, "Ce_mg_per_L": c,
                                   "Ce_over_qe": c / q if q > 1e-12 else None})
            if c > 0 and q > 0:
                iso_lin_F_rows.append({"Series": lbl,
                                       "ln_Ce": np.log(c), "ln_qe": np.log(q)})

    st.session_state["iso_results"] = {
        "co_correction_applied": apply_corr,
        "files":  [e["name"] for e in iso_entries],
        "raw_df": raw_df,
        "avg_df": avg_df,
        "series": all_series_results,
        "csv_iso":       avg_df[
            (["pH", "Ads", "Cads (g/L)", "N (mg/L)", "Reps", "Cₑ (mg/L)", "Ce_std", "%R", "%R std", "qₑ (mg/g)", "qe_std"]
             if multi_ph else
             ["Ads", "Cads (g/L)", "N (mg/L)", "Reps", "Cₑ (mg/L)", "Ce_std", "%R", "%R std", "qₑ (mg/g)", "qe_std"])
        ].copy(),
        "csv_lin_lang":  pd.DataFrame(iso_lin_L_rows),
        "csv_lin_freun": pd.DataFrame(iso_lin_F_rows),
        "extra_models":  extra_iso_results,
    }
