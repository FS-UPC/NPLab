"""
tabs/kinetics.py — Kinetics tab.

Filename convention: AKN{N}{AdsCode}{AdsConc}t{time}R{rep}.txt
  e.g.  AKN50BTO2t05R1.txt  →  N=50, AdsCode=BTO, AdsConc=2, t=5 min, R=1

Grouping logic
──────────────
  If the kinetics folder path contains the subfolder "M07":
    Series key : (N, AdsCode, AdsConc, Parity)
                 Parity = "odd"  for R1, R3, R5, …
                 Parity = "even" for R2, R4, R6, …

  Otherwise (all other folders):
    Series key : (N, AdsCode, AdsConc)
                 All replicas are averaged together.

  Per series  : group data points by t, average q(t) within the group.
  Fit         : pseudo-1st and 2nd order models fitted independently per series.
"""

import os
import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from utils import (
    list_files,
    read_spectral_file,
    read_spectral_bytes,
    parse_filename,
    count_at_wavelength,
    invert_calibration,
    apply_co_correction,
    model_first_order,
    model_second_order,
    fit_kinetic,
    linear_first_order,
    linear_second_order,
)

# ── Colour palette ─────────────────────────────────────────────────────────────

_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]

def series_color(idx: int) -> str:
    return _PALETTE[idx % len(_PALETTE)]

# ── Parity helper ──────────────────────────────────────────────────────────────

def _use_parity(kinetics_dir: str) -> bool:
    """Return True if any component of the path is exactly 'M07' (case-insensitive)."""
    parts = re.split(r"[/\\]", os.path.normcase(kinetics_dir))
    return "M07" in parts

# ── Main render ────────────────────────────────────────────────────────────────

def render(kinetics_dir: str) -> None:

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

    use_parity = _use_parity(kinetics_dir)
    if use_parity:
        st.info("ℹ️ **M07 mode**: replicates are split by parity (odd/even).")

    st.divider()

    # ── File selector ──────────────────────────────────────────────────────────
    st.subheader("📂 Input Files — Kinetics")

    source_kin = st.radio("File source:", ["📁 Folder", "⬆️ Upload"],
                          horizontal=True, key="kin_source")

    kin_entries = []  # list of {"name": str, "df": DataFrame}

    if source_kin == "📁 Folder":
        all_kin = list_files(kinetics_dir, r"AKN\d")
        if not all_kin:
            st.warning(f"No AK* files found in **{kinetics_dir}**.")
            return
        saved_kin = st.session_state.get("kin_selected_files", all_kin)
        kin_sel = st.multiselect(
            "Choose Kinetics (AK) files:",
            options=all_kin,
            default=[f for f in saved_kin if f in all_kin] or all_kin,
            key="kin_file_selector",
        )
        st.session_state["kin_selected_files"] = kin_sel
        if not kin_sel:
            st.info("Select at least one kinetics file.")
            return
        for fname in kin_sel:
            df = read_spectral_file(os.path.join(kinetics_dir, fname))
            kin_entries.append({"name": fname, "df": df})
    else:
        uploaded_kin = st.file_uploader(
            "Upload AKN kinetics files (.txt)",
            type="txt", accept_multiple_files=True, key="kin_uploader",
        )
        if not uploaded_kin:
            st.info("Upload at least one AKN kinetics file.")
            return
        st.session_state["kin_selected_files"] = [f.name for f in uploaded_kin]
        for f in uploaded_kin:
            df = read_spectral_bytes(f)
            kin_entries.append({"name": f.name, "df": df})

    # ── Optional Co correction (from Control CTRL...C000... samples) ──────────
    df_co = st.session_state.get("co_correction_df")
    co_available = df_co is not None and not df_co.empty
    apply_corr = False
    with st.expander("🧪 Optional: Co correction (initial concentration)", expanded=False):
        st.markdown(
            "If a special **CTRL ...C000...** control sample was processed in "
            "the **Control** tab, you can replace the nominal N used in the "
            "qe, %R and q(t) calculations with its measured value."
        )
        if co_available:
            apply_corr = st.checkbox(
                "Apply Co correction to N (nominal concentration)",
                value=False, key="kin_apply_co_corr",
            )
        else:
            st.caption(
                "No Co-correction data found. Process a CTRL...C000... file "
                "in the **Control** tab to enable this option."
            )

    # ── Parse & compute q(t) for every file ───────────────────────────────────
    raw_rows = []
    skipped  = []
    for e in kin_entries:
        meta = parse_filename(e["name"])
        if meta is None:
            skipped.append(f"{e['name']} (cannot parse filename)")
            continue
        c508 = count_at_wavelength(e["df"], target=_peak_wl)
        if c508 is None:
            skipped.append(f"{e['name']} (no data near {_peak_wl:.0f} nm)")
            continue
        C_t = invert_calibration(cal, c508, crange)

        N_used = meta["N"]
        if apply_corr and co_available:
            corrected = apply_co_correction(df_co, meta["pH"], meta["N"])
            if corrected is not None:
                N_used = corrected

        q_t = (N_used - C_t) / meta["AdsConc"]
        raw_rows.append({
            "Filename":        e["name"],
            "pH":              meta["pH"],
            "N":               meta["N"],
            "AdsCode":         meta["AdsCode"],
            "AdsConc":         meta["AdsConc"],
            "t (min)":         meta["t"],
            "R":               meta["R"],
            f"Counts @ {_peak_wl:.0f} nm": round(c508, 2),
            "Cₜ (mg/L)":       round(C_t, 4),
            "N used (mg/L)":   round(N_used, 4),
            "%R":              round((N_used - C_t) / N_used * 100, 4) if N_used > 1e-12 else None,
            "q(t) (mg/g)":     round(q_t, 6),
        })

    for s in skipped:
        st.warning(f"Skipped: {s}")

    if not raw_rows:
        st.error("No valid kinetics data could be processed.")
        return

    raw_df = pd.DataFrame(raw_rows).sort_values(
        ["N", "AdsCode", "AdsConc", "t (min)", "R"]
    ).reset_index(drop=True)

    # ── Parity column (only used when in M07 mode) ─────────────────────────────
    raw_df["Parity"] = raw_df["R"].apply(lambda r: "odd" if r % 2 != 0 else "even")

    # ── Raw data table ─────────────────────────────────────────────────────────
    with st.expander("📋 Raw per-file data", expanded=False):
        st.dataframe(raw_df, width='stretch')

    # ── Average & build series keys ────────────────────────────────────────────
    ph_values = set(raw_df["pH"].dropna().unique())
    multi_ph  = len(ph_values) > 1

    if multi_ph:
        st.info(f"🧪 Detected **{len(ph_values)} pH values**: {sorted(ph_values)}.")

    if use_parity:
        grp_cols = (["pH", "N", "AdsCode", "AdsConc", "Parity", "t (min)"]
                    if multi_ph else
                    ["N", "AdsCode", "AdsConc", "Parity", "t (min)"])
    else:
        grp_cols = (["pH", "N", "AdsCode", "AdsConc", "t (min)"]
                    if multi_ph else
                    ["N", "AdsCode", "AdsConc", "t (min)"])


    avg_df = (
        raw_df.groupby(grp_cols, as_index=False)
        .agg(
            **{
                "q (mg/g)":       ("q(t) (mg/g)", "mean"),
                "q std":          ("q(t) (mg/g)", "std"),
                "Reps":           ("q(t) (mg/g)", "count"),
                "Cₜ (mg/L)":     ("Cₜ (mg/L)",   "mean"),
                "Cₜ std (mg/L)": ("Cₜ (mg/L)",   "std"),
                "%R":             ("%R",           "mean"),
                "%R std":         ("%R",           "std"),
            }
        )
        .reset_index()
    )

    expander_label = (
        "📋 Averaged data (mean per parity group)"
        if use_parity else
        "📋 Averaged data (mean per series)"
    )
    with st.expander(expander_label, expanded=False):
        col_order = (
            ("N", "AdsCode", "AdsConc", "Parity", "t (min)", "Reps", "Cₜ (mg/L)", "Cₜ std (mg/L)", "%R", "%R std", "q (mg/g)", "q std")
            if use_parity else
            ("N", "AdsCode", "AdsConc", "t (min)", "Reps", "Cₜ (mg/L)", "Cₜ std (mg/L)", "%R", "%R std", "q (mg/g)", "q std")
        )
        st.dataframe(avg_df.round(5), width='stretch',
                     column_order=col_order)

    # ── Series keys ────────────────────────────────────────────────────────────
    if use_parity:
        key_cols = ["pH", "N", "AdsCode", "AdsConc", "Parity"] if multi_ph else ["N", "AdsCode", "AdsConc", "Parity"]
    else:
        key_cols = ["pH", "N", "AdsCode", "AdsConc"] if multi_ph else ["N", "AdsCode", "AdsConc"]

    series_keys = (
        avg_df[key_cols]
        .drop_duplicates()
        .sort_values(key_cols)
        .reset_index(drop=True)
    )

    st.divider()
    st.subheader("Kinetic Model Fitting  —  per series")

    tpc1, tpc2 = st.columns(2)
    with tpc1:
        st.markdown("**Pseudo-1st-order:**")
        st.latex(r"q(t)=q_e(1-e^{-k_1 t})")
    with tpc2:
        st.markdown("**Pseudo-2nd-order:**")
        st.latex(r"q(t)=\frac{k_2 q_e^2 t}{1+k_2 q_e t}")

    fig_qt   = go.Figure()
    fig_lin1 = go.Figure()
    fig_lin2 = go.Figure()

    all_series_results = []

    for idx, row in series_keys.iterrows():
        N_val  = row["N"]
        ads    = row["AdsCode"]
        adsc   = row["AdsConc"]
        parity = row["Parity"] if use_parity else None

        # Build boolean mask for this series
        mask = (
            (avg_df["N"]       == N_val) &
            (avg_df["AdsCode"] == ads)   &
            (avg_df["AdsConc"] == adsc)
        )
        if use_parity:
            mask &= (avg_df["Parity"] == parity)

        sub   = avg_df[mask].sort_values("t (min)")
        t_arr = sub["t (min)"].values.astype(float)
        q_arr = sub["q (mg/g)"].values.astype(float)
        q_std = sub["q std"].fillna(0).values.astype(float)

        color = series_color(idx)
        if use_parity:
            label = (f"pH={row['pH']} | N={N_val} | {ads}{adsc:.4g} | {parity} reps"
            if multi_ph else
            f"N={N_val} | {ads}{adsc:.4g} | {parity} reps")
        else:
            label = (f"pH={row['pH']} | N={N_val} | {ads}{adsc:.4g}"
            if multi_ph else
            f"N={N_val} | {ads}{adsc:.4g}")

        # fits
        qe_guess = max(q_arr.max() * 1.1, 1e-6)
        popt1, r2_1, mae1 = fit_kinetic(
            model_first_order,  t_arr, q_arr,
            p0=[qe_guess, 0.01], bounds=([0, 0], [np.inf, np.inf]),
        )
        popt2, r2_2, mae2 = fit_kinetic(
            model_second_order, t_arr, q_arr,
            p0=[qe_guess, 0.001], bounds=([0, 0], [np.inf, np.inf]),
        )

        t_fit = np.linspace(0, t_arr.max() * 1.05, 500)

        # experimental points
        has_err = bool((q_std > 0).any())
        fig_qt.add_trace(go.Scatter(
            x=t_arr, y=q_arr,
            error_y=dict(type="data", array=q_std, visible=has_err),
            mode="markers",
            name=f"{label} — Data",
            marker=dict(color=color, size=9, symbol="circle"),
        ))

        if popt1 is not None:
            qe1, k1 = popt1
            fig_qt.add_trace(go.Scatter(
                x=t_fit, y=model_first_order(t_fit, qe1, k1), mode="lines",
                name=f"{label} — Pseudo-1st-order",
                line=dict(color=color, width=2),
            ))
        if popt2 is not None:
            qe2, k2 = popt2
            fig_qt.add_trace(go.Scatter(
                x=t_fit, y=model_second_order(t_fit, qe2, k2), mode="lines",
                name=f"{label} — Pseudo-2nd-order",
                line=dict(color=color, width=2, dash="dash"),
            ))

        # linearised 1st order
        if popt1 is not None:
            qe1, k1 = popt1
            sl1, ic1, qe_safe = linear_first_order(t_arr, q_arr, qe1)
            y_lin1  = np.log(qe_safe - q_arr)
            t_range = np.linspace(t_arr.min(), t_arr.max(), 200)
            fig_lin1.add_trace(go.Scatter(
                x=t_arr, y=y_lin1, mode="markers",
                name=label, marker=dict(color=color, size=8),
            ))
            if sl1 is not None:
                fig_lin1.add_trace(go.Scatter(
                    x=t_range, y=ic1 + sl1 * t_range, mode="lines",
                    line=dict(color=color, width=2), showlegend=False,
                ))

        # linearised 2nd order
        if popt2 is not None:
            mask2 = q_arr > 1e-12
            t_m2  = t_arr[mask2]
            y_m2  = t_m2 / q_arr[mask2]
            sl2, ic2 = linear_second_order(t_arr, q_arr)
            t_range2 = np.linspace(t_m2.min(), t_m2.max(), 200)
            fig_lin2.add_trace(go.Scatter(
                x=t_m2, y=y_m2, mode="markers",
                name=label, marker=dict(color=color, size=8),
            ))
            if sl2 is not None:
                fig_lin2.add_trace(go.Scatter(
                    x=t_range2, y=ic2 + sl2 * t_range2, mode="lines",
                    line=dict(color=color, width=2), showlegend=False,
                ))

        # per-series parameter panel
        with st.expander(f"📊  Series: **{label}**", expanded=True):
            pc1, pc2 = st.columns(2)
            with pc1:
                st.markdown("**Pseudo-1st-order:**")
                if popt1 is not None:
                    qe1, k1 = popt1
                    st.markdown(f"- qₑ = **{qe1:.5f}** mg/g")
                    st.markdown(f"- k₁ = **{k1:.6f}** min⁻¹")
                    st.markdown(f"- R² = **{r2_1:.5f}**")
                    st.markdown(f"- MAE = **{mae1:.5f}** mg/g")
                else:
                    st.warning("Fit did not converge.")
            with pc2:
                st.markdown("**Pseudo-2nd-order:**")
                if popt2 is not None:
                    qe2, k2 = popt2
                    st.markdown(f"- qₑ = **{qe2:.5f}** mg/g")
                    st.markdown(f"- k₂ = **{k2:.6f}** g·mg⁻¹·min⁻¹")
                    st.markdown(f"- R² = **{r2_2:.5f}**")
                    st.markdown(f"- MAE = **{mae2:.5f}** mg/g")
                else:
                    st.warning("Fit did not converge.")

        # collect for download
        res = {"label": label, "N": N_val, "AdsCode": ads, "AdsConc": adsc,
               "parity": parity, "t": t_arr, "q": q_arr}
        if popt1 is not None:
            qe1, k1 = popt1
            res["1st"] = {"qe": qe1, "k1": k1, "r2": r2_1, "mae": mae1}
        if popt2 is not None:
            qe2, k2 = popt2
            res["2nd"] = {"qe": qe2, "k2": k2, "r2": r2_2, "mae": mae2}
        all_series_results.append(res)

    # ── Combined plots ─────────────────────────────────────────────────────────
    st.divider()
    st.subheader("q vs. time (all series)")
    fig_qt.update_layout(
        xaxis_title="t  (min)", yaxis_title="q  (mg/g)",
        legend=dict(), template="plotly_white", height=520,
    )
    st.plotly_chart(fig_qt, width='stretch')

    st.subheader("Linearised Forms")
    lc1, lc2 = st.columns(2)
    with lc1:
        fig_lin1.update_layout(title="1st order: ln(qₑ − q) vs. t",
                               xaxis_title="t (min)", yaxis_title="ln(qₑ − q(t))",
                               template="plotly_white", height=380)
        st.plotly_chart(fig_lin1, width='stretch')
    with lc2:
        fig_lin2.update_layout(title="2nd order: t/q vs. t",
                               xaxis_title="t (min)", yaxis_title="t / q  (min·g·mg⁻¹)",
                               template="plotly_white", height=380)
        st.plotly_chart(fig_lin2, width='stretch')


    # ══════════════════════════════════════════════════════════════════════════
    # Extra kinetic models — user-selectable
    # ══════════════════════════════════════════════════════════════════════════

    # ── Local model & helper functions ────────────────────────────────────────

    def _nth_order(t, qe, kn, n):
        """Pseudo-nth order integrated: qt = qe − (qe^(1−n) − kn·(1−n)·t)^(1/(1−n))."""
        t = np.asarray(t, dtype=float)
        n = float(n)
        if abs(n - 1.0) < 1e-4:
            return qe * (1.0 - np.exp(-kn * t))
        if abs(n - 2.0) < 1e-4:
            return kn * qe**2 * t / (1.0 + kn * qe * t)
        base = qe ** (1.0 - n) - kn * (1.0 - n) * t
        return qe - np.maximum(base, 1e-300) ** (1.0 / (1.0 - n))

    def _elovich(t, alpha, beta):
        return (1.0 / beta) * np.log(np.maximum(1.0 + alpha * beta * np.asarray(t, float), 1e-300))

    def _avrami(t, qe, kAV, nAV):
        return qe * (1.0 - np.exp(-np.maximum(kAV * np.asarray(t, float), 0.0) ** nAV))

    def _fit_extra_kin(fn, t, q, p0, bounds):
        from scipy.optimize import curve_fit
        from sklearn.metrics import r2_score, mean_absolute_error as _mae
        try:
            popt, _ = curve_fit(fn, t, q, p0=p0, bounds=bounds, maxfev=30_000)
            pred = fn(t, *popt)
            return popt, float(r2_score(q, pred)), float(_mae(q, pred))
        except Exception:
            return None, float("nan"), float("nan")

    def _ipd_linear(t_arr, q_arr):
        """OLS of q vs √t  →  (kid, C, r2, mae)."""
        from sklearn.linear_model import LinearRegression
        from sklearn.metrics import r2_score, mean_absolute_error as _mae
        X = np.sqrt(t_arr).reshape(-1, 1)
        reg = LinearRegression().fit(X, q_arr)
        pred = reg.predict(X)
        return (float(reg.coef_[0]), float(reg.intercept_),
                float(r2_score(q_arr, pred)), float(_mae(q_arr, pred)))

    def _fit_segmented_ipd(t_arr, q_arr):
        """Best two-segment linear fit in √t space. Returns dict or None."""
        from sklearn.linear_model import LinearRegression
        from sklearn.metrics import r2_score, mean_absolute_error as _mae
        n = len(t_arr)
        if n < 5:
            return None
        sqt = np.sqrt(t_arr)
        best_sse, best = np.inf, None
        for sp in range(2, n - 2):
            r1 = LinearRegression().fit(sqt[:sp].reshape(-1, 1), q_arr[:sp])
            r2_ = LinearRegression().fit(sqt[sp:].reshape(-1, 1), q_arr[sp:])
            sse = (np.sum((q_arr[:sp]  - r1.predict(sqt[:sp].reshape(-1, 1))) ** 2) +
                   np.sum((q_arr[sp:]  - r2_.predict(sqt[sp:].reshape(-1, 1))) ** 2))
            if sse < best_sse:
                best_sse = sse
                best = (sp, float(r1.coef_[0]), float(r1.intercept_),
                        float(r2_.coef_[0]), float(r2_.intercept_))
        if best is None:
            return None
        sp, k1, c1, k2, c2 = best
        pred_all = np.concatenate([k1 * sqt[:sp] + c1, k2 * sqt[sp:] + c2])
        return {
            "k1": k1, "C1": c1, "k2": k2, "C2": c2,
            "t1": float(t_arr[sp]),
            "r2":  float(r2_score(q_arr, pred_all)),
            "mae": float(_mae(q_arr, pred_all)),
        }

    # ── Model registry (name → display metadata) ──────────────────────────────
    _EXTRA_KIN_MODELS = {
        "Pseudo-nth-order": {
            "latex":        r"q_t = q_e - \!\left(q_e^{1-n} - k_n(1-n)\,t\right)^{\!\frac{1}{1-n}}",
            "params_label": "qₑ, kₙ, n",
        },
        "Elovich": {
            "latex":        r"q_t = \frac{1}{\beta}\ln\!\left(1+\alpha\beta t\right)",
            "params_label": "α, β",
        },
        "Avrami fraccionary": {
            "latex":        r"q_t = q_e\!\left[1 - e^{-(k_{AV}t)^{n_{AV}}}\right]",
            "params_label": "qₑ, kAV, nAV",
        },
        "IntraParticle Diffusion (IPD)": {
            "latex":        r"q_t = k_{id}\,t^{1/2} + C",
            "params_label": "kid, C",
        },
        "Segmented IPD": {
            "latex":        r"q_t = k_1\,t^{1/2}+C_1\;(\text{phase 1}), \quad q_t = k_2\,t^{1/2}+C_2\;(\text{phase 2})",
            "params_label": "k₁, C₁, k₂, C₂, t₁",
        },
    }

    # ── UI: heading + equation reference ──────────────────────────────────────
    st.divider()
    st.subheader("🔬 Additional kinetics models")
    st.markdown(
        "**PFO** and **PSO** models are the default (parameters obtained in the previous section). "
        "Select here any additional model to fit:"
    )

    with st.expander("📋 Available additional models", expanded=False):
        for _mn, _mi in _EXTRA_KIN_MODELS.items():
            _ca, _cb = st.columns([1, 2])
            _ca.markdown(f"**{_mn}**")
            _ca.caption(f"parameters: {_mi['params_label']}")
            _cb.latex(_mi["latex"])

    selected_extra_kin = st.multiselect(
        "Additional models to fit:",
        options=list(_EXTRA_KIN_MODELS.keys()),
        default=st.session_state.get("_kin_extra_sel", []),
        key="_kin_extra_sel",
    )

    extra_kin_results = []   # persisted to session_state["kin_results"]["extra_models"]

    if selected_extra_kin:
        st.markdown("---")
        for _si, _sr in series_keys.iterrows():
            _N    = _sr["N"]
            _ads  = _sr["AdsCode"]
            _adsc = _sr["AdsConc"]
            _par  = _sr["Parity"] if use_parity else None
            _ph   = _sr.get("pH") if multi_ph else None

            _mask = (
                (avg_df["N"]       == _N)    &
                (avg_df["AdsCode"] == _ads)  &
                (avg_df["AdsConc"] == _adsc)
            )
            if use_parity:
                _mask &= (avg_df["Parity"] == _par)
            if multi_ph:
                _mask &= (avg_df["pH"] == _ph)

            _sub   = avg_df[_mask].sort_values("t (min)")
            _tarr  = _sub["t (min)"].values.astype(float)
            _qarr  = _sub["q (mg/g)"].values.astype(float)

            if use_parity:
                _slbl = (f"pH={_ph} | N={_N} | {_ads}{_adsc:.4g} | {_par} reps"
                         if multi_ph else f"N={_N} | {_ads}{_adsc:.4g} | {_par} reps")
            else:
                _slbl = (f"pH={_ph} | N={_N} | {_ads}{_adsc:.4g}"
                         if multi_ph else f"N={_N} | {_ads}{_adsc:.4g}")

            _color  = series_color(_si)
            _qe_ref = max(_qarr.max() * 1.1, 1e-6)
            _tplot  = np.linspace(max(_tarr.min() * 0.5, 0.1), _tarr.max() * 1.05, 500)

            with st.expander(f"📊  Series: **{_slbl}**", expanded=True):
                _fig_qt  = go.Figure()   # q vs t  (non-IPD models)
                _fig_ipd = None          # q vs √t (IPD models)

                _fig_qt.add_trace(go.Scatter(
                    x=_tarr, y=_qarr, mode="markers", name="Data",
                    marker=dict(color=_color, size=9, symbol="circle"),
                ))

                _has_non_ipd = any(
                    m not in ("IntraParticle Diffusion (IPD)", "Segmented IPD")
                    for m in selected_extra_kin
                )
                _has_ipd = any(
                    m in ("IntraParticle Diffusion (IPD)", "Segmented IPD")
                    for m in selected_extra_kin
                )

                for _mn in selected_extra_kin:
                    _rec = {
                        "model":   _mn,
                        "label":   _slbl,
                        "N":       _N,
                        "AdsCode": _ads,
                        "AdsConc": _adsc,
                    }

                    # ── Pseudo-nth order ───────────────────────────────────────
                    if _mn == "Pseudo-nth-order":
                        _popt, _r2, _mae = _fit_extra_kin(
                            _nth_order, _tarr, _qarr,
                            p0=[_qe_ref, 0.01, 1.5],
                            bounds=([0.0, 1e-9, 0.01], [np.inf, np.inf, 10.0]),
                        )
                        if _popt is not None:
                            _qe, _kn, _n = _popt
                            _rec["params"] = {"qe": _qe, "kn": _kn, "n": _n,
                                              "h":  _kn * _qe ** _n}
                            _rec["r2"], _rec["mae"] = _r2, _mae
                            _fig_qt.add_trace(go.Scatter(
                                x=_tplot, y=_nth_order(_tplot, _qe, _kn, _n),
                                mode="lines", name=_mn, line=dict(width=2),
                            ))
                            st.markdown(
                                f"**{_mn}** — "
                                f"qₑ=`{_qe:.5f}` mg/g | kₙ=`{_kn:.6f}` | n=`{_n:.4f}` | "
                                f"h=kₙ·qₑⁿ=`{_kn*_qe**_n:.5f}` | "
                                f"R²=`{_r2:.5f}` | MAE=`{_mae:.5f}`"
                            )
                        else:
                            st.warning(f"**{_mn}**: fit did not converged.")
                            _rec["params"] = {}
                            _rec["r2"] = _rec["mae"] = float("nan")

                    # ── Elovich ────────────────────────────────────────────────
                    elif _mn == "Elovich":
                        _popt, _r2, _mae = _fit_extra_kin(
                            _elovich, _tarr, _qarr,
                            p0=[max(_qarr[0], 1e-3) if len(_qarr) else 1.0, 1.0],
                            bounds=([1e-9, 1e-9], [np.inf, np.inf]),
                        )
                        if _popt is not None:
                            _a, _b = _popt
                            _rec["params"] = {"alpha": _a, "beta": _b}
                            _rec["r2"], _rec["mae"] = _r2, _mae
                            _fig_qt.add_trace(go.Scatter(
                                x=_tplot, y=_elovich(_tplot, _a, _b),
                                mode="lines", name=_mn, line=dict(width=2, dash="dot"),
                            ))
                            st.markdown(
                                f"**{_mn}** — "
                                f"α=`{_a:.6f}` mg·g⁻¹·min⁻¹ | β=`{_b:.6f}` g·mg⁻¹ | "
                                f"R²=`{_r2:.5f}` | MAE=`{_mae:.5f}`"
                            )
                        else:
                            st.warning(f"**{_mn}**: fit did not converged.")
                            _rec["params"] = {}
                            _rec["r2"] = _rec["mae"] = float("nan")

                    # ── Avrami fraccionario ────────────────────────────────────
                    elif _mn == "Avrami fraccionary":
                        _popt, _r2, _mae = _fit_extra_kin(
                            _avrami, _tarr, _qarr,
                            p0=[_qe_ref, 0.1, 1.0],
                            bounds=([0.0, 1e-9, 0.01], [np.inf, np.inf, 10.0]),
                        )
                        if _popt is not None:
                            _qe, _kAV, _nAV = _popt
                            _rec["params"] = {"qe": _qe, "kAV": _kAV, "nAV": _nAV}
                            _rec["r2"], _rec["mae"] = _r2, _mae
                            _fig_qt.add_trace(go.Scatter(
                                x=_tplot, y=_avrami(_tplot, _qe, _kAV, _nAV),
                                mode="lines", name=_mn, line=dict(width=2, dash="dashdot"),
                            ))
                            st.markdown(
                                f"**{_mn}** — "
                                f"qₑ=`{_qe:.5f}` mg/g | kAV=`{_kAV:.6f}` | nAV=`{_nAV:.4f}` | "
                                f"R²=`{_r2:.5f}` | MAE=`{_mae:.5f}`"
                            )
                        else:
                            st.warning(f"**{_mn}**: fit did not converged.")
                            _rec["params"] = {}
                            _rec["r2"] = _rec["mae"] = float("nan")

                    # ── IPD ────────────────────────────────────────────────────
                    elif _mn == "IntraParticle Diffusion (IPD)":
                        _kid, _C, _r2, _mae = _ipd_linear(_tarr, _qarr)
                        _rec["params"] = {"kid": _kid, "C": _C}
                        _rec["r2"], _rec["mae"] = _r2, _mae
                        if _fig_ipd is None:
                            _fig_ipd = go.Figure()
                            _fig_ipd.add_trace(go.Scatter(
                                x=np.sqrt(_tarr), y=_qarr, mode="markers",
                                name="Data", marker=dict(color=_color, size=9),
                            ))
                        _sq_plot = np.sqrt(_tplot)
                        _fig_ipd.add_trace(go.Scatter(
                            x=_sq_plot, y=_kid * _sq_plot + _C,
                            mode="lines", name=_mn, line=dict(width=2),
                        ))
                        st.markdown(
                            f"**{_mn}** — "
                            f"kid=`{_kid:.6f}` mg·g⁻¹·min⁻½ | C=`{_C:.5f}` mg/g | "
                            f"R²=`{_r2:.5f}` | MAE=`{_mae:.5f}`"
                        )

                    # ── IPD Segmentada ─────────────────────────────────────────
                    elif _mn == "Segmented IPD":
                        _seg = _fit_segmented_ipd(_tarr, _qarr)
                        if _seg:
                            _rec["params"] = {k: v for k, v in _seg.items()
                                              if k not in ("r2", "mae")}
                            _rec["r2"], _rec["mae"] = _seg["r2"], _seg["mae"]
                            _sp = np.searchsorted(_tarr, _seg["t1"])
                            _sq1, _sq2 = np.sqrt(_tarr[:_sp]), np.sqrt(_tarr[_sp:])
                            if _fig_ipd is None:
                                _fig_ipd = go.Figure()
                                _fig_ipd.add_trace(go.Scatter(
                                    x=np.sqrt(_tarr), y=_qarr, mode="markers",
                                    name="Data", marker=dict(color=_color, size=9),
                                ))
                            _fig_ipd.add_trace(go.Scatter(
                                x=_sq1, y=_seg["k1"]*_sq1+_seg["C1"],
                                mode="lines", name="IPD Seg. — phase 1",
                                line=dict(width=2),
                            ))
                            _fig_ipd.add_trace(go.Scatter(
                                x=_sq2, y=_seg["k2"]*_sq2+_seg["C2"],
                                mode="lines", name="IPD Seg. — phase 2",
                                line=dict(width=2, dash="dash"),
                            ))
                            st.markdown(
                                f"**{_mn}** — "
                                f"k₁=`{_seg['k1']:.5f}` | C₁=`{_seg['C1']:.4f}` | "
                                f"k₂=`{_seg['k2']:.5f}` | C₂=`{_seg['C2']:.4f}` | "
                                f"t₁=`{_seg['t1']:.2f}` min | "
                                f"R²=`{_seg['r2']:.5f}` | MAE=`{_seg['mae']:.5f}`"
                            )
                        else:
                            st.warning(f"**{_mn}**: more than ≥5 points per series are needed.")
                            _rec["params"] = {}
                            _rec["r2"] = _rec["mae"] = float("nan")

                    extra_kin_results.append(_rec)

                # ── Plots ──────────────────────────────────────────────────────
                if _has_non_ipd:
                    _fig_qt.update_layout(
                        xaxis_title="t  (min)", yaxis_title="q  (mg/g)",
                        template="plotly_white", height=380,
                        legend=dict(),
                    )
                    st.plotly_chart(_fig_qt, width='stretch')

                if _fig_ipd is not None:
                    _fig_ipd.update_layout(
                        title="IntraParticle Diffusion: q vs √t",
                        xaxis_title="√t  (min^½)", yaxis_title="q  (mg/g)",
                        template="plotly_white", height=380,
                    )
                    st.plotly_chart(_fig_ipd, width='stretch')





    # ── Persist for Download tab ───────────────────────────────────────────────

    # CSV para "q vs time" — ya existe en avg_df, solo renombramos para claridad
    # CSV para formas linealizadas
    lin1_rows, lin2_rows = [], []
    for s in all_series_results:
        t_arr = s["t"]
        q_arr = s["q"]
        lbl   = s["label"]
        qe1_val = s.get("1st", {}).get("qe")
        for ti, qi in zip(t_arr, q_arr):
            lin2_rows.append({"Series": lbl, "t_min": ti,
                              "t_over_q": ti / qi if qi > 1e-12 else None})
            if qe1_val is not None:
                safe = max(qe1_val - qi, 1e-12)
                lin1_rows.append({"Series": lbl, "t_min": ti, "ln_qe_minus_q": np.log(safe)})

    st.session_state["kin_results"] = {
        "use_parity": use_parity,
        "co_correction_applied": apply_corr,
        "files":     [e["name"] for e in kin_entries],
        "df":        raw_df,
        "avg_df":    avg_df,
        "series":    all_series_results,
        "csv_qt":    avg_df[["N", "AdsCode", "AdsConc", "t (min)", "Cₜ (mg/L)", "Cₜ std (mg/L)", "%R", "%R std", "q (mg/g)", "q std"]].copy(),
        "csv_lin1":  pd.DataFrame(lin1_rows),
        "csv_lin2":  pd.DataFrame(lin2_rows),
        "extra_models": extra_kin_results,
    }
