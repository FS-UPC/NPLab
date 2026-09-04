"""
tabs/download.py — Download Results tab.

Assembles a human-readable TXT report from session_state and offers it
as a Streamlit download button.  Also provides per-dataset CSV downloads,
including extra kinetic / isotherm models when selected by the user.
"""

from datetime import datetime
import pandas as pd
import streamlit as st


# ── Formatting helpers ─────────────────────────────────────────────────────────

def _section(title: str) -> str:
    bar = "=" * 60
    return f"\n{bar}\n  {title}\n{bar}\n"


def _subsection(title: str) -> str:
    return f"\n--- {title} ---\n"


def _format_pretty_csv(df: pd.DataFrame) -> bytes:
    """
    Returns a fixed-width, human-readable text table.
    - Floats are rounded to 4 decimal places.
    - Every column is padded so header and values are left-aligned.
    """
    # Round all numeric columns to 4 decimal places
    df = df.copy()
    for col in df.select_dtypes(include="number").columns:
        df[col] = df[col].round(4)

    # Convert everything to strings so we can measure widths
    str_df = df.astype(str)

    # Column width = max of (header length, longest value in column)
    col_widths = {
        col: max(len(col), int(str_df[col].str.len().fillna(0).max()))
        for col in str_df.columns
    }
    def pad_row(values):
        return "  ".join(str(v).ljust(col_widths[col]) for col, v in zip(str_df.columns, values))

    header = pad_row(str_df.columns)
    separator = "  ".join("-" * col_widths[col] for col in str_df.columns)
    rows = [pad_row(row) for _, row in str_df.iterrows()]

    table = "\n".join([header, separator] + rows)
    return table.encode("utf-8")


def _csv_button(label: str, df, filename: str) -> None:
    """Renders two side-by-side download buttons for the same dataset:
    - Left  : pretty fixed-width text table  (*_pretty.csv)
    - Right : standard comma-separated CSV   (*.csv)
    """
    if df is None or (hasattr(df, "empty") and df.empty):
        st.caption(f"_{label}: no data yet_")
        return

    # Derive pretty filename by inserting "_pretty" before the extension
    base, _ext = filename.rsplit(".", 1) if "." in filename else (filename, "csv")
    pretty_filename = f"{base}_pretty.{_ext}"

    col_pretty, col_std = st.columns(2)
    with col_pretty:
        st.download_button(
            label=f"⬇️  {label} (formatted)",
            data=_format_pretty_csv(df),
            file_name=pretty_filename,
            mime="text/plain",
        )
    with col_std:
        st.download_button(
            label=f"⬇️  {label} (CSV)",
            data=df.to_csv(index=False).encode("utf-8"),
            file_name=filename,
            mime="text/csv",
        )

def _extra_models_to_df(extra_list: list) -> pd.DataFrame | None:
    """Flatten the extra_models list into a tidy DataFrame for CSV export."""
    if not extra_list:
        return None
    rows = []
    for rec in extra_list:
        row = {"Model": rec.get("model", ""), "Series": rec.get("label", "")}
        # Series-specific ID columns (present in kinetics but not isotherms, and vice-versa)
        for col in ("N", "AdsCode", "AdsConc", "Ads", "Cads"):
            if col in rec:
                row[col] = rec[col]
        row.update(rec.get("params", {}))
        row["R2"]  = rec.get("r2",  float("nan"))
        row["MAE"] = rec.get("mae", float("nan"))
        rows.append(row)
    return pd.DataFrame(rows) if rows else None


# ── Report builder ─────────────────────────────────────────────────────────────

def build_report() -> str:
    """Return the full report as a plain-text string."""
    lines = []
    ts = datetime.now().strftime("%Y-%m-%d  %H:%M:%S")

    lines.append("NANOPLASTIC ADSORPTION APP — RESULTS REPORT")
    lines.append(f"Generated: {ts}")

    # ── Calibration ────────────────────────────────────────────────────────────
    lines.append(_section("CALIBRATION"))

    lines.append(_subsection("Calibration settings"))
    if "calibration_model" in st.session_state:
        lines.append(f"  Calibration wavelength   : {st.session_state.get('calibration_wavelength', 508.0):.0f} nm")
        lines.append(f"  Polynomial degree        : {st.session_state['calibration_poly_deg']}")
        lines.append(f"  Equation (counts = f(conc)): {st.session_state.get('calibration_equation', 'N/A')}")
        lines.append(f"  R²                       : {st.session_state['calibration_r2']:.6f}")
        lines.append(f"  MAE                      : {st.session_state['calibration_mae']:.6f}")
        _peak_wl = st.session_state.get("calibration_wavelength", 508.0)  
        cr = st.session_state.get("calibration_conc_range")
        if cr:
            lines.append(f"  Concentration range used : {cr[0]:.4g} – {cr[1]:.4g}")
    else:
        lines.append("  ⚠  Calibration has not been set (button 'Set calibration' not pressed).")

    # ── Blank / LOD / LOQ ─────────────────────────────────────────────────────
    lines.append(_section("BLANK FILES — LOD / LOQ"))

    lod_data = st.session_state.get("lod_loq", {})
    if lod_data:
        lines.append(_subsection("Blank signal statistics"))
        lines.append(f"  Number of blanks (n)               : {lod_data['n']}")
        lines.append(f"  Blank mean ȳ_B  (cps)              : {lod_data['y_mean']:.4f}")
        lines.append(f"  Sample std dev sB  (cps)           : {lod_data['y_std']:.6f}")
        lines.append(f"  Calibration slope b  (cps·L/mg)    : {lod_data['b']:.6f}")
        lines.append("")
        lines.append(_subsection("Results"))
        lines.append(f"  LOD = 3  × sB / b = {lod_data['LOD']:.6f}  mg/L")
        lines.append(f"  LOQ = 10 × sB / b = {lod_data['LOQ']:.6f}  mg/L")
    else:
        lines.append(
            "  (LOD/LOQ not yet calculated — load ≥ 10 BLK files in the "
            "Blank tab with a linear calibration active.)"
        )
 

    # ── Drift ──────────────────────────────────────────────────────────────────
    lines.append(_section("DRIFT MONITORING (DRF)"))

    drf = st.session_state.get("drift_results", {})
    if not drf:
        lines.append(
            "  (Drift tab not yet processed — load DRF files and run the analysis.)"
        )
    else:
        wl_win = drf.get("wl_window", (495, 525))

        # ── Type-1 ────────────────────────────────────────────────────────
        lines.append(_subsection("Type 1 — DRF without nanoplastics (N = 0)"))
        t1_avg = drf.get("type1_avg")
        if t1_avg is not None and not t1_avg.empty:
            lines.append(
                f"  Wavelength window for median: "
                f"{wl_win[0]:.0f} – {wl_win[1]:.0f} nm"
            )
            lines.append("")
            lines.append(
                f"  {'t':>4}  {'n_rep':>5}  {'Mean_of_medians (cps)':>22}"
                f"  {'Std (cps)':>10}"
                + (f"  {'Conc (mg/L)':>12}" if "Conc (mg/L)" in t1_avg.columns else "")
            )
            lines.append("  " + "-" * 60)
            for _, row in t1_avg.iterrows():
                conc_str = (
                    f"  {row['Conc (mg/L)']:>12.5f}"
                    if "Conc (mg/L)" in t1_avg.columns
                    else ""
                )
                lines.append(
                    f"  {int(row['t']):>4}  {int(row['n_replicas']):>5}"
                    f"  {row['Mean_of_medians']:>22.4f}"
                    f"  {row['Std_of_medians']:>10.4f}"
                    + conc_str
                )

            ts1 = drf.get("type1_trend_sig", {})
            tc1 = drf.get("type1_trend_conc", {})
            lines.append("")
            lines.append("  Signal drift (cps):")
            lines.append(
                f"    slope = {ts1.get('slope', float('nan')):.6f} cps/step  |  "
                f"intercept = {ts1.get('intercept', float('nan')):.4f}  |  "
                f"R² = {ts1.get('r2', float('nan')):.4f}"
            )
            lines.append("  Concentration drift (mg/L):")
            lines.append(
                f"    slope = {tc1.get('slope', float('nan')):.6f} mg/L/step  |  "
                f"intercept = {tc1.get('intercept', float('nan')):.4f}  |  "
                f"R² = {tc1.get('r2', float('nan')):.4f}"
            )
        else:
            lines.append("  (no Type-1 DRF data)")

        # ── Type-2 ────────────────────────────────────────────────────────
        lines.append(_subsection("Type 2 — DRF with nanoplastics (N > 0)"))
        t2_avg = drf.get("type2_avg")
        if t2_avg is not None and not t2_avg.empty:
            for N_val in sorted(t2_avg["N (mg/L)"].unique()):
                sub = t2_avg[t2_avg["N (mg/L)"] == N_val].sort_values("t")
                lines.append(f"\n  N = {N_val} mg/L (nominal):")
                lines.append(
                    f"    {'t':>4}  {'n_rep':>5}  {f'Avg counts @ {_peak_wl:.0f} nm':>20}"
                    f"  {'Std':>10}"
                )
                lines.append("    " + "-" * 44)
                for _, row in sub.iterrows():
                    lines.append(
                        f"    {int(row['t']):>4}  {int(row['n_replicas']):>5}"
                        f"  {row['Avg_counts']:>20.4f}  {row['Std_counts']:>10.4f}"
                    )

                ts2 = drf.get("type2_trend_sig", {}).get(N_val, {})
                tc2 = drf.get("type2_trend_conc", {}).get(N_val, {})
                lines.append(f"    Signal trend:  slope = {ts2.get('slope', float('nan')):.6f} cps/step"
                             f"  |  R² = {ts2.get('r2', float('nan')):.4f}")
                lines.append(f"    Conc trend:    slope = {tc2.get('slope', float('nan')):.6f} mg/L/step"
                             f"  |  R² = {tc2.get('r2', float('nan')):.4f}")

            dev_df = drf.get("type2_deviation")
            if dev_df is not None and not dev_df.empty:
                lines.append("")
                lines.append("  Deviation of recovered concentration from nominal:")
                lines.append(
                    f"  {'N_nom':>8}  {'t':>4}  {'Counts':>12}"
                    f"  {'Calc (mg/L)':>12}  {'Dev (mg/L)':>12}  {'Dev %':>8}"
                )
                lines.append("  " + "-" * 62)
                for _, row in dev_df.iterrows():
                    lines.append(
                        f"  {row['N_nominal (mg/L)']:>8.2f}  {int(row['t']):>4}"
                        f"  {row['Avg_counts (cps)']:>12.4f}"
                        f"  {row['Conc_calc (mg/L)']:>12.5f}"
                        f"  {row['Dev_from_nom (mg/L)']:>12.5f}"
                        f"  {row['Dev_%']:>8.3f}"
                    )
        else:
            lines.append("  (no Type-2 DRF data)")

    # ── Control ────────────────────────────────────────────────────────────────
    lines.append(_section("CONTROL EXPERIMENTS (CTRL)"))

    ctrl = st.session_state.get("ctrl_results", {})
    if not ctrl:
        lines.append(
            "  (Control tab not yet processed — load CTRL files and run the analysis.)"
        )
    else:
        wl_win_c = ctrl.get("wl_window", (495, 525))

        # ── Type-1 ────────────────────────────────────────────────────────
        lines.append(_subsection("Type 1 — CTRL without nanoplastics (N = 0)"))
        t1_avg = ctrl.get("type1_avg")
        if t1_avg is not None and not t1_avg.empty:
            lines.append(
                f"  Wavelength window for median: "
                f"{wl_win_c[0]:.0f} – {wl_win_c[1]:.0f} nm"
            )
            lines.append(
                f"  {'Ads_series':<24}  {'C-code':>7}  "
                f"{'Mean_median (cps)':>18}  {'Std (cps)':>10}  {'n_rep':>5}"
            )
            lines.append("  " + "-" * 70)
            for _, row in t1_avg.iterrows():
                lines.append(
                    f"  {str(row['Ads_series']):<24}  "
                    f"{str(row['Ccode']):>7}  "
                    f"{float(row['Mean_median']):>18.4f}  "
                    f"{float(row['Std_median']):>10.4f}  "
                    f"{int(row['n_replicas']):>5}"
                )
        else:
            lines.append("  (no Type-1 CTRL data)")

        # ── Type-2 ────────────────────────────────────────────────────────
        lines.append(_subsection("Type 2 — CTRL with nanoplastics (N > 0)"))
        t2_avg = ctrl.get("type2_avg")
        if t2_avg is not None and not t2_avg.empty:
            has_conc = "Conc_calc (mg/L)" in t2_avg.columns
            hdr = (
                f"  {'N_series':<18}  {'Ads_series':<20}  {'C-code':>7}"
                f"  {'Avg_counts':>12}  {'Std':>10}"
            )
            if has_conc:
                hdr += f"  {'Conc (mg/L)':>12}"
            lines.append(hdr)
            lines.append("  " + "-" * (80 + (13 if has_conc else 0)))
            for _, row in t2_avg.iterrows():
                line = (
                    f"  {str(row['N_series']):<18}  "
                    f"{str(row['Ads_series']):<20}  "
                    f"{str(row['Ccode']):>7}  "
                    f"{float(row['Avg_counts']):>12.4f}  "
                    f"{float(row['Std_counts']):>10.4f}"
                )
                if has_conc:
                    line += f"  {float(row['Conc_calc (mg/L)']):>12.5f}"
                lines.append(line)

            dev_df = ctrl.get("type2_deviation")
            if dev_df is not None and not dev_df.empty:
                lines.append("")
                lines.append("  Deviation from nominal concentration:")
                lines.append(
                    f"  {'N_nom':>8}  {'Ads_series':<20}  {'C-code':>7}"
                    f"  {'Calc (mg/L)':>12}  {'Dev (mg/L)':>12}  {'Dev %':>8}"
                )
                lines.append("  " + "-" * 76)
                for _, row in dev_df.iterrows():
                    lines.append(
                        f"  {float(row['N_nominal (mg/L)']):>8.2f}  "
                        f"{str(row['Ads_series']):<20}  "
                        f"{str(row['C-code']):>7}  "
                        f"{float(row['Conc_calc (mg/L)']):>12.5f}  "
                        f"{float(row['Dev_from_nom (mg/L)']):>12.5f}  "
                        f"{float(row['Dev_%']):>8.3f}"
                    )
        else:
            lines.append("  (no Type-2 CTRL data)")

        # ── Co correction ────────────────────────────────────────────────
        co_df = ctrl.get("co_correction")
        if co_df is not None and not co_df.empty:
            lines.append(_subsection("Co correction — special control (C000)"))
            lines.append(
                f"  {'pH':>6}  {'N_nominal':>10}  {'N_measured':>11}"
                f"  {'Ratio':>8}  {'Dev_%':>8}"
            )
            lines.append("  " + "-" * 50)
            for _, row in co_df.iterrows():
                ph_str = f"{row['pH']:>6.4g}" if pd.notna(row["pH"]) else f"{'–':>6}"
                lines.append(
                    f"  {ph_str}  {row['N_nominal']:>10.2f}  {row['N_measured']:>11.4f}"
                    f"  {row['Ratio (measured/nominal)']:>8.4f}  {row['Deviation (%)']:>8.3f}"
                )

        # ── Cross-type summary ─────────────────────────────────────────────
        cross = ctrl.get("cross_summary")
        if cross is not None and not cross.empty:
            lines.append(_subsection("Cross-type condition summary"))
            lines.append(
                f"  {'Type':<8}  {'Ads/N_series':<24}  {'C-code':>7}"
                f"  {'Value':>12}  {'Std':>10}  {'Unit'}"
            )
            lines.append("  " + "-" * 80)
            for _, row in cross.iterrows():
                lines.append(
                    f"  {str(row['Type']):<8}  "
                    f"{str(row['Ads_series']):<24}  "
                    f"{str(row['Ccode']):>7}  "
                    f"{float(row['Value']):>12.4f}  "
                    f"{float(row['Std']):>10.4f}  "
                    f"{str(row['Unit'])}"
                )

    # ── Kinetics ───────────────────────────────────────────────────────────────
    lines.append(_section("KINETICS"))

    kin = st.session_state.get("kin_results", {})

    lines.append(_subsection("Parameters"))
    if kin:
        lines.append(f"  Formula : q(t) = (N − Cₜ) / AdsConc  [mg/g]")
        lines.append(f"  N and AdsConc are read from each filename.")
        co_on = kin.get("co_correction_applied", False)
        lines.append(
            f"  Co correction (from Control C000) applied : {'YES' if co_on else 'no'}"
        )
    else:
        lines.append("  (kinetics tab not yet processed)")

    lines.append(_subsection("Processed data (t, Ct, q(t))"))
    df = kin.get("df")
    if df is not None:
        lines.append(f"  {'Filename':<35} {'t (min)':>8} {'Ct (mg/L)':>12} {'q(t) (mg/g)':>14}")
        lines.append("  " + "-" * 73)
        for _, row in df.iterrows():
            lines.append(
                f"  {row['Filename']:<35} {row['t (min)']:>8.1f} "
                f"{row['Cₜ (mg/L)']:>12.4f} {row['q(t) (mg/g)']:>14.6f}"
            )
    else:
        lines.append("  (no data)")

    lines.append(_subsection("Averaged data (mean Cₜ and q per series)"))
    avg_df_rep = kin.get("avg_df")
    if avg_df_rep is not None and not avg_df_rep.empty:
        lines.append(
            f"  {'N':>6} {'AdsCode':>10} {'AdsConc':>8} {'t (min)':>8} {'Reps':>5}"
            f" {'Cₜ (mg/L)':>12} {'Cₜ std':>10} {'q (mg/g)':>12} {'q std':>10}"
        )
        lines.append("  " + "-" * 87)
        for _, row in avg_df_rep.iterrows():
            lines.append(
                f"  {row['N']:>6.1f} {str(row['AdsCode']):>10} {row['AdsConc']:>8.2f}"
                f" {row['t (min)']:>8.1f} {int(row['Reps']):>5}"
                f" {row['Cₜ (mg/L)']:>12.4f} {row.get('Cₜ std (mg/L)', float('nan')):>10.4f}"
                f" {row['q (mg/g)']:>12.6f} {row.get('q std', float('nan')):>10.6f}"
            )
    else:
        lines.append("  (no data)")

    lines.append(_subsection("Default kinetic models  —  Pseudo-1st & 2nd order  (per series)"))
    series_list = kin.get("series", [])
    if series_list:
        for s in series_list:
            lines.append(f"\n  Series: {s['label']}")
            lines.append(f"    N={s['N']}  AdsCode={s['AdsCode']}  AdsConc={s['AdsConc']} g/L")
            m1 = s.get("1st")
            if m1:
                lines.append(f"    1st-order: qₑ={m1['qe']:.6f} mg/g  k₁={m1['k1']:.8f} min⁻¹  "
                             f"R²={m1['r2']:.6f}  MAE={m1['mae']:.6f} mg/g")
            else:
                lines.append("    1st-order: fit did not converge")
            m2 = s.get("2nd")
            if m2:
                lines.append(f"    2nd-order: qₑ={m2['qe']:.6f} mg/g  k₂={m2['k2']:.8f} g·mg⁻¹·min⁻¹  "
                             f"R²={m2['r2']:.6f}  MAE={m2['mae']:.6f} mg/g")
            else:
                lines.append("    2nd-order: fit did not converge")
    elif kin:
        lines.append("  (kinetics tab not yet processed)")
    else:
        lines.append("  (not available)")

    # ── Extra kinetic models ───────────────────────────────────────────────────
    lines.append(_subsection("Additional kinetic models  (user-selected)"))
    extra_kin = kin.get("extra_models", [])
    if extra_kin:
        # Group by model name
        seen_models: dict[str, list] = {}
        for rec in extra_kin:
            seen_models.setdefault(rec["model"], []).append(rec)
        for model_name, recs in seen_models.items():
            lines.append(f"\n  ── Model: {model_name} ──")
            for rec in recs:
                lines.append(f"    Series : {rec['label']}")
                params = rec.get("params", {})
                if params:
                    for p_name, p_val in params.items():
                        try:
                            lines.append(f"      {p_name:<10} = {float(p_val):.6f}")
                        except (TypeError, ValueError):
                            lines.append(f"      {p_name:<10} = {p_val}")
                r2_  = rec.get("r2",  float("nan"))
                mae_ = rec.get("mae", float("nan"))
                lines.append(f"      R²         = {r2_:.6f}")
                lines.append(f"      MAE        = {mae_:.6f}")
    else:
        lines.append("  (no additional kinetic models selected)")

    # ── Isotherms ──────────────────────────────────────────────────────────────
    lines.append(_section("ISOTHERMS"))

    iso = st.session_state.get("iso_results", {})

    lines.append(_subsection("Fixed parameters"))
    if iso:
        lines.append(f"  Calibration wavelength = {_peak_wl:.0f} nm")
        lines.append(f"  qₑ formula    = (N − Cₑ) / C_ads   [mg/g]")
        lines.append(f"  R_L formula   = 1 / (1 + K_L · N_max)  — N_max per series")
        co_on_iso = iso.get("co_correction_applied", False)
        lines.append(
            f"  Co correction (from Control C000) applied : {'YES' if co_on_iso else 'no'}"
        )
    else:
        lines.append("  (not available)")

    lines.append(_subsection("Raw data summary"))
    raw_df = iso.get("raw_df")
    if raw_df is not None and not raw_df.empty:
        header = f"  {'Filename':<35} {'N':>8} {'Ads':>8} {'Cads':>5} {'R':>3} {'Counts':>12} {'Ce (mg/L)':>10} {'qe (mg/g)':>12}"
        lines.append(header)
        lines.append("  " + "-" * 92)
        for _, row in raw_df.iterrows():
            lines.append(
                f"  {str(row['Filename']):<35} {row['N (mg/L)']:>8.2f} "
                f"{str(row['Ads']):>8} {int(row['Cads (g/L)']):>5} {int(row['R']):>3} "
                f"{row[f'Counts @ {_peak_wl:.0f} nm']:>12.2f} {row['Cₑ (mg/L)']:>10.4f} "
                f"{row['qₑ (mg/g)']:>12.6f}"
            )
    else:
        lines.append("  (no data)")

    lines.append(_subsection("Averaged data (mean Cₑ and qₑ per N × Ads × Cads)"))
    avg_df_iso = iso.get("avg_df")
    if avg_df_iso is not None and not avg_df_iso.empty:
        has_ph = "pH" in avg_df_iso.columns
        hdr = (
            f"  {'pH':>6} {'N (mg/L)':>10} {'Ads':>8} {'Cads (g/L)':>10} {'Reps':>5}"
            f" {'Cₑ (mg/L)':>12} {'Ce std':>10} {'qₑ (mg/g)':>12} {'qe std':>10}"
            if has_ph else
            f"  {'N (mg/L)':>10} {'Ads':>8} {'Cads (g/L)':>10} {'Reps':>5}"
            f" {'Cₑ (mg/L)':>12} {'Ce std':>10} {'qₑ (mg/g)':>12} {'qe std':>10}"
        )
        lines.append(hdr)
        lines.append("  " + "-" * (93 if has_ph else 85))
        for _, row in avg_df_iso.iterrows():
            ph_str = f"  {row['pH']:>6.4g}" if has_ph else ""
            lines.append(
                f"{ph_str}  {row['N (mg/L)']:>10.2f} {str(row['Ads']):>8}"
                f" {row['Cads (g/L)']:>10.2f} {int(row['Reps']):>5}"
                f" {row['Cₑ (mg/L)']:>12.4f} {row.get('Ce_std', float('nan')):>10.4f}"
                f" {row['qₑ (mg/g)']:>12.6f} {row.get('qe_std', float('nan')):>10.6f}"
            )
    else:
        lines.append("  (no data)")

    lines.append(_subsection("Default isotherm models  —  Langmuir & Freundlich  (per series)"))
    if not iso.get("series"):
        lines.append("  (isotherms tab not yet processed)")
    else:
        for s in iso["series"]:
            # Check if there is actually data to report
            if len(s['Ce']) == 0:
                lines.append(f"\n  ── Series: {s['label']} ──")
                lines.append("      No data points available.")
                continue # Skip to the next series

            lines.append(f"\n  ── Series: {s['label']} ──")
            lines.append(f"     Data points : {len(s['Ce'])} (after averaging replicates)")
            lines.append(f"     Ce range    : {s['Ce'].min():.4g} – {s['Ce'].max():.4g} mg/L")
            lines.append(f"     qe range    : {s['qe'].min():.4g} – {s['qe'].max():.4g} mg/g")

            lm = s.get("langmuir")
            lines.append("")
            lines.append("     Langmuir model:  qe = qmax·KL·Ce / (1 + KL·Ce)")
            if lm:
                RL = lm["RL"]
                interp = ("unfavorable" if RL > 1 else
                          "linear"      if abs(RL - 1) < 1e-4 else
                          "favorable"   if RL > 0 else "irreversible")
                lines.append(f"       q_max  = {lm['qmax']:.6f}  mg/g")
                lines.append(f"       K_L    = {lm['KL']:.8f}  L/mg")
                lines.append(f"       C₀ (N_max for R_L) = {s.get('C0_ref', '?'):.4g} mg/L")
                lines.append(f"       R_L    = {lm['RL']:.6f}  → {interp}")
                lines.append(f"       R²     = {lm['r2']:.6f}")
                lines.append(f"       MAE    = {lm['mae']:.6f}  mg/g")
            else:
                lines.append("       Fit did not converge.")

            fm = s.get("freundlich")
            lines.append("")
            lines.append("     Freundlich model:  qe = KF·Ce^(1/n)")
            if fm:
                fav = ("favorable" if fm["inv_n"] < 1
                       else "linear" if abs(fm["inv_n"] - 1) < 1e-3
                       else "unfavorable")
                lines.append(f"       K_F    = {fm['KF']:.6f}  (mg/g)·(L/mg)^(1/n)")
                lines.append(f"       n      = {fm['n']:.6f}")
                lines.append(f"       1/n    = {fm['inv_n']:.6f}  → {fav}")
                lines.append(f"       R²     = {fm['r2']:.6f}")
                lines.append(f"       MAE    = {fm['mae']:.6f}  mg/g")
            else:
                lines.append("       Fit did not converge.")

    # ── Extra isotherm models ──────────────────────────────────────────────────
    lines.append(_subsection("Additional isotherm models  (user-selected)"))
    extra_iso = iso.get("extra_models", [])
    if extra_iso:
        seen_iso: dict[str, list] = {}
        for rec in extra_iso:
            seen_iso.setdefault(rec["model"], []).append(rec)
        for model_name, recs in seen_iso.items():
            lines.append(f"\n  ── Model: {model_name} ──")
            for rec in recs:
                lines.append(f"    Series : {rec['label']}")
                params = rec.get("params", {})
                if params:
                    for p_name, p_val in params.items():
                        try:
                            lines.append(f"      {p_name:<10} = {float(p_val):.6f}")
                        except (TypeError, ValueError):
                            lines.append(f"      {p_name:<10} = {p_val}")
                r2_  = rec.get("r2",  float("nan"))
                mae_ = rec.get("mae", float("nan"))
                lines.append(f"      R²         = {r2_:.6f}")
                lines.append(f"      MAE        = {mae_:.6f}")
    else:
        lines.append("  (no additional isotherm models selected)")

    lines.append("\n" + "=" * 60)
    lines.append("END OF REPORT")
    lines.append("=" * 60 + "\n")

    return "\n".join(lines)


# ── Main render ────────────────────────────────────────────────────────────────

def render() -> None:
    """Main entry-point; call from app.py inside the Download Results tab context."""
    st.subheader("📥 Download Results")
    st.markdown(
        "This tab assembles a plain-text summary of everything computed so far "
        "— calibration settings, kinetic parameters and the processed data table — "
        "and lets you download it as a `.txt` file."
    )

    # ── Status summary ─────────────────────────────────────────────────────────
    has_cal  = "calibration_model" in st.session_state
    has_blk  = "lod_loq"          in st.session_state
    has_drf  = "drift_results"    in st.session_state
    has_ctrl = "ctrl_results"     in st.session_state
    has_kin  = "kin_results"      in st.session_state
    has_iso  = "iso_results"      in st.session_state

    col1, col2, col3, col4, col5, col6 = st.columns(6)
    with col1:
        if has_cal:
            deg = st.session_state["calibration_poly_deg"]
            st.success(f"✅ Calibration (deg {deg})")
        else:
            st.warning("⚠️ Calibration not set")
    with col2:
        if has_blk:
            lod_d = st.session_state["lod_loq"]
            st.success(
                f"✅ Blank — LOD {lod_d['LOD']:.4f} / "
                f"LOQ {lod_d['LOQ']:.4f} mg/L"
            )
        else:
            st.warning("⚠️ Blank not processed")
    with col3:
        if has_drf:
            drf_d = st.session_state["drift_results"]
            n1 = len(drf_d.get("type1_per_file", pd.DataFrame()))
            n2 = len(drf_d.get("type2_per_file", pd.DataFrame()))
            st.success(f"✅ Drift (T1: {n1}, T2: {n2})")
        else:
            st.warning("⚠️ Drift not processed")
    with col4:
        if has_ctrl:
            ctrl_d = st.session_state["ctrl_results"]
            n1c = len(ctrl_d.get("type1_per_file", pd.DataFrame()))
            n2c = len(ctrl_d.get("type2_per_file", pd.DataFrame()))
            st.success(f"✅ Control (T1: {n1c}, T2: {n2c})")
        else:
            st.warning("⚠️ Control not processed")
    with col5:
        if has_kin:
            n_pts = len(st.session_state["kin_results"].get("df", []))
            extra_kin_n = len(st.session_state["kin_results"].get("extra_models", []))
            label_kin = f"✅ Kinetics ({n_pts} pts"
            label_kin += f", {extra_kin_n} extra)" if extra_kin_n else ")"
            st.success(label_kin)
        else:
            st.warning("⚠️ Kinetics not processed")
    with col6:
        if has_iso:
            n_ptsiso = len(st.session_state["iso_results"].get("raw_df", []))
            extra_iso_n = len(st.session_state["iso_results"].get("extra_models", []))
            label_iso = f"✅ Isotherms ({n_ptsiso} pts"
            label_iso += f", {extra_iso_n} extra)" if extra_iso_n else ")"
            st.success(label_iso)
        else:
            st.warning("⚠️ Isotherms not processed")

    st.divider()

    # ── Preview ────────────────────────────────────────────────────────────────
    report_text = build_report()

    with st.expander("👁  Preview report", expanded=False):
        st.code(report_text, language="")

    # ── Download TXT ───────────────────────────────────────────────────────────
    ts_file  = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"Nanoplastic_adsorption_results_{ts_file}.txt"

    st.download_button(
        label="⬇️  Download results as TXT",
        data=report_text.encode("utf-8"),
        file_name=filename,
        mime="text/plain",
        type="primary",
        width='stretch',
    )
    st.caption(f"File will be saved as: `{filename}`")

    # ── CSV Downloads ──────────────────────────────────────────────────────────
    st.divider()
    st.subheader("📊 Download chart data as CSV")

    # ── Calibration ────────────────────────────────────────────────────────────
    st.markdown("**Calibration**")
    _csv_button(
        "Calibration curve data",
        st.session_state.get("cal_curve_df"),
        "calibration_curve.csv",
    )

    # ── Blank / LOD / LOQ ──────────────────────────────────────────────────────
    st.markdown("**Blank (BLK) — LOD / LOQ**")
    lod_data_dl = st.session_state.get("lod_loq", {})
    _csv_button(
        "Blank BLK — counts per file",
        lod_data_dl.get("blk_df"),
        "blank_BLK_counts.csv",
    )
    if lod_data_dl:
        import pandas as _pd
        _lod_summary = _pd.DataFrame([{
            "n_blanks":      lod_data_dl["n"],
            "y_mean_cps":    lod_data_dl["y_mean"],
            "sB_cps":        lod_data_dl["y_std"],
            "b_cps_per_mgL": lod_data_dl["b"],
            "LOD_mg_per_L":  lod_data_dl["LOD"],
            "LOQ_mg_per_L":  lod_data_dl["LOQ"],
        }])
        _csv_button(
            "LOD / LOQ — summary",
            _lod_summary,
            "blank_LOD_LOQ_summary.csv",
        )

    # ── Drift ──────────────────────────────────────────────────────────────────
    st.markdown("**Drift (DRF)**")
    drf_d = st.session_state.get("drift_results", {})
    _csv_button(
        "Type-1 per-file medians",
        drf_d.get("type1_per_file"),
        "drift_type1_per_file.csv",
    )
    _csv_button(
        "Type-1 mean of medians per time-point",
        drf_d.get("type1_avg"),
        "drift_type1_avg_per_timepoint.csv",
    )
    _csv_button(
        "Type-2 per-file counts @ 508 nm",
        drf_d.get("type2_per_file"),
        "drift_type2_per_file.csv",
    )
    _csv_button(
        "Type-2 averaged counts per (N, t) group",
        drf_d.get("type2_avg"),
        "drift_type2_avg_per_group.csv",
    )
    _csv_button(
        "Type-2 deviation from nominal concentration",
        drf_d.get("type2_deviation"),
        "drift_type2_deviation.csv",
    )

    # ── Control ────────────────────────────────────────────────────────────────
    st.markdown("**Control (CTRL)**")
    ctrl_d = st.session_state.get("ctrl_results", {})
    _csv_button(
        "Type-1 per-file medians",
        ctrl_d.get("type1_per_file"),
        "control_type1_per_file.csv",
    )
    _csv_button(
        "Type-1 mean of medians per condition",
        ctrl_d.get("type1_avg"),
        "control_type1_avg_per_condition.csv",
    )
    _csv_button(
        "Type-2 per-file counts @ 508 nm",
        ctrl_d.get("type2_per_file"),
        "control_type2_per_file.csv",
    )
    _csv_button(
        "Type-2 averaged counts & concentration per condition",
        ctrl_d.get("type2_avg"),
        "control_type2_avg_per_condition.csv",
    )
    _csv_button(
        "Type-2 deviation from nominal concentration",
        ctrl_d.get("type2_deviation"),
        "control_type2_deviation.csv",
    )
    _csv_button(
        "Cross-type condition summary",
        ctrl_d.get("cross_summary"),
        "control_cross_type_summary.csv",
    )
    _csv_button(
        "Co correction — measured vs. nominal (C000)",
        ctrl_d.get("co_correction"),
        "control_co_correction.csv",
    )

    # ── Kinetics — default models ──────────────────────────────────────────────
    st.markdown("**Kinetics — default models (PFO & PSO)**")
    kin = st.session_state.get("kin_results", {})
    _csv_button("q vs. time (all series)",  kin.get("csv_qt"),   "kinetics_q_vs_time.csv")
    _csv_button("Linearised 1st order",     kin.get("csv_lin1"), "kinetics_linear_1st_order.csv")
    _csv_button("Linearised 2nd order",     kin.get("csv_lin2"), "kinetics_linear_2nd_order.csv")

    # ── Kinetics — extra models ────────────────────────────────────────────────
    extra_kin = kin.get("extra_models", [])
    if extra_kin:
        st.markdown("**Kinetics — additional models**")
        _csv_button(
            "All additional kinetic model parameters",
            _extra_models_to_df(extra_kin),
            "kinetics_extra_models_params.csv",
        )
        # One CSV per model type
        model_names_kin = sorted({r["model"] for r in extra_kin})
        for m_name in model_names_kin:
            subset = [r for r in extra_kin if r["model"] == m_name]
            safe_name = m_name.lower().replace(" ", "_").replace("(", "").replace(")", "")
            _csv_button(
                f"  {m_name} — parameters per series",
                _extra_models_to_df(subset),
                f"kinetics_{safe_name}_params.csv",
            )

    # ── Isotherms — default models ─────────────────────────────────────────────
    st.markdown("**Isotherms — default models (Langmuir & Freundlich)**")
    iso = st.session_state.get("iso_results", {})
    _csv_button("Averaged data — Cₑ & qₑ per (N × Ads × Cads)", iso.get("csv_iso"), "isotherms_averaged_data.csv")
    _csv_button("Linearised Langmuir (Ce/qe vs Ce)",     iso.get("csv_lin_lang"),  "isotherms_linear_langmuir.csv")
    _csv_button("Linearised Freundlich (ln qe vs ln Ce)",iso.get("csv_lin_freun"), "isotherms_linear_freundlich.csv")

    # ── Isotherms — extra models ───────────────────────────────────────────────
    extra_iso = iso.get("extra_models", [])
    if extra_iso:
        st.markdown("**Isotherms — additional models**")
        _csv_button(
            "All additional isotherm model parameters",
            _extra_models_to_df(extra_iso),
            "isotherms_extra_models_params.csv",
        )
        model_names_iso = sorted({r["model"] for r in extra_iso})
        for m_name in model_names_iso:
            subset = [r for r in extra_iso if r["model"] == m_name]
            safe_name = m_name.lower().replace(" ", "_").replace("(", "").replace(")", "").replace("-", "")
            _csv_button(
                f"  {m_name} — parameters per series",
                _extra_models_to_df(subset),
                f"isotherms_{safe_name}_params.csv",
            )
