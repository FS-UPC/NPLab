"""
tabs/visualization.py — Visualization tab.

Lets the user load one or more spectral files (from folder or upload) and displays:
  1. Counts vs Wavelength
  2. log₁₀(Counts) vs Wavelength
Legend label = filename. Each file gets a distinct colour.

Accepted file formats:
  • TXT  — standard format with text header, data as "wavelength,counts" pairs
  • ASC  — headerless, 4 integer fields per line separated by comma
            (comma acts as both column separator and decimal point)
  • AKN  — same binary layout as ASC, different extension
"""

import os

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from utils import list_files, read_spectral_bytes, read_spectral_file

# Distinct colour palette (cycles if more files than colours)
_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]

# Extensions shown in the uploader widget
_UPLOAD_EXTENSIONS = ["txt", "asc", "akn"]


def render(data_dir: str = "") -> None:
    st.subheader("📂 Input Files — Visualization")

    source = st.radio(
        "File source:",
        ["📁 Folder", "⬆️ Upload"],
        horizontal=True,
        key="viz_source",
    )

    entries = []   # list of {"name": str, "df": DataFrame}

    # ── Folder mode ───────────────────────────────────────────────────────────
    if source == "📁 Folder":
        all_files = list_files(data_dir, r".*")   # all .txt / .asc / .akn
        if not all_files:
            st.warning(
                f"No spectral files (.txt, .asc, .akn) found in **{data_dir}**. "
                "Choose a different subfolder in the sidebar or switch to Upload mode."
            )
            return

        selected = st.multiselect(
            "Choose spectral files to visualize:",
            options=all_files,
            default=[],
            key="viz_folder_selector",
        )
        if not selected:
            st.info("Select at least one file to visualize.")
            return

        for fname in selected:
            df = read_spectral_file(os.path.join(data_dir, fname))
            entries.append({"name": fname, "df": df})

    # ── Upload mode ───────────────────────────────────────────────────────────
    else:
        uploaded = st.file_uploader(
            "Upload spectral files (.txt, .asc, .akn)",
            type=_UPLOAD_EXTENSIONS,
            accept_multiple_files=True,
            key="viz_uploader",
        )
        if not uploaded:
            st.info("Upload at least one spectral file to visualize.")
            return

        for f in uploaded:
            df = read_spectral_bytes(f)
            entries.append({"name": f.name, "df": df})

    # ── Filter empty / unreadable files ───────────────────────────────────────
    bad = [e["name"] for e in entries if e["df"].empty]
    for name in bad:
        st.warning(f"Skipped (no readable data): {name}")

    entries = [e for e in entries if not e["df"].empty]
    if not entries:
        st.error("None of the selected files contained readable spectral data.")
        return

    st.success(f"{len(entries)} file(s) loaded.")
    st.divider()

    # ── Build both figures ────────────────────────────────────────────────────
    fig_lin = go.Figure()
    fig_log = go.Figure()

    for i, e in enumerate(entries):
        color          = _PALETTE[i % len(_PALETTE)]
        df             = e["df"]
        name           = e["name"]
        counts_clipped = df["counts"].clip(lower=1e-10)

        shared = dict(
            x=df["wavelength"],
            mode="lines",
            name=name,
            line=dict(color=color, width=2),
        )

        fig_lin.add_trace(go.Scatter(y=df["counts"],                 **shared))
        fig_log.add_trace(go.Scatter(y=np.log10(counts_clipped),     **shared))

    common_layout = dict(
        xaxis_title="Wavelength (nm)",
        template="plotly_white",
        height=420,
        legend=dict(),
    )

    fig_lin.update_layout(
        yaxis_title="Counts",
        title="Counts vs Wavelength",
        **common_layout,
    )
    fig_log.update_layout(
        yaxis_title="log₁₀(Counts)",
        title="log₁₀(Counts) vs Wavelength",
        **common_layout,
    )

    st.plotly_chart(fig_lin, width='stretch')
    st.plotly_chart(fig_log, width='stretch')
