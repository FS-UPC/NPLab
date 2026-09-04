"""
app.py — main entry point for NPLab (NanoPlastic Adsorption Analysis).

Run with:
    streamlit run app.py

Handles login, the sidebar data source (folder picker or ZIP upload) and the
tab layout; each tab is rendered by its own module.

Module structure:
    app.py                ← this file (login, sidebar, tab layout)
    utils.py              ← shared helpers (spectrum I/O, filename parsers,
                            calibration inversion, kinetic models, colours)
    tabs/visualization.py ← Visualization tab
    tabs/calibration.py   ← Calibration tab
    tabs/blank.py         ← Blank / LOD / LOQ tab
    tabs/drift.py         ← Laser drift monitoring tab
    tabs/control.py       ← Control experiments tab (incl. Co correction)
    tabs/kinetics.py      ← Kinetics tab
    tabs/isotherms.py     ← Isotherms tab
    tabs/download.py      ← Download Results tab
    tabs/help.py          ← Help tab
"""

import os
import re
import zipfile
import tempfile
import shutil

import streamlit as st

# Import tab renderers
from tabs.visualization import render as render_visualization
from tabs.calibration   import render as render_calibration
from tabs.blank         import render as render_blank
from tabs.drift         import render as render_drift
from tabs.control       import render as render_control
from tabs.kinetics      import render as render_kinetics
from tabs.isotherms     import render as render_isotherms
from tabs.download      import render as render_download
from tabs.help          import render as render_help
from utils import folder_picker

# ── Page config ────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="NPLab",
    page_icon="🔬",
    layout="wide",
)

# ── Authentication ─────────────────────────────────────────────────────────────

def check_credentials(username: str, password: str) -> str | None:
    users = st.secrets.get("credentials", {}).get("usernames", {})
    user = users.get(username)
    if user and user["password"] == password:
        return user["name"]
    return None

if st.session_state.get("authenticated") is not True:
    st.title("🔬 NPLab - NanoPlastic Adsorption Analysis")
    st.subheader("Login")
    username = st.text_input("User")
    password = st.text_input("Password", type="password")
    if st.button("Login", type="primary"):
        name = check_credentials(username, password)
        if name:
            st.session_state["authenticated"] = True
            st.session_state["user_name"] = name
            st.session_state["username"] = username
            st.rerun()
        else:
            st.error("❌ Incorrect user or password.")
    st.stop()


# ── From here, user is authenticated ──────────────────────────────────────────
with st.sidebar:
    st.write(f"👤 {st.session_state['user_name']}")
    if st.button("Close session"):
        # Clean up any ZIP temp directory before logging out
        _tmp = st.session_state.get("_zip_tmpdir")
        if _tmp and os.path.isdir(_tmp):
            shutil.rmtree(_tmp, ignore_errors=True)
        st.session_state.clear()
        st.rerun()

st.title("🔬 NPLab - NanoPlastic Adsorption Analysis")


# ── ZIP helper ─────────────────────────────────────────────────────────────────

# Filename patterns used to classify files by tab (for the summary only)
_FILE_PATTERNS = {
    "Calibration": re.compile(r"(CCN|HCCN|LCCN)\d", re.IGNORECASE),
    "Kinetics":    re.compile(r"AKN\d",              re.IGNORECASE),
    "Isotherms":   re.compile(r"AIN\d",              re.IGNORECASE),
    "Drift":       re.compile(r"DRF\d",              re.IGNORECASE),
    "Control":     re.compile(r"CTRL\d",             re.IGNORECASE),
    "Blank":       re.compile(r"(BLK|BLANK)\d",      re.IGNORECASE),
}
_SPECTRAL_EXT = (".txt", ".asc", ".akn")


def _extract_zip(uploaded_zip) -> tuple[str, dict[str, int]]:
    """
    Extract *uploaded_zip* (a Streamlit UploadedFile) into a fresh temp
    directory, flattening any sub-folder structure so all spectral files
    end up in the root of that directory.

    Returns:
        tmpdir  – absolute path to the temp directory
        summary – {tab_name: file_count} for the sidebar status badge
    """
    tmpdir = tempfile.mkdtemp(prefix="nplab_zip_")
    with zipfile.ZipFile(uploaded_zip, "r") as zf:
        for member in zf.infolist():
            # Skip directories and macOS/Windows metadata entries
            if member.is_dir():
                continue
            fname = os.path.basename(member.filename)
            if not fname or fname.startswith((".", "__")):
                continue
            if not fname.lower().endswith(_SPECTRAL_EXT):
                continue
            dest = os.path.join(tmpdir, fname)
            # Avoid overwriting if two sub-folders contain a same-named file
            if os.path.exists(dest):
                base, ext = os.path.splitext(fname)
                counter = 1
                while os.path.exists(dest):
                    dest = os.path.join(tmpdir, f"{base}_{counter}{ext}")
                    counter += 1
            with zf.open(member) as src, open(dest, "wb") as dst:
                dst.write(src.read())

    # Build summary
    all_files = [
        f for f in os.listdir(tmpdir)
        if f.lower().endswith(_SPECTRAL_EXT)
    ]
    summary: dict[str, int] = {}
    for tab, pat in _FILE_PATTERNS.items():
        n = sum(1 for f in all_files if pat.search(f))
        if n:
            summary[tab] = n
    unclassified = len(all_files) - sum(summary.values())
    if unclassified:
        summary["Other"] = unclassified

    return tmpdir, summary


# ── Sidebar ────────────────────────────────────────────────────────────────────

current_user = st.session_state.get("username", "")
default_folder = "data-samples"

with st.sidebar:
    st.header("⚙️ Data Source")

    data_source = st.radio(
        "Input mode",
        ["📁 Folder", "🗜️ ZIP file"],
        horizontal=True,
        key="data_source_mode",
    )

    # ── Mode A: folder picker (original behaviour) ─────────────────────────
    if data_source == "📁 Folder":
        # If we're switching back from ZIP mode, clean up the old temp dir
        _old_tmp = st.session_state.pop("_zip_tmpdir", None)
        if _old_tmp and os.path.isdir(_old_tmp):
            shutil.rmtree(_old_tmp, ignore_errors=True)
        st.session_state.pop("_zip_summary", None)

        EXPERIMENT_DIR = folder_picker("Experiment data folder", default_folder, "exp")

    # ── Mode B: ZIP upload ─────────────────────────────────────────────────
    else:
        uploaded_zip = st.file_uploader(
            "Upload experiment ZIP",
            type="zip",
            key="zip_uploader",
            help=(
                "Pack all your spectral files (.txt / .asc / .akn) "
                "into a single ZIP. Sub-folders are supported — "
                "files will be flattened automatically."
            ),
        )

        if uploaded_zip is not None:
            # Re-extract only when a new file is uploaded (compare by name+size)
            zip_fingerprint = (uploaded_zip.name, uploaded_zip.size)
            if st.session_state.get("_zip_fingerprint") != zip_fingerprint:
                # Remove previous temp dir if any
                _old_tmp = st.session_state.get("_zip_tmpdir")
                if _old_tmp and os.path.isdir(_old_tmp):
                    shutil.rmtree(_old_tmp, ignore_errors=True)

                with st.spinner("Extracting ZIP…"):
                    tmpdir, summary = _extract_zip(uploaded_zip)

                st.session_state["_zip_tmpdir"]      = tmpdir
                st.session_state["_zip_summary"]     = summary
                st.session_state["_zip_fingerprint"] = zip_fingerprint

            EXPERIMENT_DIR = st.session_state["_zip_tmpdir"]

            # Show a concise file-count badge per tab
            summary = st.session_state.get("_zip_summary", {})
            if summary:
                st.success(f"✅ **{uploaded_zip.name}** loaded")
                total = sum(summary.values())
                st.caption(f"{total} spectral file(s) found.")
            else:
                st.warning("No recognised spectral files found in ZIP.")

        else:
            # No ZIP uploaded yet — fall back to the default folder so tabs
            # don't crash, and inform the user.
            EXPERIMENT_DIR = default_folder
            st.info("Upload a ZIP file to load your experiment data.")


# ── Tabs ───────────────────────────────────────────────────────────────────────
tab_viz, tab_cal, tab_blk, tab_drf, tab_ctrl, tab_kin, tab_iso, tab_dl, tab_help = st.tabs(
    ["Visualization", "Calibration", "Blank", "Drift", "Control",
     "Kinetics", "Isotherms", "Download Results", "Help"]
)

with tab_viz:
    render_visualization(EXPERIMENT_DIR)

with tab_cal:
    render_calibration(EXPERIMENT_DIR)

with tab_blk:
    render_blank(EXPERIMENT_DIR)

with tab_drf:
    render_drift(EXPERIMENT_DIR)

with tab_ctrl:
    render_control(EXPERIMENT_DIR)

with tab_kin:
    render_kinetics(EXPERIMENT_DIR)

with tab_iso:
    render_isotherms(EXPERIMENT_DIR)

with tab_dl:
    render_download()

with tab_help:
    render_help()
