#!/usr/bin/env python3
"""
gutMAP Interactive Pipeline Dashboard
======================================
Streamlit app that loads a single HJH + merged trace file, lets the user
select batches from the trace's ``batch`` column, and generates a live report.

Usage:
    streamlit run analysis/dashboard.py
"""

import io
import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from figures import (
    FULL_DATASET,
    STEP_COLOR,
    STEP_ORDER,
    build_html,
    build_size_cost_model,
    extrapolate_size_weighted,
    fig_billed_bar,
    fig_billed_pie,
    fig_billed_vs_cpu,
    fig_concurrency,
    fig_cpu_share,
    fig_cpu_share_pie,
    fig_efficiency_box,
    fig_efficiency_weighted,
    fig_io,
    fig_mem_waste,
    fig_memory_box,
    fig_model_diagnostics_table,
    fig_outlier_table,
    fig_per_sample_cpu,
    fig_rss_box,
    fig_scaling,
    fig_scaling_size_weighted,
    fig_size_distribution,
    fig_size_vs_cost_scatter,
    fig_status_overall,
    fig_status_per_step,
    fig_summary_table,
    load_full_sizes,
    load_hjh,
    load_trace,
    _read_trace_csv,
)

DATA = Path(__file__).resolve().parents[2] / "data"  # repo/data

# --- Defaults (latest mega-run) ---
_DEFAULT_DIR = DATA / "raw"
DEFAULT_HJH = _DEFAULT_DIR / "hpc-jobs-history.csv"
DEFAULT_TRACE = max(_DEFAULT_DIR.glob("merged_*.csv"), default=None, key=lambda p: p.name)


# ---------------------------------------------------------------------------
# Loaders (cached)
# ---------------------------------------------------------------------------

@st.cache_data
def cached_load_hjh(path: str) -> pd.DataFrame:
    return load_hjh(Path(path))


@st.cache_data
def cached_load_hjh_bytes(data: bytes, name: str) -> pd.DataFrame:
    return load_hjh(io.BytesIO(data))


@st.cache_data
def cached_load_trace(path: str) -> pd.DataFrame:
    return load_trace(Path(path))


@st.cache_data
def cached_load_trace_bytes(data: bytes, name: str) -> pd.DataFrame:
    return load_trace(io.BytesIO(data))


@st.cache_data
def _parse_trace_all(data) -> pd.DataFrame:
    """Load all trace rows keeping batch, native_id, name, status.
    CACHED rows are dropped: they reference jobs from prior runs and would
    distort cost/efficiency stats and status charts for the current run."""
    df = _read_trace_csv(data)
    df = df[df["status"] != "CACHED"].copy()
    df["pipeline_step"] = df["name"].str.split(r"\s*\(", n=1).str[0].str.strip()
    df["sra_id"] = df["name"].str.extract(r"([A-Z]{3}\d+)\)?$")
    df = df[df["pipeline_step"].isin(STEP_ORDER)].copy()
    df["pipeline_step"] = pd.Categorical(
        df["pipeline_step"], categories=STEP_ORDER, ordered=True,
    )
    return df


@st.cache_data
def cached_load_trace_all(path: str) -> pd.DataFrame:
    return _parse_trace_all(path)


@st.cache_data
def cached_load_trace_all_bytes(data: bytes, name: str) -> pd.DataFrame:
    return _parse_trace_all(io.BytesIO(data))


@st.cache_data
def _parse_batch_map(data) -> pd.DataFrame:
    """Load just native_id + batch for filtering. Drops CACHED rows so that
    the (native_id -> batch) lookup never pulls a prior-run SLURM job into
    a later batch's filtered HJH (which would double-count its cost)."""
    df = _read_trace_csv(data)
    df = df[df["status"] != "CACHED"].copy()
    df = df[["native_id", "batch"]].copy()
    return df


@st.cache_data
def cached_load_batch_map(path: str) -> pd.DataFrame:
    return _parse_batch_map(path)


@st.cache_data
def cached_load_batch_map_bytes(data: bytes, name: str) -> pd.DataFrame:
    return _parse_batch_map(io.BytesIO(data))


@st.cache_data
def cached_load_full_sizes() -> pd.DataFrame | None:
    if FULL_DATASET.exists():
        return load_full_sizes(FULL_DATASET)
    return None


def _batch_sort_key(name: str) -> int:
    """Extract start number from batch_s123-456 for natural sort."""
    m = re.search(r"batch_s(\d+)", name)
    return int(m.group(1)) if m else 0


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def filter_by_batches(
    hjh: pd.DataFrame,
    trace: pd.DataFrame,
    trace_all: pd.DataFrame,
    batch_map: pd.DataFrame,
    selected: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    keep_native = (
        batch_map[batch_map["batch"].isin(selected)]["native_id"]
        .dropna().astype(str).str.strip()
    )
    keep_set = set(keep_native)

    trace_f = trace[trace["native_id"].astype(str).str.strip().isin(keep_set)].copy()
    hjh_f = hjh[hjh["job_id"].astype(str).isin(keep_set)].copy()

    if "native_id" in trace_all.columns:
        trace_all_f = trace_all[
            trace_all["native_id"].astype(str).str.strip().isin(keep_set)
        ].copy()
    else:
        # Fallback: filter by name matching
        trace_all_f = trace_all[
            trace_all["name"].isin(trace_f["name"])
        ].copy()

    return hjh_f, trace_f, trace_all_f


# ---------------------------------------------------------------------------
# Extra figures (dashboard-only)
# ---------------------------------------------------------------------------

def fig_cpu_concurrency(hjh: pd.DataFrame) -> go.Figure:
    """Total CPUs in use over time (similar to fig_concurrency but sums cores)."""
    df = hjh[["pipeline_step", "end_time", "wall_used[s]", "cores"]].dropna().copy()
    df["start_time"] = df["end_time"] - df["wall_used[s]"]
    t_min, t_max = df["start_time"].min(), df["end_time"].max()

    resolution = max(1, int((t_max - t_min) / 500))
    starts = df["start_time"].values
    ends = df["end_time"].values
    cores = df["cores"].values

    times = np.arange(int(t_min), int(t_max), resolution)
    cpu_counts = np.zeros(len(times), dtype=np.float64)
    for i, t in enumerate(times):
        mask = (starts <= t) & (ends >= t)
        cpu_counts[i] = cores[mask].sum()

    hours_offset = (times - t_min) / 3600
    fig = go.Figure(go.Scatter(
        x=hours_offset.tolist(), y=cpu_counts.tolist(),
        mode="lines", fill="tozeroy", line_color="#e67e22",
    ))
    fig.update_layout(
        title="CPU cores in use over time (all steps combined)",
        xaxis_title="Elapsed time (h)", yaxis_title="Total CPU cores",
        height=350,
    )
    return fig


# ---------------------------------------------------------------------------
# KPI + figures
# ---------------------------------------------------------------------------

def compute_kpis(
    hjh: pd.DataFrame,
    size_model_available: bool = False,
    size_weighted_total: float = 0,
    pct_extrap: float = 0,
    n_target_size: int = 0,
) -> list[tuple[str, str, str]]:
    n_samples = hjh["sra_id"].nunique()
    total_cpu_h = hjh["cpu_hours"].sum()
    total_billed_h = hjh["billed_core_hours"].sum()
    total_wall_h = hjh["wall_hours"].sum()
    n_jobs = len(hjh)
    cpu_sum = hjh["cpu_hours"].sum()
    mean_eff = (
        (hjh["efficiency[%]"] * hjh["cpu_hours"]).sum() / cpu_sum
        if cpu_sum > 0 else 0
    )

    kpis = [
        ("Samples", f"{n_samples:,}", ""),
        ("Total jobs", f"{n_jobs:,}", "SLURM jobs"),
        ("Billed core-hours", f"{total_billed_h:,.0f}", "TRESBillingWeights"),
        ("Total core-hours", f"{total_cpu_h:,.0f}", "CPU used"),
        ("Total wall hours", f"{total_wall_h:,.0f}", "all jobs summed"),
        ("Weighted mean efficiency", f"{mean_eff:.1f}%", "CPU eff"),
        (
            "Projected billed core-h @ 100k",
            f"{total_billed_h / n_samples * 100_000:,.0f}" if n_samples > 0 else "N/A",
            "linear extrapolation",
        ),
    ]
    if size_model_available:
        kpis.append((
            f"Size-weighted projected billed core-h ({n_target_size:,} samples)",
            f"{size_weighted_total:,.0f}",
            f"{pct_extrap:.1f}% from extrapolation",
        ))
    return kpis


def build_all_figures(
    hjh: pd.DataFrame,
    trace: pd.DataFrame | None,
    trace_all: pd.DataFrame | None,
    size_model_available: bool = False,
    models_df: pd.DataFrame | None = None,
    hjh_with_size: pd.DataFrame | None = None,
    extrap_df: pd.DataFrame | None = None,
    full_sizes: pd.DataFrame | None = None,
    n_target_size: int = 0,
    bins_df: pd.DataFrame | None = None,
) -> dict[str, go.Figure]:
    n_samples = hjh["sra_id"].nunique()
    figs = {
        "summary_table": fig_summary_table(hjh),
        "billed_bar": fig_billed_bar(hjh),
        "billed_pie": fig_billed_pie(hjh),
        "billed_vs_cpu": fig_billed_vs_cpu(hjh),
        "cpu_share_bar": fig_cpu_share(hjh),
        "cpu_share_pie": fig_cpu_share_pie(hjh),
        "efficiency_box": fig_efficiency_box(hjh),
        "efficiency_weighted": fig_efficiency_weighted(hjh),
        "memory_box": fig_memory_box(hjh),
        "mem_waste": fig_mem_waste(hjh),
        "concurrency": fig_concurrency(hjh),
        "cpu_concurrency": fig_cpu_concurrency(hjh),
        "per_sample_cpu": fig_per_sample_cpu(hjh),
        "outlier_table": fig_outlier_table(hjh),
        "scaling": fig_scaling(hjh, n_samples),
    }

    if trace is not None and len(trace) > 0:
        figs["io"] = fig_io(trace)
        figs["rss_box"] = fig_rss_box(trace)
    else:
        figs["io"] = go.Figure().update_layout(title="I/O — trace file not available")
        figs["rss_box"] = go.Figure().update_layout(title="RSS — trace file not available")

    if trace_all is not None and len(trace_all) > 0:
        figs["status_per_step"] = fig_status_per_step(trace_all)
        figs["status_overall"] = fig_status_overall(trace_all)
    else:
        figs["status_per_step"] = go.Figure().update_layout(
            title="Job status — trace not available")
        figs["status_overall"] = go.Figure().update_layout(
            title="Job status — trace not available")

    if size_model_available:
        figs["scaling_size_weighted"] = fig_scaling_size_weighted(
            hjh, extrap_df, n_samples, n_target_size)
        figs["model_diagnostics"] = fig_model_diagnostics_table(models_df, extrap_df)
        figs["size_vs_cost_scatter"] = fig_size_vs_cost_scatter(hjh_with_size, models_df, bins_df)
        figs["size_distribution"] = fig_size_distribution(hjh_with_size, full_sizes)

    return figs


# ---------------------------------------------------------------------------
# Streamlit UI
# ---------------------------------------------------------------------------

st.set_page_config(page_title="gutMAP Dashboard", layout="wide")
st.title("gutMAP Pipeline Performance Report")

# ---- Sidebar: file pickers (compact) + batch selector ----
with st.sidebar:
    # File pickers — small popover buttons with icon
    fc1, fc2 = st.columns(2)
    with fc1.popover("\U0001F4C2 HJH"):
        hjh_upload = st.file_uploader("HJH CSV", type=["csv"], key="hjh_up",
                                      label_visibility="collapsed")
    with fc2.popover("\U0001F4C2 Trace"):
        trace_upload = st.file_uploader("Trace CSV", type=["csv"], key="trace_up",
                                        label_visibility="collapsed")

    # Resolve data sources: uploaded files take priority over defaults
    hjh_source = trace_source = None
    hjh_label = trace_label = ""

    if hjh_upload is not None:
        hjh_source = ("upload", hjh_upload.getvalue(), hjh_upload.name)
        hjh_label = hjh_upload.name
    elif DEFAULT_HJH and DEFAULT_HJH.exists():
        hjh_source = ("path", str(DEFAULT_HJH), None)
        hjh_label = DEFAULT_HJH.name

    if trace_upload is not None:
        trace_source = ("upload", trace_upload.getvalue(), trace_upload.name)
        trace_label = trace_upload.name
    elif DEFAULT_TRACE and DEFAULT_TRACE.exists():
        trace_source = ("path", str(DEFAULT_TRACE), None)
        trace_label = DEFAULT_TRACE.name

    if not hjh_source:
        st.error("No HJH file")
        st.stop()
    if not trace_source:
        st.error("No trace file")
        st.stop()

    st.caption(f"{hjh_label}  \n{trace_label}")

    # Load data (cached)
    if hjh_source[0] == "upload":
        hjh_full = cached_load_hjh_bytes(hjh_source[1], hjh_source[2])
    else:
        hjh_full = cached_load_hjh(hjh_source[1])

    if trace_source[0] == "upload":
        trace_full = cached_load_trace_bytes(trace_source[1], trace_source[2])
        trace_all_full = cached_load_trace_all_bytes(trace_source[1], trace_source[2])
        batch_map = cached_load_batch_map_bytes(trace_source[1], trace_source[2])
    else:
        trace_full = cached_load_trace(trace_source[1])
        trace_all_full = cached_load_trace_all(trace_source[1])
        batch_map = cached_load_batch_map(trace_source[1])

    st.divider()
    st.header("Batches")

    if "batch" in batch_map.columns:
        all_batches = sorted(
            batch_map["batch"].dropna().unique().tolist(),
            key=_batch_sort_key,
        )
    else:
        all_batches = ["all"]

    col_a, col_b = st.columns(2)
    if col_a.button("Select all"):
        st.session_state["batch_selector"] = all_batches
    if col_b.button("Clear all"):
        st.session_state["batch_selector"] = []

    selected = st.multiselect(
        "Select batches to include",
        options=all_batches,
        default=all_batches[:1],
        key="batch_selector",
    )

    st.divider()
    generate = st.button("Generate report", type="primary", use_container_width=True)


# ---- Main area ----
if generate and selected:
    with st.spinner("Filtering and computing..."):
        # Filter by selected batches
        if all_batches == ["all"]:
            hjh, trace, trace_all = hjh_full, trace_full, trace_all_full
        else:
            hjh, trace, trace_all = filter_by_batches(
                hjh_full, trace_full, trace_all_full, batch_map, selected,
            )

        if len(hjh) == 0:
            st.warning("No HJH rows match the selected batches.")
            st.stop()

        st.caption(
            f"Showing {hjh['sra_id'].nunique():,} samples, "
            f"{len(hjh):,} jobs from {len(selected)} batch(es)"
        )

        # Size-weighted model
        full_sizes = cached_load_full_sizes()
        size_model_available = False
        models_df = hjh_with_size = extrap_df = bins_df = None
        n_target_size = 0
        size_weighted_total = pct_extrap = 0

        if full_sizes is not None:
            n_target_size = full_sizes["sra_id"].nunique()
            models_df, hjh_with_size, bins_df = build_size_cost_model(hjh, full_sizes)
            if len(models_df) > 0:
                extrap_df = extrapolate_size_weighted(models_df, full_sizes, hjh)
                size_weighted_total = extrap_df["total_predicted"].sum()
                pct_extrap = (
                    extrap_df["extrapolated_cost"].sum()
                    / extrap_df["total_predicted"].sum() * 100
                )
                size_model_available = True

        # KPIs
        kpis = compute_kpis(
            hjh, size_model_available, size_weighted_total, pct_extrap, n_target_size,
        )

        # Render KPI cards — split into rows of 4
        for row_start in range(0, len(kpis), 4):
            row_kpis = kpis[row_start : row_start + 4]
            kpi_cols = st.columns(len(row_kpis))
            for col, (label, value, sub) in zip(kpi_cols, row_kpis):
                col.metric(label=label, value=value, help=sub if sub else None)

        # Build figures
        figs = build_all_figures(
            hjh, trace, trace_all,
            size_model_available, models_df, hjh_with_size,
            extrap_df, full_sizes, n_target_size, bins_df,
        )

    # ---- Section 1: Overview ----
    st.header("1 — Overview")
    st.plotly_chart(figs["summary_table"], use_container_width=True)

    # ---- Section 2: Job completion status ----
    st.header("2 — Job completion status")
    c1, c2 = st.columns(2)
    c1.plotly_chart(figs["status_per_step"], use_container_width=True)
    c2.plotly_chart(figs["status_overall"], use_container_width=True)

    # ---- Section 3: Billed core-hours ----
    st.header("3 — Billed core-hours")
    c1, c2 = st.columns(2)
    c1.plotly_chart(figs["billed_bar"], use_container_width=True)
    c2.plotly_chart(figs["billed_pie"], use_container_width=True)
    st.plotly_chart(figs["billed_vs_cpu"], use_container_width=True)

    # ---- Section 4: Core-hours (CPU used) ----
    st.header("4 — Core-hours (CPU used)")
    c1, c2 = st.columns(2)
    c1.plotly_chart(figs["cpu_share_bar"], use_container_width=True)
    c2.plotly_chart(figs["cpu_share_pie"], use_container_width=True)

    # ---- Section 5: CPU efficiency ----
    st.header("5 — CPU efficiency")
    st.plotly_chart(figs["efficiency_box"], use_container_width=True)
    st.plotly_chart(figs["efficiency_weighted"], use_container_width=True)

    # ---- Section 6: Memory ----
    st.header("6 — Memory")
    st.plotly_chart(figs["memory_box"], use_container_width=True)
    st.plotly_chart(figs["mem_waste"], use_container_width=True)

    # ---- Section 7: I/O and RSS ----
    st.header("7 — I/O and RSS (Nextflow trace)")
    c1, c2 = st.columns(2)
    c1.plotly_chart(figs["io"], use_container_width=True)
    c2.plotly_chart(figs["rss_box"], use_container_width=True)

    # ---- Section 8: Concurrency ----
    st.header("8 — Concurrency")
    st.plotly_chart(figs["concurrency"], use_container_width=True)
    st.plotly_chart(figs["cpu_concurrency"], use_container_width=True)

    # ---- Section 9: Sample-level variability ----
    st.header("9 — Sample-level variability")
    c1, c2 = st.columns(2)
    c1.plotly_chart(figs["per_sample_cpu"], use_container_width=True)
    c2.plotly_chart(figs["outlier_table"], use_container_width=True)

    # ---- Section 10: 100k scaling (count-based) ----
    st.header("10 — 100k scaling projection (count-based)")
    st.plotly_chart(figs["scaling"], use_container_width=True)

    # ---- Section 11: Size-weighted scaling (conditional) ----
    if size_model_available:
        st.header("11 — Size-weighted scaling projection")
        st.caption(
            "Projects costs to the full target dataset by modeling per-step cost as a "
            "linear function of input file size (cost = a + b * size_GB), "
            "then applying each step's model to every sample in the full dataset. "
            "Costs predicted outside the observed size range are flagged as extrapolation."
        )
        st.plotly_chart(figs["scaling_size_weighted"], use_container_width=True)
        st.plotly_chart(figs["model_diagnostics"], use_container_width=True)
        c1, c2 = st.columns(2)
        c1.plotly_chart(figs["size_vs_cost_scatter"], use_container_width=True)
        c2.plotly_chart(figs["size_distribution"], use_container_width=True)

    # --- HTML Export (removable) ---
    html_report = build_html(figs, kpis)
    st.sidebar.download_button(
        "Export as HTML",
        data=html_report,
        file_name="pipeline_report.html",
        mime="text/html",
    )
    # --- End HTML Export ---

elif not selected:
    st.info("Select at least one batch from the sidebar, then click **Generate report**.")
