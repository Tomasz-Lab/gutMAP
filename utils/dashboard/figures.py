"""
gutMAP dashboard figure library.

Loaders, cost models, and Plotly figure builders consumed by dashboard.py
(and build_html, used by the dashboard's HTML-export button). Not an entry
point — import from it, don't run it.
"""

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE = Path(__file__).parent
ROOT = HERE.parents[1]            # repo root (utils/dashboard -> utils -> repo)
DATA = ROOT / "data"
FULL_DATASET  = DATA / "input" / "input_dataset.csv"

STEP_ORDER = [
    "CHECK_AND_CONVERT_READS",
    "FASTP_QC",
    "REMOVE_HUMAN_READS",
    "MEGAHIT_ASSEMBLY",
    "SYLPH_TAXONOMY",
    "BAKTA_ANNOTATION",
    "BAKTA_BATCH",
    "CREATE_CATALOG_AND_INDEX",
    "ALIGN_AND_QUANTIFY_READS",
    "CALCULATE_TPM_AND_ANNOTATE",
    "MAP_FOR_BINNING",
    "METABAT2_BINNING",
    "CHECKM_QA",
    "FILTER_AND_ANNOTATE",
]

# Steps where one SLURM job processes multiple samples.
# Used to compute correct per-sample cost for scaling / extrapolation.
BATCH_MULTIPLIER = {"BAKTA_BATCH": 20}

COLORS = px.colors.qualitative.Bold
STEP_COLOR = {s: COLORS[i % len(COLORS)] for i, s in enumerate(STEP_ORDER)}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_size(val: str) -> float:
    """'1.4 GB' -> bytes as float. Returns NaN on failure."""
    if pd.isna(val):
        return float("nan")
    val = str(val).strip()
    m = re.match(r"^([0-9.]+)\s*([A-Za-z]+)$", val)
    if not m:
        return float("nan")
    num, unit = float(m.group(1)), m.group(2).upper()
    units = {"B": 1, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12,
             "KIB": 1024, "MIB": 1024**2, "GIB": 1024**3, "TIB": 1024**4}
    return num * units.get(unit, float("nan"))


def _parse_duration(val: str) -> float:
    """'1m 3s' / '2.4s' / '1h 2m 3s' -> seconds as float."""
    if pd.isna(val):
        return float("nan")
    val = str(val).strip()
    total = 0.0
    for amount, unit in re.findall(r"([0-9.]+)\s*([a-z]+)", val):
        amount = float(amount)
        if unit.startswith("d"):
            total += amount * 86400
        elif unit.startswith("h"):
            total += amount * 3600
        elif unit.startswith("m"):
            total += amount * 60
        elif unit.startswith("s"):
            total += amount
        elif unit.startswith("ms"):
            total += amount / 1000
    return total if total > 0 else float("nan")


PARTITION_MEM_WEIGHT = {
    "plgrid":        0.0005,    # cpu=1,mem=0.0005M  (per MB declared)
    "plgrid-bigmem": 0.00025,   # cpu=1,mem=0.00025M
}
DEFAULT_MEM_WEIGHT = 0.0005


def load_hjh(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Drop non-pipeline steps (e.g. 'bash' wrapper jobs)
    df = df[df["pipeline_step"].isin(STEP_ORDER)].copy()
    df["pipeline_step"] = pd.Categorical(df["pipeline_step"], categories=STEP_ORDER, ordered=True)
    # Coerce to numeric — '--' placeholders become NaN
    for col in ("mem_usage[%]", "efficiency[%]", "cpu_used[s]", "wall_used[s]"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["cpu_hours"]    = df["cpu_used[s]"] / 3600
    df["wall_hours"]   = df["wall_used[s]"] / 3600
    df["decl_mem_GB"]  = df["decl_mem[B]"] / 1e9
    df["actual_mem_GB"] = df["decl_mem_GB"] * df["mem_usage[%]"] / 100
    df["wasted_mem_GB"] = df["decl_mem_GB"] - df["actual_mem_GB"]
    # Billed core-hours: TRESBillingWeights cpu=1, mem weight per partition
    mem_weight = df["partition"].map(PARTITION_MEM_WEIGHT).fillna(DEFAULT_MEM_WEIGHT)
    cpu_billing  = df["cores"]
    mem_billing  = df["decl_mem[B]"] / 1e6 * mem_weight   # bytes → MB, then weight
    df["billing_rate"]     = cpu_billing.where(cpu_billing >= mem_billing, mem_billing)
    df["billed_core_hours"] = df["billing_rate"] * df["wall_used[s]"] / 3600
    return df


def _read_trace_csv(path, **kwargs) -> pd.DataFrame:
    """Read a Nextflow trace CSV that may have unquoted commas in the name column.

    BAKTA_BATCH rows list sra_ids in the name field (e.g. "BAKTA_BATCH (... SRR1, SRR2, ...)"),
    producing more fields than the header. This reader detects that and merges the
    overflow back into the name column.
    """
    import csv, io

    # Read raw text
    if isinstance(path, (str, Path)):
        text = Path(path).read_text(encoding="utf-8")
    else:  # BytesIO / file-like
        raw = path.read()
        text = raw.decode("utf-8") if isinstance(raw, bytes) else raw
        if hasattr(path, "seek"):
            path.seek(0)

    reader = csv.reader(io.StringIO(text))
    header = next(reader)
    n_cols = len(header)
    name_idx = header.index("name")

    rows = []
    for fields in reader:
        if len(fields) > n_cols:
            # Merge overflow fields back into name
            overflow = n_cols - len(fields)  # negative
            name_parts = fields[name_idx : name_idx + 1 - overflow]
            merged_name = ",".join(name_parts)
            fields = fields[:name_idx] + [merged_name] + fields[name_idx + 1 - overflow :]
        rows.append(fields)

    return pd.DataFrame(rows, columns=header, **kwargs)


def load_trace(path: Path) -> pd.DataFrame:
    df = _read_trace_csv(path)
    df = df[df["status"] == "COMPLETED"].copy()

    # Parse string columns
    df["cpu_pct"]       = pd.to_numeric(df["%cpu"].str.replace("%", "", regex=False), errors="coerce")
    df["peak_rss_GB"]   = df["peak_rss"].apply(_parse_size) / 1e9
    df["peak_vmem_GB"]  = df["peak_vmem"].apply(_parse_size) / 1e9
    df["rchar_GB"]      = df["rchar"].apply(_parse_size) / 1e9
    df["wchar_GB"]      = df["wchar"].apply(_parse_size) / 1e9
    df["realtime_s"]    = df["realtime"].apply(_parse_duration)
    df["duration_s"]    = df["duration"].apply(_parse_duration)

    # Extract step name from trace 'name' column: "FASTP_QC (fastp: ERR123)"
    df["pipeline_step"] = df["name"].str.split(r"\s*\(", n=1).str[0].str.strip()
    df["sra_id"]        = df["name"].str.extract(r"([A-Z]{3}\d+)\)?$")
    df = df[df["pipeline_step"].isin(STEP_ORDER)].copy()
    df["pipeline_step"] = pd.Categorical(df["pipeline_step"], categories=STEP_ORDER, ordered=True)
    return df


def load_full_sizes(path: Path) -> pd.DataFrame:
    """Load full dataset with file sizes. Returns DataFrame with sra_id and size_GB."""
    df = pd.read_csv(path)
    df["size_GB"] = pd.to_numeric(df["Bytes"], errors="coerce") / 1e9
    return df


def build_size_cost_model(
    hjh: pd.DataFrame, full_sizes: pd.DataFrame, n_bins: int = 20,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Per-step binned regression: bin samples by file size, average cost per bin,
    then OLS on bin means: mean_billed_core_hours = intercept + slope * mean_size_GB.
    Returns (models_df, hjh_with_size, bins_df).
    """
    hjh_with_size = hjh.merge(
        full_sizes[["sra_id", "size_GB"]], on="sra_id", how="left",
    )
    model_rows = []
    all_bins = []
    for step in STEP_ORDER:
        sub = hjh_with_size[
            (hjh_with_size["pipeline_step"] == step)
            & hjh_with_size["size_GB"].notna()
            & hjh_with_size["billed_core_hours"].notna()
        ]
        if len(sub) < 3:
            continue
        # Bin by size
        sub = sub.copy()
        sub["size_bin"] = pd.cut(sub["size_GB"], bins=n_bins)
        binned = sub.groupby("size_bin", observed=True).agg(
            mean_size_GB=("size_GB", "mean"),
            mean_cost=("billed_core_hours", "mean"),
            std_cost=("billed_core_hours", "std"),
            count=("billed_core_hours", "count"),
        ).dropna(subset=["mean_size_GB", "mean_cost"]).reset_index()
        binned["pipeline_step"] = step
        all_bins.append(binned)

        if len(binned) < 2:
            continue
        x = binned["mean_size_GB"].values
        y = binned["mean_cost"].values
        slope, intercept = np.polyfit(x, y, 1)
        y_pred = intercept + slope * x
        ss_res = ((y - y_pred) ** 2).sum()
        ss_tot = ((y - y.mean()) ** 2).sum()
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        model_rows.append({
            "pipeline_step": step,
            "intercept": intercept,
            "slope": slope,
            "r_squared": r2,
            "n_obs": len(sub),
            "n_bins": len(binned),
            "size_min_GB": float(sub["size_GB"].min()),
            "size_max_GB": float(sub["size_GB"].max()),
        })
    models_df = pd.DataFrame(model_rows)
    bins_df = pd.concat(all_bins, ignore_index=True) if all_bins else pd.DataFrame()
    return models_df, hjh_with_size, bins_df


def extrapolate_size_weighted(
    models_df: pd.DataFrame, full_sizes: pd.DataFrame,
    hjh: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Apply per-step binned-regression models to every sample in full_sizes.
    For batch steps (BATCH_MULTIPLIER), falls back to count-based extrapolation
    since they lack per-sample sra_ids for size regression.
    Returns per-step totals with interpolation vs extrapolation breakdown."""
    sizes = full_sizes[["sra_id", "size_GB"]].drop_duplicates(subset="sra_id").copy()
    n_target = len(sizes)
    median_size = sizes["size_GB"].median()
    sizes["size_GB"] = sizes["size_GB"].fillna(median_size)

    modeled_steps = set(models_df["pipeline_step"])

    results = []
    for _, model in models_df.iterrows():
        step = model["pipeline_step"]
        a, b = model["intercept"], model["slope"]
        predicted = (a + b * sizes["size_GB"]).clip(lower=0)
        in_range = (sizes["size_GB"] >= model["size_min_GB"]) & (
            sizes["size_GB"] <= model["size_max_GB"]
        )
        total = predicted.sum()
        interp = predicted[in_range].sum()
        extrap = predicted[~in_range].sum()
        results.append({
            "pipeline_step": step,
            "total_predicted": total,
            "interpolated_cost": interp,
            "extrapolated_cost": extrap,
            "pct_extrapolated": extrap / total * 100 if total > 0 else 0,
            "r_squared": model["r_squared"],
        })

    # Add count-based fallback for batch steps missing from size regression
    if hjh is not None:
        for step, mult in BATCH_MULTIPLIER.items():
            if step in modeled_steps:
                continue
            sub = hjh[hjh["pipeline_step"] == step]
            if len(sub) == 0:
                continue
            total_billed = sub["billed_core_hours"].sum()
            n_samples_processed = len(sub) * mult
            per_sample = total_billed / n_samples_processed if n_samples_processed > 0 else 0
            total_predicted = per_sample * n_target
            results.append({
                "pipeline_step": step,
                "total_predicted": total_predicted,
                "interpolated_cost": total_predicted,
                "extrapolated_cost": 0,
                "pct_extrapolated": 0,
                "r_squared": float("nan"),
            })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Plot builders
# ---------------------------------------------------------------------------

# Standard margins: enough room for angled tick labels (bottom) and
# outside-bar text labels (top).
_M_BAR  = dict(t=90,  b=160, l=70, r=30)   # bars with outside text + angled x
_M_BOX  = dict(t=70,  b=160, l=70, r=30)   # box/violin with angled x
_M_SCAT = dict(t=70,  b=60,  l=70, r=30)   # scatter (no angled x)
_M_HIST = dict(t=70,  b=60,  l=70, r=30)


def _no_clip(fig: go.Figure) -> go.Figure:
    """Prevent x-axis labels and outside-bar/scatter text from being clipped."""
    fig.update_xaxes(automargin=True)
    # cliponaxis is only valid for scatter and bar traces, not box/violin
    fig.update_traces(cliponaxis=False, selector={"type": "bar"})
    fig.update_traces(cliponaxis=False, selector={"type": "scatter"})
    return fig


def fig_cpu_share(hjh: pd.DataFrame) -> go.Figure:
    agg = hjh.groupby("pipeline_step", observed=True)["cpu_hours"].sum().reset_index()
    agg = agg.sort_values("cpu_hours", ascending=False)
    fig = px.bar(
        agg, x="pipeline_step", y="cpu_hours",
        color="pipeline_step", color_discrete_map=STEP_COLOR,
        labels={"cpu_hours": "Total core-hours", "pipeline_step": ""},
        title="Total core-hours per pipeline step",
        text=agg["cpu_hours"].round(1),
        height=520,
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(showlegend=False, xaxis_tickangle=-40, margin=_M_BAR)
    return _no_clip(fig)


def fig_cpu_share_pie(hjh: pd.DataFrame) -> go.Figure:
    agg = hjh.groupby("pipeline_step", observed=True)["cpu_hours"].sum().reset_index()
    fig = px.pie(
        agg, names="pipeline_step", values="cpu_hours",
        color="pipeline_step", color_discrete_map=STEP_COLOR,
        title="Core-hour share per step",
        hole=0.35,
        height=520,
    )
    fig.update_traces(textinfo="percent+label")
    fig.update_layout(margin=dict(t=70, b=30, l=30, r=30))
    return fig


def fig_efficiency_box(hjh: pd.DataFrame) -> go.Figure:
    fig = px.box(
        hjh.sort_values("pipeline_step"),
        x="pipeline_step", y="efficiency[%]",
        color="pipeline_step", color_discrete_map=STEP_COLOR,
        points="outliers",
        labels={"efficiency[%]": "CPU efficiency (%)", "pipeline_step": ""},
        title="CPU efficiency distribution per step (unweighted — shows spread)",
        height=520,
    )
    fig.update_layout(showlegend=False, xaxis_tickangle=-40, margin=_M_BOX)
    fig.add_hline(y=100, line_dash="dot", line_color="gray", annotation_text="100%")
    return _no_clip(fig)


def fig_efficiency_weighted(hjh: pd.DataFrame) -> go.Figure:
    """Where efficiency improvements pay off most.

    Allocated = cores * wall_h (== cpu_used[s]/3600 since SLURM's CPUTime field
    is elapsed*cores, not actual usage). Real CPU consumption is only encoded
    in efficiency[%], so:
        used_h   = allocated_h * efficiency[%] / 100
        wasted_h = allocated_h - used_h
    Sorted by wasted desc — top row = biggest absolute opportunity."""
    df = hjh.dropna(subset=["cores", "wall_used[s]", "efficiency[%]"]).copy()
    df["allocated_h"] = df["cores"] * df["wall_used[s]"] / 3600
    df["used_h"]      = df["allocated_h"] * df["efficiency[%]"] / 100
    df["wasted_h"]    = df["allocated_h"] - df["used_h"]

    agg = df.groupby("pipeline_step", observed=True).agg(
        allocated=("allocated_h", "sum"),
        used     =("used_h", "sum"),
        wasted   =("wasted_h", "sum"),
    ).reset_index()
    agg["weighted_eff_pct"] = (agg["used"] / agg["allocated"] * 100).fillna(0)
    # Drop categorical dtype so sort on a regular column works without surprises;
    # ascending order on y-axis means the biggest 'wasted' lands at the top of the bar.
    agg["pipeline_step"] = agg["pipeline_step"].astype(str)
    agg = agg.sort_values("wasted", ascending=True).reset_index(drop=True)
    y_order = agg["pipeline_step"].tolist()

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=agg["pipeline_step"], x=agg["used"], name="Used (cpu-h)",
        orientation="h", marker_color="#2ecc71",
        text=agg["used"].round(0), textposition="inside",
    ))
    fig.add_trace(go.Bar(
        y=agg["pipeline_step"], x=agg["wasted"], name="Wasted (cpu-h)",
        orientation="h", marker_color="#e74c3c",
        text=[f"save {w:,.0f} cpu-h ({100 - e:.0f}% waste, eff {e:.0f}%)"
              for e, w in zip(agg["weighted_eff_pct"], agg["wasted"])],
        textposition="outside",
    ))
    fig.update_layout(
        barmode="stack",
        title="Weighted CPU efficiency — where 1% improvement saves the most core-hours",
        xaxis_title="Core-hours (allocated = cores × wall)",
        yaxis_title="",
        yaxis=dict(categoryorder="array", categoryarray=y_order),
        height=520,
        margin=dict(t=70, b=60, l=200, r=240),
    )
    return _no_clip(fig)



def fig_memory_box(hjh: pd.DataFrame) -> go.Figure:
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=["Actual peak memory (GB)", "Memory utilisation (%)"])
    for i, step in enumerate(STEP_ORDER):
        sub = hjh[hjh["pipeline_step"] == step]
        color = STEP_COLOR.get(step, "#888")
        fig.add_trace(go.Box(y=sub["actual_mem_GB"], name=step,
                             marker_color=color, showlegend=False,
                             boxpoints="outliers"), row=1, col=1)
        fig.add_trace(go.Box(y=sub["mem_usage[%]"], name=step,
                             marker_color=color, showlegend=i == 0,
                             boxpoints="outliers"), row=1, col=2)
    fig.update_layout(title="Memory: actual usage and utilisation per step",
                      height=560, legend_title="Step", margin=_M_BOX)
    fig.update_xaxes(tickangle=-40, automargin=True)
    return fig


def fig_mem_waste(hjh: pd.DataFrame) -> go.Figure:
    agg = hjh.groupby("pipeline_step", observed=True).agg(
        declared=("decl_mem_GB", "mean"),
        actual=("actual_mem_GB", "mean"),
        wasted=("wasted_mem_GB", "mean"),
        max_actual=("actual_mem_GB", "max"),
    ).reset_index()
    fig = go.Figure()
    fig.add_trace(go.Bar(name="Actual (mean GB)", x=agg["pipeline_step"],
                         y=agg["actual"], marker_color="steelblue"))
    fig.add_trace(go.Bar(name="Wasted (mean GB)", x=agg["pipeline_step"],
                         y=agg["wasted"], marker_color="salmon"))
    fig.add_trace(go.Scatter(
        name="Max actual (GB)", x=agg["pipeline_step"], y=agg["max_actual"],
        mode="markers+text", marker=dict(symbol="diamond", size=10, color="#e74c3c"),
        text=agg["max_actual"].round(1), textposition="top center",
        textfont=dict(size=9, color="#e74c3c"),
    ))
    fig.update_layout(barmode="stack",
                      title="Declared vs actual memory per step (mean + max per job)",
                      xaxis_tickangle=-40, yaxis_title="Memory (GB)",
                      legend_title="", height=520, margin=_M_BAR)
    return _no_clip(fig)


def fig_per_sample_cpu(hjh: pd.DataFrame) -> go.Figure:
    per_sample = hjh.groupby("sra_id")["cpu_hours"].sum().reset_index()
    per_sample = per_sample.sort_values("cpu_hours", ascending=False)
    fig = px.histogram(
        per_sample, x="cpu_hours", nbins=60,
        labels={"cpu_hours": "Total core-hours per sample"},
        title="Per-sample total core-hours distribution",
        height=460,
    )
    q95 = per_sample["cpu_hours"].quantile(0.95)
    fig.add_vline(x=per_sample["cpu_hours"].median(), line_dash="dash",
                  annotation_text="median", line_color="steelblue")
    fig.add_vline(x=q95, line_dash="dot",
                  annotation_text="p95", line_color="tomato")
    fig.update_layout(margin=_M_HIST)
    return fig


def fig_scaling(hjh: pd.DataFrame, n_current: int, n_target: int = 100_000) -> go.Figure:
    scale = n_target / n_current
    agg = hjh.groupby("pipeline_step", observed=True).agg(
        cpu_hours=("cpu_hours", "sum"),
        n_jobs=("job_id", "count"),
    ).reset_index()
    # For batch steps, scale by actual samples processed (n_jobs * multiplier)
    def _step_scale(row):
        mult = BATCH_MULTIPLIER.get(row["pipeline_step"], 0)
        if mult:
            n_processed = row["n_jobs"] * mult
            return n_target / n_processed if n_processed > 0 else 0
        return scale
    agg["projected_cpu_hours"] = agg.apply(
        lambda r: r["cpu_hours"] * _step_scale(r), axis=1)
    agg = agg.sort_values("projected_cpu_hours", ascending=False)
    fig = px.bar(
        agg, x="pipeline_step", y="projected_cpu_hours",
        color="pipeline_step", color_discrete_map=STEP_COLOR,
        title=f"Projected core-hours at {n_target:,} samples (scaled from {n_current:,})",
        labels={"projected_cpu_hours": "Projected core-hours", "pipeline_step": ""},
        text=agg["projected_cpu_hours"].apply(lambda x: f"{x/1000:.1f}k"),
        height=520,
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(showlegend=False, xaxis_tickangle=-40, margin=_M_BAR)
    return _no_clip(fig)


def fig_io(trace: pd.DataFrame) -> go.Figure:
    agg = trace.groupby("pipeline_step", observed=True).agg(
        read_TB=("rchar_GB", lambda x: x.sum() / 1000),
        write_TB=("wchar_GB", lambda x: x.sum() / 1000),
    ).reset_index()
    fig = go.Figure()
    fig.add_trace(go.Bar(name="Read (TB total)", x=agg["pipeline_step"],
                         y=agg["read_TB"], marker_color="steelblue"))
    fig.add_trace(go.Bar(name="Write (TB total)", x=agg["pipeline_step"],
                         y=agg["write_TB"], marker_color="darkorange"))
    fig.update_layout(barmode="group", title="Total I/O per step (TB)",
                      xaxis_tickangle=-40, yaxis_title="TB",
                      legend_title="", height=520, margin=_M_BAR)
    return _no_clip(fig)


def fig_rss_box(trace: pd.DataFrame) -> go.Figure:
    fig = px.box(
        trace.sort_values("pipeline_step"),
        x="pipeline_step", y="peak_rss_GB",
        color="pipeline_step", color_discrete_map=STEP_COLOR,
        points="outliers",
        labels={"peak_rss_GB": "Peak RSS (GB)", "pipeline_step": ""},
        title="Peak RSS memory per step (from Nextflow trace)",
        height=520,
    )
    fig.update_layout(showlegend=False, xaxis_tickangle=-40, margin=_M_BOX)
    return _no_clip(fig)


def fig_concurrency(hjh: pd.DataFrame) -> go.Figure:
    """Compute running job count over time from end_time + wall_used[s]."""
    df = hjh[["pipeline_step", "end_time", "wall_used[s]"]].dropna().copy()
    df["start_time"] = df["end_time"] - df["wall_used[s]"]
    t_min, t_max = df["start_time"].min(), df["end_time"].max()

    resolution = max(1, int((t_max - t_min) / 500))
    times = range(int(t_min), int(t_max), resolution)
    counts = [((df["start_time"] <= t) & (df["end_time"] >= t)).sum() for t in times]

    hours_offset = [(t - t_min) / 3600 for t in times]
    fig = go.Figure(go.Scatter(x=hours_offset, y=counts, mode="lines",
                               fill="tozeroy", line_color="steelblue"))
    fig.update_layout(title="Pipeline concurrency over time (all steps combined)",
                      xaxis_title="Elapsed time (h)", yaxis_title="Running jobs",
                      height=350)
    return fig


def fig_summary_table(hjh: pd.DataFrame) -> go.Figure:
    agg = hjh.groupby("pipeline_step", observed=True).agg(
        n_jobs=("job_id", "count"),
        cpu_hours=("cpu_hours", "sum"),
        billed_core_hours=("billed_core_hours", "sum"),
        wall_hours=("wall_hours", "sum"),
        decl_cores=("cores", "median"),
        mean_mem_util=("mem_usage[%]", "mean"),
        mean_decl_mem_GB=("decl_mem_GB", "mean"),
        mean_actual_mem_GB=("actual_mem_GB", "mean"),
        _eff_x_cpu=("efficiency[%]", lambda x: (x * hjh.loc[x.index, "cpu_hours"]).sum()),
    ).reset_index()

    agg["mean_eff"] = agg["_eff_x_cpu"] / agg["cpu_hours"]
    agg = agg.drop(columns=["_eff_x_cpu"])

    total_cpu = agg["cpu_hours"].sum()
    agg["cpu_share_%"] = (agg["cpu_hours"] / total_cpu * 100).round(1)
    total_billed = agg["billed_core_hours"].sum()
    agg["billed_share_%"] = (agg["billed_core_hours"] / total_billed * 100).round(1)
    agg = agg.sort_values("billed_core_hours", ascending=False)

    fig = go.Figure(go.Table(
        header=dict(
            values=["Step", "Jobs", "Billed core-h", "Billed share %",
                    "Core-hours", "Core-hour share %",
                    "Weighted mean eff %", "Decl cores",
                    "Mean mem util %", "Decl mem (GB)", "Actual mem (GB)"],
            fill_color="steelblue", font_color="white", align="left",
        ),
        cells=dict(
            values=[
                agg["pipeline_step"],
                agg["n_jobs"],
                agg["billed_core_hours"].round(1),
                agg["billed_share_%"],
                agg["cpu_hours"].round(1),
                agg["cpu_share_%"],
                agg["mean_eff"].round(1),
                agg["decl_cores"].astype(int),
                agg["mean_mem_util"].round(1),
                agg["mean_decl_mem_GB"].round(1),
                agg["mean_actual_mem_GB"].round(1),
            ],
            align="left",
        )
    ))
    fig.update_layout(
        title="Per-step summary",
        height=450,
    )
    return fig


def fig_outlier_table(hjh: pd.DataFrame, top_n: int = 20) -> go.Figure:
    per_sample = hjh.groupby("sra_id").agg(
        total_cpu_h=("cpu_hours", "sum"),
        total_wall_h=("wall_hours", "sum"),
        n_steps=("pipeline_step", "nunique"),
    ).reset_index().sort_values("total_cpu_h", ascending=False).head(top_n)

    fig = go.Figure(go.Table(
        header=dict(
            values=["SRA ID", "Total core-hours", "Total wall hours", "Steps completed"],
            fill_color="steelblue", font_color="white", align="left",
        ),
        cells=dict(
            values=[
                per_sample["sra_id"],
                per_sample["total_cpu_h"].round(2),
                per_sample["total_wall_h"].round(2),
                per_sample["n_steps"],
            ],
            align="left",
        )
    ))
    fig.update_layout(title=f"Top {top_n} most expensive samples", height=550)
    return fig


def fig_billed_bar(hjh: pd.DataFrame) -> go.Figure:
    """Stacked bar: CPU base vs memory premium in billed core-hours per step.
    Memory premium only contributes when mem_rate > cores for that job;
    rows where cores dominate add zero to the premium total (clip lower=0)."""
    mem_weight = hjh["partition"].map(PARTITION_MEM_WEIGHT).fillna(DEFAULT_MEM_WEIGHT)
    mem_rate = hjh["decl_mem[B]"] / 1e6 * mem_weight
    cpu_base_h = hjh["cores"] * hjh["wall_used[s]"] / 3600
    mem_extra_h = ((mem_rate - hjh["cores"]).clip(lower=0)) * hjh["wall_used[s]"] / 3600

    tmp = hjh[["pipeline_step"]].copy()
    tmp["cpu_base"] = cpu_base_h
    tmp["mem_extra"] = mem_extra_h

    agg = tmp.groupby("pipeline_step", observed=True).agg(
        cpu_base=("cpu_base", "sum"),
        mem_extra=("mem_extra", "sum"),
    ).reset_index()
    agg["total"] = agg["cpu_base"] + agg["mem_extra"]
    agg = agg.sort_values("total", ascending=False)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=agg["pipeline_step"], y=agg["cpu_base"], name="CPU (cores)",
        marker_color="#3498db", text=agg["cpu_base"].round(0), textposition="inside",
    ))
    fig.add_trace(go.Bar(
        x=agg["pipeline_step"], y=agg["mem_extra"], name="Memory premium",
        marker_color="#e74c3c", text=agg["mem_extra"].round(0), textposition="inside",
    ))
    # Total label above each bar
    for _, row in agg.iterrows():
        fig.add_annotation(
            x=row["pipeline_step"], y=row["total"],
            text=f"{row['total']:.0f}",
            showarrow=False, yanchor="bottom", yshift=4,
            font=dict(size=10),
        )
    fig.update_layout(
        barmode="stack",
        title="Billed core-hours per step — CPU vs memory premium (mem weight per partition)",
        xaxis_tickangle=-40,
        yaxis_title="Billed core-hours",
        height=520,
        margin=_M_BAR,
    )
    return _no_clip(fig)


def fig_billed_pie(hjh: pd.DataFrame) -> go.Figure:
    """Pie chart of billed core-hour share per step. Groups <1% slices into Other."""
    agg = hjh.groupby("pipeline_step", observed=True)["billed_core_hours"].sum().reset_index()
    total = agg["billed_core_hours"].sum()
    agg["share"] = agg["billed_core_hours"] / total * 100
    # Group tiny slices into "Other"
    big = agg[agg["share"] >= 1.0].copy()
    small_total = agg.loc[agg["share"] < 1.0, "billed_core_hours"].sum()
    if small_total > 0:
        other = pd.DataFrame({"pipeline_step": ["Other"], "billed_core_hours": [small_total]})
        big = pd.concat([big, other], ignore_index=True)
    colors = [STEP_COLOR.get(s, "#aaa") for s in big["pipeline_step"]]
    fig = go.Figure(go.Pie(
        labels=big["pipeline_step"],
        values=big["billed_core_hours"],
        marker_colors=colors,
        hole=0.35,
        textinfo="percent+label",
    ))
    fig.update_layout(
        title="Billed core-hour share per step",
        height=520,
        margin=dict(t=70, b=30, l=30, r=30),
    )
    return fig


def fig_billed_vs_cpu(hjh: pd.DataFrame) -> go.Figure:
    """Grouped bar comparing billed core-hours vs actual core-hours per step."""
    agg = hjh.groupby("pipeline_step", observed=True).agg(
        billed=("billed_core_hours", "sum"),
        actual=("cpu_hours", "sum"),
    ).reset_index()
    agg["overhead_%"] = ((agg["billed"] - agg["actual"]) / agg["actual"] * 100).round(1)
    agg = agg.sort_values("billed", ascending=False)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=agg["pipeline_step"], y=agg["billed"], name="Billed core-h",
        marker_color="#e67e22", text=agg["billed"].round(0), textposition="outside",
    ))
    fig.add_trace(go.Bar(
        x=agg["pipeline_step"], y=agg["actual"], name="Actual core-h (CPU used)",
        marker_color="#3498db", text=agg["actual"].round(0), textposition="outside",
    ))
    # Annotate overhead %
    for _, row in agg.iterrows():
        fig.add_annotation(
            x=row["pipeline_step"], y=max(row["billed"], row["actual"]),
            text=f"+{row['overhead_%']}%" if row["overhead_%"] > 0 else f"{row['overhead_%']}%",
            showarrow=False, yanchor="bottom", yshift=20,
            font=dict(size=9, color="#888"),
        )
    fig.update_layout(
        barmode="group",
        title="Billed vs actual core-hours per step (overhead %)",
        xaxis_tickangle=-40,
        yaxis_title="Core-hours",
        height=560,
        margin=_M_BAR,
    )
    return _no_clip(fig)


_STATUS_COLORS = {"COMPLETED": "#2ecc71", "FAILED": "#e74c3c"}


def fig_status_per_step(trace_all: pd.DataFrame) -> go.Figure:
    """Stacked bar: COMPLETED vs FAILED counts per pipeline step."""
    agg = (
        trace_all.groupby(["pipeline_step", "status"], observed=True)
        .size()
        .reset_index(name="count")
    )
    # Totals for failure-rate annotation
    totals = agg.groupby("pipeline_step", observed=True)["count"].sum()
    failed = agg[agg["status"] == "FAILED"].set_index("pipeline_step")["count"]
    rate = (failed / totals * 100).fillna(0).round(1)

    fig = px.bar(
        agg, x="pipeline_step", y="count", color="status",
        color_discrete_map=_STATUS_COLORS,
        barmode="stack",
        labels={"count": "Jobs", "pipeline_step": "", "status": "Status"},
        title="Job completion status per pipeline step",
        height=520,
        text="count",
    )
    fig.update_traces(textposition="inside", textfont_size=11)
    # Annotate failure rate above each bar
    for step, r in rate.items():
        if r > 0:
            fig.add_annotation(
                x=step, y=totals[step],
                text=f"{r}% fail",
                showarrow=False,
                yanchor="bottom", yshift=4,
                font=dict(size=10, color="#c0392b"),
            )
    fig.update_layout(showlegend=True, xaxis_tickangle=-40, margin=_M_BAR)
    return _no_clip(fig)


def fig_status_overall(trace_all: pd.DataFrame) -> go.Figure:
    """Pie chart of overall COMPLETED vs FAILED across the whole pipeline."""
    agg = trace_all["status"].value_counts().reset_index()
    agg.columns = ["status", "count"]
    colors = [_STATUS_COLORS.get(s, "#888") for s in agg["status"]]
    fig = go.Figure(go.Pie(
        labels=agg["status"],
        values=agg["count"],
        marker_colors=colors,
        hole=0.35,
        textinfo="label+percent+value",
    ))
    fig.update_layout(
        title="Overall job status (whole pipeline)",
        height=520,
        margin=dict(t=70, b=30, l=30, r=150),
    )
    return fig


def fig_scaling_size_weighted(
    hjh: pd.DataFrame,
    extrap_df: pd.DataFrame,
    n_current: int,
    n_target: int,
) -> go.Figure:
    """Grouped bar: count-based linear vs size-weighted projection per step."""
    scale = n_target / n_current
    linear = hjh.groupby("pipeline_step", observed=True)["billed_core_hours"].sum().reset_index()
    linear["linear_projected"] = linear["billed_core_hours"] * scale
    merged = linear.merge(
        extrap_df[["pipeline_step", "total_predicted", "pct_extrapolated"]],
        on="pipeline_step", how="outer",
    ).sort_values("linear_projected", ascending=False, na_position="last")

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=merged["pipeline_step"], y=merged["linear_projected"],
        name="Count-based linear", marker_color="#3498db",
        text=merged["linear_projected"].apply(
            lambda x: f"{x / 1000:.1f}k" if pd.notna(x) else ""),
        textposition="outside",
    ))
    fig.add_trace(go.Bar(
        x=merged["pipeline_step"], y=merged["total_predicted"],
        name="Size-weighted", marker_color="#e67e22",
        text=merged["total_predicted"].apply(
            lambda x: f"{x / 1000:.1f}k" if pd.notna(x) else ""),
        textposition="outside",
    ))
    for _, row in merged.iterrows():
        pct = row.get("pct_extrapolated", 0)
        if pd.notna(pct) and pct > 1:
            fig.add_annotation(
                x=row["pipeline_step"],
                y=max(row.get("linear_projected", 0) or 0,
                      row.get("total_predicted", 0) or 0),
                text=f"{pct:.0f}% extrap.",
                showarrow=False, yanchor="bottom", yshift=25,
                font=dict(size=9, color="#c0392b"),
            )
    fig.update_layout(
        barmode="group",
        title=f"Projected billed core-hours: count-based vs size-weighted ({n_target:,} samples)",
        xaxis_tickangle=-40, yaxis_title="Projected billed core-hours",
        height=560, margin=_M_BAR,
    )
    return _no_clip(fig)


def fig_model_diagnostics_table(
    models_df: pd.DataFrame, extrap_df: pd.DataFrame,
) -> go.Figure:
    """Table with per-step binned-regression diagnostics."""
    merged = models_df.merge(
        extrap_df[["pipeline_step", "pct_extrapolated"]],
        on="pipeline_step", how="left",
    )
    fig = go.Figure(go.Table(
        header=dict(
            values=["Step", "Intercept (h)", "Slope (h/GB)", "R\u00b2",
                    "N obs", "N bins", "Size range (GB)", "% cost from extrapolation"],
            fill_color="steelblue", font_color="white", align="left",
        ),
        cells=dict(values=[
            merged["pipeline_step"],
            merged["intercept"].round(4),
            merged["slope"].round(4),
            merged["r_squared"].round(3),
            merged["n_obs"],
            merged["n_bins"],
            merged.apply(
                lambda r: f"{r['size_min_GB']:.2f} \u2013 {r['size_max_GB']:.2f}", axis=1),
            merged["pct_extrapolated"].round(1),
        ], align="left"),
    ))
    fig.update_layout(title="Per-step binned-regression diagnostics", height=450)
    return fig


def fig_size_vs_cost_scatter(
    hjh_with_size: pd.DataFrame, models_df: pd.DataFrame,
    bins_df: pd.DataFrame | None = None,
) -> go.Figure:
    """Scatter: raw data (faded) + bin averages (bold) + regression lines."""
    plot_df = hjh_with_size.dropna(subset=["size_GB", "billed_core_hours"])
    # Raw data points — very faded
    fig = px.scatter(
        plot_df, x="size_GB", y="billed_core_hours",
        color="pipeline_step", color_discrete_map=STEP_COLOR,
        hover_data=["sra_id"],
        labels={"size_GB": "Input file size (GB)",
                "billed_core_hours": "Billed core-hours"},
        title="Per-job cost vs input file size (binned regression)",
        opacity=0.15, height=560,
    )
    # Bin averages — bold markers
    if bins_df is not None and len(bins_df) > 0:
        for step in bins_df["pipeline_step"].unique():
            sb = bins_df[bins_df["pipeline_step"] == step]
            fig.add_trace(go.Scatter(
                x=sb["mean_size_GB"], y=sb["mean_cost"],
                mode="markers",
                marker=dict(size=10, color=STEP_COLOR.get(step, "#888"),
                            line=dict(width=1, color="white")),
                name=f"{step} (bin avg, n={sb['count'].sum()})",
                showlegend=True,
            ))
    # Regression lines
    for _, m in models_df.iterrows():
        x_range = np.linspace(m["size_min_GB"], m["size_max_GB"], 50)
        y_range = m["intercept"] + m["slope"] * x_range
        fig.add_trace(go.Scatter(
            x=x_range, y=y_range, mode="lines",
            line=dict(dash="dash", width=2,
                      color=STEP_COLOR.get(m["pipeline_step"], "#888")),
            name=f"{m['pipeline_step']} fit (R\u00b2={m['r_squared']:.2f})",
            showlegend=True,
        ))
    fig.update_layout(margin=_M_SCAT)
    return fig


def fig_size_distribution(
    hjh_with_size: pd.DataFrame, full_sizes: pd.DataFrame,
) -> go.Figure:
    """Overlaid histograms: batch vs full dataset file-size distributions."""
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=full_sizes["size_GB"].dropna(), nbinsx=80,
        name=f"Full dataset ({len(full_sizes):,})",
        opacity=0.6, histnorm="probability density",
    ))
    batch_sizes = hjh_with_size.drop_duplicates("sra_id")["size_GB"].dropna()
    fig.add_trace(go.Histogram(
        x=batch_sizes, nbinsx=80,
        name=f"Batch ({len(batch_sizes):,})",
        opacity=0.6, histnorm="probability density",
    ))
    fig.update_layout(
        barmode="overlay",
        title="Input file size distribution: batch vs full target dataset",
        xaxis_title="File size (GB)", yaxis_title="Density",
        height=460, margin=_M_HIST,
    )
    return fig


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

SECTION_STYLE = "margin: 40px 0 10px 0; border-bottom: 2px solid #ddd; padding-bottom: 6px;"
CARD_STYLE = "background:#f9f9f9; border-radius:8px; padding:18px 24px; margin:10px 0; display:inline-block; min-width:180px;"


def _kpi(label: str, value: str, sub: str = "") -> str:
    return (
        f'<div style="{CARD_STYLE}">'
        f'<div style="font-size:2em;font-weight:700;color:#1a6bad">{value}</div>'
        f'<div style="font-size:0.95em;color:#555">{label}</div>'
        f'{"<div style=&quot;font-size:0.8em;color:#888&quot;>" + sub + "</div>" if sub else ""}'
        f"</div>"
    )


def build_html(figures: dict, kpis: list[tuple]) -> str:
    # Embed plotly JS once inline with the first figure so it is always
    # available before any Plotly.newPlot() call executes.
    _js_included = [False]

    def embed(fig, cls="plot-full"):
        first = not _js_included[0]
        _js_included[0] = True
        html = fig.to_html(
            full_html=False,
            include_plotlyjs=first,          # True for first fig, False after
            config={"responsive": True},
        )
        return f'<div class="{cls}">{html}</div>'

    def section(title):
        return f'<h2 style="{SECTION_STYLE}">{title}</h2>'

    parts = ["""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>gutMAP Pipeline Report</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
         max-width: 1400px; margin: auto; padding: 20px 30px; background: #fff; }
  h1 { color: #1a3a5c; }
  h2 { color: #1a6bad; }
  .kpi-row { display: flex; flex-wrap: wrap; gap: 16px; margin: 20px 0; }
  .plot-row { display: grid; grid-template-columns: 1fr 1fr; gap: 10px; }
  .plot-half { min-width: 0; overflow: hidden; }
  .plot-full { width: 100%; }
</style>
</head>
<body>
<h1>gutMAP Pipeline Performance Report</h1>
<p style="color:#666">Analysing 1,000-sample HPC run &mdash; basis for extrapolation to 100k samples.</p>
"""]

    # KPI row
    parts.append('<div class="kpi-row">')
    for label, value, sub in kpis:
        parts.append(_kpi(label, value, sub))
    parts.append("</div>")

    parts.append(section("1 &mdash; Overview"))
    parts.append(embed(figures["summary_table"]))

    parts.append(section("2 &mdash; Job completion status"))
    parts.append('<div class="plot-row">')
    parts.append(embed(figures["status_per_step"], "plot-half"))
    parts.append(embed(figures["status_overall"], "plot-half"))
    parts.append("</div>")

    parts.append(section("3 &mdash; Billed core-hours"))
    parts.append('<div class="plot-row">')
    parts.append(embed(figures["billed_bar"], "plot-half"))
    parts.append(embed(figures["billed_pie"], "plot-half"))
    parts.append("</div>")
    parts.append(embed(figures["billed_vs_cpu"]))

    parts.append(section("4 &mdash; Core-hours (CPU used)"))
    parts.append('<div class="plot-row">')
    parts.append(embed(figures["cpu_share_bar"], "plot-half"))
    parts.append(embed(figures["cpu_share_pie"], "plot-half"))
    parts.append("</div>")

    parts.append(section("5 &mdash; CPU efficiency"))
    parts.append(embed(figures["efficiency_box"]))
    parts.append(embed(figures["efficiency_weighted"]))

    parts.append(section("6 &mdash; Memory"))
    parts.append(embed(figures["memory_box"]))
    parts.append(embed(figures["mem_waste"]))

    parts.append(section("7 &mdash; I/O and RSS (Nextflow trace)"))
    parts.append('<div class="plot-row">')
    parts.append(embed(figures["io"], "plot-half"))
    parts.append(embed(figures["rss_box"], "plot-half"))
    parts.append("</div>")

    parts.append(section("8 &mdash; Concurrency"))
    parts.append(embed(figures["concurrency"]))

    parts.append(section("9 &mdash; Sample-level variability"))
    parts.append('<div class="plot-row">')
    parts.append(embed(figures["per_sample_cpu"], "plot-half"))
    parts.append(embed(figures["outlier_table"], "plot-half"))
    parts.append("</div>")

    parts.append(section("10 &mdash; 100k scaling projection (count-based)"))
    parts.append(embed(figures["scaling"]))

    if "scaling_size_weighted" in figures:
        parts.append(section("11 &mdash; Size-weighted scaling projection"))
        parts.append(
            '<p style="color:#555; margin:0 0 10px 0">'
            "Projects costs to the full target dataset by modeling per-step cost as a "
            "linear function of input file size (<code>cost = a + b &times; size_GB</code>), "
            "then applying each step&rsquo;s model to every sample in the full dataset. "
            "Costs predicted outside the observed size range are flagged as extrapolation.</p>"
        )
        parts.append(embed(figures["scaling_size_weighted"]))
        parts.append(embed(figures["model_diagnostics"]))
        parts.append('<div class="plot-row">')
        parts.append(embed(figures["size_vs_cost_scatter"], "plot-half"))
        parts.append(embed(figures["size_distribution"], "plot-half"))
        parts.append("</div>")

    parts.append("</body></html>")
    return "\n".join(parts)


