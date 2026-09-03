# SPDX-FileCopyrightText: 2024-present Proxima Fusion GmbH <info@proximafusion.com>
#
# SPDX-License-Identifier: MIT
"""Tables, plots and the markdown report from one or more result ledgers, and the
comparison of two ledgers of the same matrix (baseline against candidate).

Figures are static PNGs drawn with matplotlib on a light surface; colours follow
one categorical order (never cycled), one sequential hue, and a fixed status
palette for outcome groups.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .cases import DEFAULT_MANIFEST, load_case, load_manifest
from .classify import CLASSES, DESCRIPTIONS
from .settings import FIELDS, LEVELS, REGIMES

# ---------------------------------------------------------------------------
# Palette (validated reference instance) and style
# ---------------------------------------------------------------------------

SERIES = [
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#4a3aa7",
    "#e34948",
]
SEQUENTIAL = [
    "#cde2fb",
    "#9ec5f4",
    "#6da7ec",
    "#3987e5",
    "#256abf",
    "#184f95",
    "#0d366b",
]
STATUS = {
    "good": "#0ca30c",
    "warning": "#fab219",
    "serious": "#ec835a",
    "critical": "#d03b3b",
}
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"
DEEMPHASIS = "#c3c2b7"

GROUP_OF_CLASS = {
    "converged": "converged",
    "max_iter_near": "still moving",
    "max_iter_slow": "still moving",
    "max_iter_stalled": "stuck",
    "max_iter_oscillating": "stuck",
    "vacuum_grid_exceeded": "aborted",
    "jacobian_75_times": "aborted",
    "bad_jacobian": "aborted",
    "early_failure": "aborted",
    "nan_residual": "aborted",
    "timeout": "aborted",
    "error_other": "aborted",
}
GROUPS = ["converged", "still moving", "stuck", "aborted"]
GROUP_COLOUR = {
    "converged": STATUS["good"],
    "still moving": STATUS["warning"],
    "stuck": STATUS["serious"],
    "aborted": STATUS["critical"],
}
GROUP_MARKER = {"converged": "o", "still moving": "^", "stuck": "s", "aborted": "x"}
REGIME_ORDER = ["vacuum", "vacuum_fixed", "beta_scaled", "beta1", "beta2_current"]


def _style() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.edgecolor": AXIS,
            "axes.labelcolor": INK_SECONDARY,
            "axes.titlecolor": INK,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.frameon": False,
            "lines.linewidth": 1.6,
            "font.family": "sans-serif",
            "figure.dpi": 110,
        }
    )


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_ledgers(results_dirs: list[Path]) -> pd.DataFrame:
    records: list[dict] = []
    for d in results_dirs:
        path = Path(d) / "results.jsonl"
        if not path.exists():
            continue
        with open(path) as f:
            for line in f:
                if line.strip():
                    rec = json.loads(line)
                    rec["results_dir"] = str(d)
                    records.append(rec)
    if not records:
        return pd.DataFrame()
    lcfs = [rec.pop("lcfs", None) for rec in records]
    df = pd.json_normalize(records, sep=".")
    df["lcfs"] = lcfs
    # The same experiment can appear in several ledgers (a rerun directory
    # overlapping the ladder); the first listed directory wins.
    df = df.drop_duplicates(subset="key", keep="first").reset_index(drop=True)
    df["group"] = df["status"].map(GROUP_OF_CLASS).fillna("aborted")
    if "log.vacuum_on_iters" in df:
        df["vacuum_on_iter"] = df["log.vacuum_on_iters"].map(
            lambda v: v[0] if isinstance(v, list) and v else np.nan
        )
    if "restart_counts.BAD_JACOBIAN" not in df:
        df["restart_counts.BAD_JACOBIAN"] = 0
    if "restart_counts.BAD_PROGRESS" not in df:
        df["restart_counts.BAD_PROGRESS"] = 0
    for c in ("restart_counts.BAD_JACOBIAN", "restart_counts.BAD_PROGRESS"):
        df[c] = df[c].fillna(0).astype(int)
    return df


def load_meta(results_dirs: list[Path]) -> list[dict]:
    metas = []
    for d in results_dirs:
        p = Path(d) / "meta.json"
        if p.exists():
            metas.append(json.loads(p.read_text()))
    return metas


def with_manifest(df: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    m = manifest.rename(columns={"ID": "case_id"})
    # Index values that share a name with a solver output get a quasr_ prefix.
    m = m.rename(
        columns={
            c: f"quasr_{c}" for c in m.columns if c in df.columns and c != "case_id"
        }
    )
    out = df.merge(m, on="case_id", how="left")
    edges = [-np.inf, 1.2, 1.9, np.inf]
    out["dist_bin"] = pd.cut(
        out["coil_dist_over_a"],
        edges,
        labels=["near (<1.2 a)", "mid (1.2-1.9 a)", "far (>1.9 a)"],
    ).astype(str)
    out["iota_bin"] = pd.cut(
        out["mean_iota"].abs(),
        [-np.inf, 0.25, 0.75, np.inf],
        labels=["iota <= 0.25", "0.25-0.75", "> 0.75"],
    ).astype(str)
    return out


# ---------------------------------------------------------------------------
# Markdown helpers
# ---------------------------------------------------------------------------


def md_table(df: pd.DataFrame, floatfmt: str = "{:.3g}", index: bool = True) -> str:
    frame = df.reset_index() if index else df
    cols = list(frame.columns)

    def fmt(v) -> str:
        if isinstance(v, float):
            if np.isnan(v):
                return ""
            return floatfmt.format(v)
        return str(v)

    lines = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(fmt(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def pct(n: float, d: float) -> str:
    return f"{100.0 * n / d:.0f}%" if d else ""


def _percentiles(s: pd.Series) -> dict[str, float]:
    s = pd.to_numeric(s, errors="coerce").dropna()
    if s.empty:
        return {"n": 0, "median": np.nan, "p90": np.nan, "max": np.nan}
    return {
        "n": len(s),
        "median": float(s.median()),
        "p90": float(s.quantile(0.9)),
        "max": float(s.max()),
    }


def ecdf(ax, values, colour, label, lw=1.6, ls="-"):
    v = np.sort(
        np.asarray(pd.to_numeric(values, errors="coerce").dropna(), dtype=float)
    )
    if len(v) == 0:
        return
    y = np.arange(1, len(v) + 1) / len(v)
    ax.step(
        v, y, where="post", color=colour, label=f"{label} (n={len(v)})", lw=lw, ls=ls
    )


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def fig_outcomes(df: pd.DataFrame, out: Path) -> str | None:
    base = df[df["level"] == "base"] if "base" in set(df["level"]) else df
    regimes = [r for r in REGIME_ORDER if r in set(base["regime"])]
    if not regimes:
        return None
    fig, ax = plt.subplots(figsize=(7.5, 0.6 * len(regimes) + 1.8))
    for i, regime in enumerate(regimes):
        sub = base[base["regime"] == regime]
        n = len(sub)
        left = 0.0
        for g in GROUPS:
            k = int((sub["group"] == g).sum())
            if k == 0:
                continue
            w = 100.0 * k / n
            ax.barh(
                i,
                w,
                left=left,
                color=GROUP_COLOUR[g],
                edgecolor=SURFACE,
                linewidth=2,
                height=0.55,
            )
            if w >= 7:
                ax.text(
                    left + w / 2,
                    i,
                    f"{k}",
                    ha="center",
                    va="center",
                    color=SURFACE if g != "still moving" else INK,
                    fontsize=8,
                )
            left += w
    ax.set_yticks(range(len(regimes)))
    ax.set_yticklabels(regimes)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of cases [%]")
    ax.grid(axis="y", visible=False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=GROUP_COLOUR[g]) for g in GROUPS]
    ax.set_title("Outcome of the free-boundary solve per regime (base level)")
    fig.legend(handles, GROUPS, ncol=4, loc="lower center", frameon=False)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    path = out / "fig_outcomes.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path.name


def fig_outcome_map(df: pd.DataFrame, out: Path) -> str | None:
    base = df[df["level"] == "base"] if "base" in set(df["level"]) else df
    regimes = [r for r in REGIME_ORDER if r in set(base["regime"])]
    if not regimes or "coil_dist_over_a" not in base:
        return None
    n = len(regimes)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 3.6), sharey=True, squeeze=False)
    for ax, regime in zip(axes[0], regimes, strict=True):
        sub = base[base["regime"] == regime]
        for g in GROUPS:
            s = sub[sub["group"] == g]
            if s.empty:
                continue
            style = {"edgecolor": SURFACE, "linewidth": 0.6} if g != "aborted" else {}
            ax.scatter(
                s["coil_dist_over_a"],
                s["aspect_ratio"],
                s=18,
                marker=GROUP_MARKER[g],
                color=GROUP_COLOUR[g],
                label=f"{g} ({len(s)})",
                alpha=0.9,
                **style,
            )
        ax.set_xscale("log")
        ax.set_xlabel("coil-surface distance / a")
        ax.set_title(regime)
        ax.legend(loc="upper left", markerscale=1.1)
    axes[0][0].set_ylabel("aspect ratio")
    fig.suptitle(
        "Where the solve succeeds: each dot is one QUASR case",
        y=1.02,
        fontsize=10,
        color=INK,
    )
    fig.tight_layout()
    path = out / "fig_outcome_map.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path.name


def fig_residual_gallery(
    df: pd.DataFrame, out: Path, regime: str = "vacuum"
) -> str | None:
    sub = (
        df[(df["regime"] == regime) & (df["level"] == "base")]
        if "fsq_trace" in df
        else pd.DataFrame()
    )
    if sub.empty:
        return None
    classes = [
        c
        for c in CLASSES
        if c in set(sub["status"]) and c not in ("timeout", "error_other")
    ]
    if not classes:
        return None
    n = len(classes)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(3.8 * ncol, 2.8 * nrow), squeeze=False, sharey=True
    )
    for k, cls in enumerate(classes):
        ax = axes[k // ncol][k % ncol]
        rows = sub[sub["status"] == cls]
        sample = rows.sample(min(12, len(rows)), random_state=0)
        for j, (_, r) in enumerate(sample.iterrows()):
            trace = r["fsq_trace"]
            if not isinstance(trace, list):
                continue
            y = np.array([np.nan if v is None else v for v in trace], dtype=float)
            x = np.arange(len(y)) * float(r.get("fsq_trace_step", 1) or 1)
            ax.plot(
                x,
                y,
                color=DEEMPHASIS if j else SERIES[0],
                lw=0.9 if j else 1.6,
                zorder=2 if j else 3,
            )
        ax.axhline(np.log10(float(rows["ftol"].iloc[0])), color=MUTED, lw=0.8, ls=":")
        if not any(isinstance(t, list) and len(t) > 1 for t in sample["fsq_trace"]):
            ax.text(
                0.5,
                0.5,
                "no iterations recorded",
                transform=ax.transAxes,
                ha="center",
                color=MUTED,
                fontsize=8,
            )
        ax.set_title(f"{cls} ({len(rows)})", fontsize=9)
        ax.set_xlabel("iteration (all stages)")
        if k % ncol == 0:
            ax.set_ylabel("log10 force residual")
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle(
        f"Force-residual histories by outcome class ({regime}, up to 12 cases each; one highlighted)",
        fontsize=10,
        color=INK,
    )
    fig.tight_layout()
    path = out / f"fig_residuals_{regime}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path.name


def fig_accuracy(df: pd.DataFrame, out: Path) -> str | None:
    vac = df[(df["regime"] == "vacuum") & (df["level"] == "base") & df["converged"]]
    if vac.empty or "bn_rms" not in vac:
        return None
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    ax = axes[0]
    ecdf(ax, vac["bn_rms"], SERIES[0], "B.n/|B| RMS on VMEC++ boundary")
    ecdf(
        ax,
        vac["mgrid_floor_rms"],
        MUTED,
        "mgrid table error on boundary (floor)",
        ls="--",
    )
    ax.set_xscale("log")
    ax.set_xlabel("relative error")
    ax.set_ylabel("fraction of converged vacuum cases")
    ax.legend(loc="lower right")
    ax.set_title("Normal-field error against the coil field")

    ax = axes[1]
    ok = vac.dropna(subset=["fl_dev_rms", "bn_rms"])
    if not ok.empty:
        ax.scatter(
            ok["bn_rms"],
            ok["fl_dev_rms"],
            s=16,
            color=SERIES[0],
            alpha=0.8,
            edgecolor=SURFACE,
            linewidth=0.5,
        )
        x = np.logspace(np.log10(ok["bn_rms"].min()), np.log10(ok["bn_rms"].max()), 10)
        scale = float(np.median(2.0 * np.pi * ok["Rmajor_p"] / ok["Aminor_p"]))
        ax.plot(
            x,
            scale * x,
            color=MUTED,
            lw=0.9,
            ls=":",
            label="B.n/|B| x 2 pi R/a (median R/a)",
        )
        ax.legend(loc="upper left")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("B.n/|B| RMS")
    ax.set_ylabel("field-line deviation RMS [a]")
    ax.set_title("Traced field lines stray by the integrated B.n")

    ax = axes[2]
    ok = vac.dropna(subset=["fl_iota", "iota_edge"])
    if not ok.empty:
        v = ok["iota_edge"].abs()
        t = ok["fl_iota"].abs()
        ax.scatter(
            v,
            t,
            s=16,
            color=SERIES[0],
            alpha=0.8,
            edgecolor=SURFACE,
            linewidth=0.5,
            label="traced vs VMEC++ edge iota",
        )
        lim = [min(v.min(), t.min()) * 0.9, max(v.max(), t.max()) * 1.1]
        ax.plot(lim, lim, color=MUTED, lw=0.8, ls=":")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.set_xlabel("|iota| at the edge, VMEC++")
    ax.set_ylabel("|iota| from traced field lines")
    ax.set_title("Rotational transform cross-check")
    fig.tight_layout()
    path = out / "fig_accuracy.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path.name


def _ladder_cases(df: pd.DataFrame) -> pd.DataFrame:
    """Vacuum runs of the ladder subset only (reruns of other cases at ladder
    levels, such as the initial-Jacobian rerun, must not enter the ladders)."""
    sub = df[df["regime"] == "vacuum"]
    if "ladder" in sub:
        sub = sub[sub["ladder"].fillna(False).astype(bool)]
    return sub


def _ladder_panel(ax, df, levels, xvals, metric, xlabel, title, floor=True, ylog=True):
    sub = _ladder_cases(df)
    sub = sub[sub["level"].isin(levels) & (sub["field"] == "f1280")]
    if sub.empty or metric not in sub:
        ax.axis("off")
        return False
    pos = dict(zip(levels, xvals, strict=True))
    for _, rows in sub.groupby("case_id"):
        g = rows.set_index("level")
        xs, ys = [], []
        for lv in levels:
            if (
                lv in g.index
                and bool(g.loc[lv, "converged"])
                and np.isfinite(g.loc[lv, metric])
            ):
                xs.append(pos[lv])
                ys.append(float(g.loc[lv, metric]))
        if xs:
            ax.plot(xs, ys, color=DEEMPHASIS, lw=0.9, marker="o", ms=3, zorder=2)
        failed = [
            pos[lv]
            for lv in levels
            if lv in g.index and not bool(g.loc[lv, "converged"])
        ]
        if failed:
            ax.scatter(failed, [np.nan] * len(failed))
    med_x, med_y, floor_y = [], [], []
    for lv in levels:
        g = sub[(sub["level"] == lv) & sub["converged"]]
        if len(g):
            med_x.append(pos[lv])
            med_y.append(float(pd.to_numeric(g[metric], errors="coerce").median()))
            floor_y.append(
                float(pd.to_numeric(g["mgrid_floor_rms"], errors="coerce").median())
                if "mgrid_floor_rms" in g
                else np.nan
            )
    if med_x:
        ax.plot(
            med_x,
            med_y,
            color=SERIES[0],
            lw=2.0,
            marker="o",
            ms=4,
            zorder=4,
            label="median of converged cases",
        )
        if floor and metric.startswith("bn"):
            ax.plot(
                med_x,
                floor_y,
                color=MUTED,
                lw=1.0,
                ls="--",
                zorder=3,
                label="median mgrid floor",
            )
    counts = [
        f"{int(((sub['level'] == lv) & sub['converged']).sum())}/{int((sub['level'] == lv).sum())}"
        for lv in levels
    ]
    ax.set_xticks(list(xvals))
    ax.set_xticklabels(
        [f"{lv}\n{c}" for lv, c in zip(levels, counts, strict=True)], fontsize=7
    )
    if ylog:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=9)
    ax.legend(loc="best")
    return True


def fig_ladders(df: pd.DataFrame, out: Path) -> str | None:
    have = set(df["level"])
    fourier = [lv for lv in ("m6", "m8", "m10", "m12") if lv in have]
    radial = [lv for lv in ("ns31", "m8", "ns71", "ns101") if lv in have]
    ftol = [lv for lv in ("ftol8", "m8", "ftol11", "ftol13") if lv in have]
    if len(fourier) < 2 and len(radial) < 2:
        return None
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7))
    _ladder_panel(
        axes[0][0],
        df,
        fourier,
        [LEVELS[lv].mpol for lv in fourier],
        "bn_rms",
        "mpol = ntor (ns = 51)",
        "B.n/|B| RMS: Fourier ladder",
    )
    _ladder_panel(
        axes[0][1],
        df,
        radial,
        [LEVELS[lv].ns_final for lv in radial],
        "bn_rms",
        "final ns (mpol = ntor = 8)",
        "B.n/|B| RMS: radial ladder",
    )
    _ladder_panel(
        axes[0][2],
        df,
        ftol,
        [-np.log10(LEVELS[lv].ftol) for lv in ftol],
        "bn_rms",
        "-log10 ftol (m8)",
        "B.n/|B| RMS: tolerance ladder",
    )
    _ladder_panel(
        axes[1][0],
        df,
        fourier,
        [LEVELS[lv].mpol for lv in fourier],
        "fl_dev_rms",
        "mpol = ntor (ns = 51)",
        "field-line deviation RMS [a]: Fourier ladder",
        floor=False,
    )
    _ladder_panel(
        axes[1][1],
        df,
        radial,
        [LEVELS[lv].ns_final for lv in radial],
        "fl_dev_rms",
        "final ns (mpol = ntor = 8)",
        "field-line deviation RMS [a]: radial ladder",
        floor=False,
    )
    _ladder_panel(
        axes[1][2],
        df,
        fourier,
        [LEVELS[lv].mpol for lv in fourier],
        "itfsq_total",
        "mpol = ntor (ns = 51)",
        "total iterations: Fourier ladder",
        floor=False,
    )
    for row in axes:
        for ax in row:
            ax.grid(True, which="major")
    fig.suptitle(
        "Resolution ladders on the ladder subset (grey: one case each; labels: converged/attempted)",
        fontsize=10,
        color=INK,
    )
    fig.tight_layout()
    path = out / "fig_ladders.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path.name


def fig_field_ladder(df: pd.DataFrame, out: Path) -> str | None:
    have = set(df["field"])
    fields = [f for f in ("f160", "f1280", "f2560") if f in have]
    sub = _ladder_cases(df)
    sub = sub[(sub["level"] == "m8") & sub["field"].isin(fields)]
    if len(fields) < 2 or sub.empty:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    pos = {f: i for i, f in enumerate(fields)}
    for metric, ax, title in (
        ("bn_rms", axes[0], "B.n/|B| RMS"),
        ("fl_dev_rms", axes[1], "field-line deviation RMS [a]"),
    ):
        for _, rows in sub.groupby("case_id"):
            g = rows[rows["converged"]].set_index("field")
            xs = [pos[f] for f in fields if f in g.index]
            ys = [float(g.loc[f, metric]) for f in fields if f in g.index]
            if xs:
                ax.plot(xs, ys, color=DEEMPHASIS, lw=0.9, marker="o", ms=3)
        med = [
            float(
                pd.to_numeric(
                    sub[(sub["field"] == f) & sub["converged"]][metric], errors="coerce"
                ).median()
            )
            for f in fields
        ]
        ax.plot(
            [pos[f] for f in fields],
            med,
            color=SERIES[0],
            lw=2.0,
            marker="o",
            ms=4,
            label="median",
        )
        if metric == "bn_rms":
            fl = [
                float(
                    pd.to_numeric(
                        sub[sub["field"] == f]["mgrid_floor_rms"], errors="coerce"
                    ).median()
                )
                for f in fields
            ]
            ax.plot(
                [pos[f] for f in fields],
                fl,
                color=MUTED,
                lw=1.0,
                ls="--",
                label="median mgrid floor",
            )
        ax.set_xticks(range(len(fields)))
        ax.set_xticklabels(
            [
                f"{f}\n{FIELDS[f].coil_points} pts, {FIELDS[f].points_rz}^2"
                for f in fields
            ],
            fontsize=8,
        )
        ax.set_yscale("log")
        ax.set_title(f"{title} vs field fidelity (m8)", fontsize=9)
        ax.legend(loc="best")
    fig.tight_layout()
    path = out / "fig_field_ladder.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path.name


def fig_speed(df: pd.DataFrame, out: Path) -> str | None:
    base = df[(df["level"] == "base") & df["converged"]]
    regimes = [r for r in REGIME_ORDER if r in set(base["regime"])]
    if not regimes:
        return None
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.6))
    for i, regime in enumerate(regimes):
        s = base[base["regime"] == regime]
        ecdf(axes[0], s["itfsq_total"], SERIES[i], regime)
        ecdf(axes[1], s["t_run"], SERIES[i], regime)
    axes[0].set_xlabel("total force-balance iterations (all stages)")
    axes[0].set_ylabel("fraction of converged cases")
    axes[0].set_title("Iterations to converge")
    axes[0].legend(loc="lower right")
    axes[1].set_xlabel("wall time of vmecpp.run [s]")
    axes[1].set_title("Wall time to converge")
    axes[1].set_xscale("log")
    axes[1].legend(loc="lower right")
    vac = base[base["regime"] == "vacuum"]
    if not vac.empty and "coil_dist_over_a" in vac:
        axes[2].scatter(
            vac["coil_dist_over_a"],
            vac["itfsq_total"],
            s=16,
            color=SERIES[0],
            alpha=0.8,
            edgecolor=SURFACE,
            linewidth=0.5,
        )
        axes[2].set_xscale("log")
        axes[2].set_xlabel("coil-surface distance / a")
        axes[2].set_ylabel("total iterations (vacuum)")
        axes[2].set_title("Iterations against coil proximity")
    fig.tight_layout()
    path = out / "fig_speed.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path.name


def fig_flow_control(df: pd.DataFrame, out: Path) -> str | None:
    have = set(df["level"])
    variants = [lv for lv in ("nvac1", "delt05") if lv in have]
    if not variants or "base" not in have:
        return None
    vac = df[df["regime"] == "vacuum"]
    base = vac[vac["level"] == "base"].set_index("case_id")
    fig, axes = plt.subplots(
        1, len(variants), figsize=(4.2 * len(variants), 3.8), squeeze=False
    )
    for ax, lv in zip(axes[0], variants, strict=True):
        v = vac[vac["level"] == lv].set_index("case_id")
        common = base.index.intersection(v.index)
        for cid in common:
            b, c = base.loc[cid], v.loc[cid]
            both = bool(b["converged"]) and bool(c["converged"])
            colour = SERIES[0] if both else STATUS["critical"]
            ax.scatter(
                b["itfsq_total"],
                c["itfsq_total"],
                s=18,
                color=colour,
                edgecolor=SURFACE,
                linewidth=0.5,
                alpha=0.9,
            )
        lim = [
            min(base["itfsq_total"].min(), v["itfsq_total"].min()) * 0.8,
            max(base["itfsq_total"].max(), v["itfsq_total"].max()) * 1.2,
        ]
        ax.plot(lim, lim, color=MUTED, lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("total iterations, base")
        ax.set_ylabel(f"total iterations, {lv}")
        nb = int(base.loc[common, "converged"].sum())
        nv = int(v.loc[common, "converged"].sum())
        ax.set_title(
            f"{lv}: {LEVELS[lv].description or 'flow-control variant'}\n"
            f"converged {nv}/{len(common)} (base {nb}/{len(common)}); "
            "red = not both converged",
            fontsize=8,
        )
    fig.tight_layout()
    path = out / "fig_flow_control.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path.name


def _cross_section_from_lcfs(lcfs: dict, phi: float, n: int = 400):
    xm = np.asarray(lcfs["xm"], dtype=float)
    xn = np.asarray(lcfs["xn"], dtype=float)
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ang = np.outer(theta, xm) - phi * xn[None, :]
    r = np.cos(ang) @ np.asarray(lcfs["rmnc"])
    z = np.sin(ang) @ np.asarray(lcfs["zmns"])
    return np.append(r, r[0]), np.append(z, z[0])


def fig_examples(
    df: pd.DataFrame, out: Path, cache_dir: Path | None = None
) -> str | None:
    vac = df[
        (df["regime"] == "vacuum") & (df["level"] == "base") & df["converged"]
    ].dropna(subset=["bn_rms"])
    if vac.empty or "lcfs" not in vac:
        return None
    order = vac.sort_values("bn_rms")
    picks = [
        ("best", order.iloc[0]),
        ("median", order.iloc[len(order) // 2]),
        ("worst", order.iloc[-1]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))
    for ax, (label, r) in zip(axes, picks, strict=True):
        lcfs = r["lcfs"]
        if not isinstance(lcfs, dict):
            ax.axis("off")
            continue
        try:
            case = load_case(int(r["case_id"]), cache_dir)
        except Exception:  # noqa: BLE001
            ax.axis("off")
            continue
        nfp = int(r["nfp"])
        for k, phi_fraction in enumerate((0.0, 0.25 / nfp)):
            sec = case.boundary.cross_section(phi_fraction, thetas=400)
            rq = np.hypot(sec[:, 0], sec[:, 1])
            zq = sec[:, 2]
            ax.plot(
                np.append(rq, rq[0]),
                np.append(zq, zq[0]),
                color=INK,
                lw=1.2,
                label="QUASR surface" if k == 0 else None,
            )
            rv, zv = _cross_section_from_lcfs(lcfs, 2 * np.pi * phi_fraction)
            ax.plot(
                rv,
                zv,
                color=SERIES[0],
                lw=1.2,
                ls="--",
                label="VMEC++ free boundary" if k == 0 else None,
            )
        ax.set_aspect("equal")
        dev = r.get("fl_dev_rms", np.nan)
        dev_text = f"{dev:.1e} a" if np.isfinite(dev) else "lines lost"
        ax.set_title(
            f"{label}: QUASR {int(r['case_id'])} (nfp={nfp})\n"
            f"B.n/|B| RMS = {r['bn_rms']:.1e}, field-line dev = {dev_text}",
            fontsize=8,
        )
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.legend(loc="upper right", fontsize=7)
    fig.suptitle(
        "Cross-sections at phi = 0 and a quarter period: converged vacuum solutions against the QUASR surface",
        fontsize=10,
        color=INK,
    )
    fig.tight_layout()
    path = out / "fig_examples.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path.name


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def table_outcomes(df: pd.DataFrame) -> str:
    base = df[df["level"] == "base"] if "base" in set(df["level"]) else df
    regimes = [r for r in REGIME_ORDER if r in set(base["regime"])]
    rows = []
    for cls in CLASSES:
        row = {"class": cls, "meaning": DESCRIPTIONS[cls]}
        any_nonzero = False
        for regime in regimes:
            sub = base[base["regime"] == regime]
            k = int((sub["status"] == cls).sum())
            any_nonzero |= k > 0
            row[regime] = f"{k} ({pct(k, len(sub))})" if len(sub) else ""
        if any_nonzero:
            rows.append(row)
    rows.append(
        {
            "class": "total",
            "meaning": "",
            **{r: str(int((base["regime"] == r).sum())) for r in regimes},
        }
    )
    return md_table(pd.DataFrame(rows), index=False)


def table_rate_by(df: pd.DataFrame, column: str, label: str) -> str:
    base = df[df["level"] == "base"] if "base" in set(df["level"]) else df
    regimes = [r for r in REGIME_ORDER if r in set(base["regime"])]
    if column not in base:
        return ""
    rows = []
    for value, sub in base.groupby(column, sort=True):
        row = {label: value}
        for regime in regimes:
            s = sub[sub["regime"] == regime]
            row[regime] = (
                f"{int(s['converged'].sum())}/{len(s)} ({pct(s['converged'].sum(), len(s))})"
                if len(s)
                else ""
            )
        rows.append(row)
    return md_table(pd.DataFrame(rows), index=False)


def table_accuracy(df: pd.DataFrame) -> str:
    vac = df[(df["regime"] == "vacuum") & (df["level"] == "base") & df["converged"]]
    metrics = [
        ("bn_rms", "B.n/|B| RMS"),
        ("bn_max", "B.n/|B| max"),
        ("fl_dev_rms", "field-line deviation RMS [a]"),
        ("fl_dev_max", "field-line deviation max [a]"),
        ("fl_lost_fraction", "fraction of traced lines lost"),
        ("b2_rms", "|B|^2 mismatch RMS"),
        ("quasr_dist_rms", "distance to QUASR surface RMS [a]"),
        ("quasr_dist_max", "distance to QUASR surface max [a]"),
        ("volume_ratio_quasr", "volume / QUASR volume"),
        ("delbsq_last", "delbsq (VMEC++ edge |B|^2 mismatch)"),
        ("mgrid_floor_rms", "mgrid table error on boundary RMS (floor)"),
        ("mgrid_floor_max", "mgrid table error on boundary max"),
        ("exact_field_check", "ground-truth self-check (max rel.)"),
    ]
    rows = []
    for col, name in metrics:
        if col in vac:
            rows.append({"metric": name, **_percentiles(vac[col])})
    if "bn_rms" in vac and "mgrid_floor_rms" in vac:
        ratio = vac["bn_rms"] / vac["mgrid_floor_rms"]
        rows.append({"metric": "B.n RMS / floor", **_percentiles(ratio)})
    if "fl_lost_quasr" in vac:
        rows.append(
            {
                "metric": "fraction of lines lost from the QUASR surface itself "
                "(computed for the cases with lost lines)",
                **_percentiles(vac["fl_lost_quasr"]),
            }
        )
    if "quasr_iota_edge" in vac and "iota_edge" in vac:
        rel = (vac["iota_edge"].abs() - vac["quasr_iota_edge"].abs()).abs() / vac[
            "quasr_iota_edge"
        ].abs()
        rows.append(
            {"metric": "edge iota: |VMEC++ - QUASR| / |QUASR|", **_percentiles(rel)}
        )
    if "fl_iota" in vac and "iota_edge" in vac:
        rel = (vac["fl_iota"].abs() - vac["iota_edge"].abs()).abs() / vac[
            "iota_edge"
        ].abs()
        rows.append(
            {"metric": "edge iota: |traced - VMEC++| / |VMEC++|", **_percentiles(rel)}
        )
    return md_table(pd.DataFrame(rows), index=False)


def confinement_summary(df: pd.DataFrame) -> str:
    """One sentence on where lost field lines come from."""
    vac = df[(df["regime"] == "vacuum") & (df["level"] == "base") & df["converged"]]
    if "fl_lost_fraction" not in vac:
        return ""
    lost = vac[
        (pd.to_numeric(vac["fl_lost_fraction"], errors="coerce").fillna(0) > 0)
        | ~np.isfinite(pd.to_numeric(vac.get("fl_dev_rms"), errors="coerce"))
    ]
    if lost.empty:
        return "No converged case lost a traced field line."
    half = lost[
        (pd.to_numeric(lost["fl_lost_fraction"], errors="coerce").fillna(1) >= 0.5)
    ]
    text = (
        f"{len(lost)} of the {len(vac)} converged vacuum cases lose at least one "
        f"of the traced field lines ({len(half)} lose half or more)."
    )
    if "fl_lost_quasr" in lost:
        q = pd.to_numeric(lost["fl_lost_quasr"], errors="coerce")
        known = q.notna().sum()
        text += (
            f" Of the {known} checked, the same tracing started on the QUASR surface "
            f"itself loses lines in {int((q > 0).sum())}: there the coil field's last "
            "closed flux surface lies inside the QUASR boundary and the lost lines "
            "say nothing about the solver."
        )
    return text


def table_ladder(df: pd.DataFrame) -> str:
    levels = [
        lv
        for lv in (
            "base",
            "m6",
            "m8",
            "m10",
            "m12",
            "ns31",
            "ns71",
            "ns101",
            "ftol8",
            "ftol11",
            "ftol13",
            "nvac1",
            "delt05",
        )
        if lv in set(df["level"])
    ]
    ladder_ids = (
        set(df.loc[df["level"] != "base", "case_id"])
        if "ladder" not in df
        else set(df.loc[df["ladder"].fillna(False).astype(bool), "case_id"])
    )
    vac = df[
        (df["regime"] == "vacuum")
        & df["case_id"].isin(ladder_ids)
        & (df["field"] == "f1280")
    ]
    rows = []
    for lv in levels:
        s = vac[vac["level"] == lv]
        if s.empty:
            continue
        c = s[s["converged"]]
        L = LEVELS[lv]
        rows.append(
            {
                "level": lv,
                "mpol": L.mpol,
                "ns": L.ns_final,
                "ftol": f"{L.ftol:.0e}",
                "converged": f"{len(c)}/{len(s)}",
                "B.n RMS median": float(c["bn_rms"].median())
                if "bn_rms" in c and len(c)
                else np.nan,
                "floor median": float(c["mgrid_floor_rms"].median())
                if "mgrid_floor_rms" in c and len(c)
                else np.nan,
                "fl dev RMS median [a]": float(c["fl_dev_rms"].median())
                if "fl_dev_rms" in c and len(c)
                else np.nan,
                "iterations median": float(c["itfsq_total"].median())
                if len(c)
                else np.nan,
                "wall time median [s]": float(c["t_run"].median())
                if len(c)
                else np.nan,
                "ms / iteration median": float(c["t_per_iter_ms"].median())
                if "t_per_iter_ms" in c and len(c)
                else np.nan,
            }
        )
    return md_table(pd.DataFrame(rows), index=False)


LADDERS = {
    "Fourier (ns = 51)": ("m6", "m8", "m10", "m12"),
    "radial (mpol = 8)": ("ns31", "m8", "ns71", "ns101"),
    "tolerance (m8)": ("ftol8", "m8", "ftol11", "ftol13"),
}


def table_ladder_paired(df: pd.DataFrame) -> str:
    """Per ladder, the cases that converged at every level of it, and the
    median of each metric over that common subset at each level. Unlike the
    per-level medians, these compare like with like."""
    vac = _ladder_cases(df)
    vac = vac[vac["field"] == "f1280"]
    rows = []
    for name, levels in LADDERS.items():
        if not all(lv in set(vac["level"]) for lv in levels):
            continue
        ok = None
        for lv in levels:
            ids = set(vac.loc[(vac["level"] == lv) & vac["converged"], "case_id"])
            ok = ids if ok is None else ok & ids
        ok = ok or set()
        for lv in levels:
            s_lv = vac[(vac["level"] == lv) & vac["case_id"].isin(ok)]
            rows.append(
                {
                    "ladder": name,
                    "level": lv,
                    "cases converged at every level": len(ok),
                    "B.n RMS median": float(s_lv["bn_rms"].median())
                    if len(s_lv)
                    else np.nan,
                    "B.n RMS max": float(s_lv["bn_rms"].max()) if len(s_lv) else np.nan,
                    "fl dev RMS median [a]": float(s_lv["fl_dev_rms"].median())
                    if len(s_lv)
                    else np.nan,
                    "iterations median": float(s_lv["itfsq_total"].median())
                    if len(s_lv)
                    else np.nan,
                }
            )
    return md_table(pd.DataFrame(rows), index=False) if rows else ""


def table_speed(df: pd.DataFrame) -> str:
    base = df[(df["level"] == "base") & df["converged"]]
    regimes = [r for r in REGIME_ORDER if r in set(base["regime"])]
    rows = []
    for regime in regimes:
        s = base[base["regime"] == regime]
        rows.append(
            {
                "regime": regime,
                "converged": len(s),
                "iterations median": float(s["itfsq_total"].median()),
                "iterations p90": float(s["itfsq_total"].quantile(0.9)),
                "wall time median [s]": float(s["t_run"].median()),
                "wall time p90 [s]": float(s["t_run"].quantile(0.9)),
                "ms / iteration median": float(s["t_per_iter_ms"].median())
                if "t_per_iter_ms" in s
                else np.nan,
                "vacuum switched on at iteration (median)": float(
                    s["vacuum_on_iter"].median()
                )
                if "vacuum_on_iter" in s
                else np.nan,
                "bad-Jacobian restarts (median)": float(
                    s["restart_counts.BAD_JACOBIAN"].median()
                ),
                "mgrid build median [s]": float(s["t_mgrid"].median()),
            }
        )
    return md_table(pd.DataFrame(rows), index=False)


def table_finite_beta(df: pd.DataFrame) -> str:
    base = df[df["level"] == "base"]
    vac = base[(base["regime"] == "vacuum") & base["converged"]].set_index("case_id")
    rows = []
    beta_regimes = [
        r
        for r in REGIME_ORDER
        if r in set(base["regime"])
        and r in REGIMES
        and (REGIMES[r].target_beta > 0 or REGIMES[r].shift_target is not None)
    ]
    for regime in beta_regimes:
        s = base[(base["regime"] == regime)]
        c = s[s["converged"]].set_index("case_id")
        common = c.index.intersection(vac.index)
        shift = (c.loc[common, "raxis_phi0"] - vac.loc[common, "raxis_phi0"]) / c.loc[
            common, "Aminor_p"
        ]
        beta_ratio = c["betatotal"] / c["target_beta"].replace(0, np.nan)
        row = {
            "regime": regime,
            "converged": f"{len(c)}/{len(s)}",
            "target beta median (all cases)": float(s["target_beta"].median())
            if "target_beta" in s
            else np.nan,
            "target beta median (converged)": float(c["target_beta"].median())
            if "target_beta" in c and len(c)
            else np.nan,
            "achieved / target beta (median)": float(beta_ratio.median())
            if len(c)
            else np.nan,
            "axis shift / a vs vacuum (median)": float(shift.median())
            if len(common)
            else np.nan,
            "outboard shift in": f"{int((shift > 0).sum())}/{len(common)}",
            "runs aborted by grid overrun": int(
                (s["status"] == "vacuum_grid_exceeded").sum()
            ),
            "runs aborted by Jacobian": int((s["status"] == "jacobian_75_times").sum()),
        }
        if regime == "beta2_current" and len(c):
            row["ctor / prescribed (median)"] = float(
                (c["ctor"] / c["curtor"]).median()
            )
        rows.append(row)
    return md_table(pd.DataFrame(rows), index=False) if rows else ""


def table_control(df: pd.DataFrame) -> str:
    """Free-boundary outcome against the fixed-boundary control, per case."""
    base = df[df["level"] == "base"]
    free = base[base["regime"] == "vacuum"].set_index("case_id")
    fixed = base[base["regime"] == "vacuum_fixed"].set_index("case_id")
    common = free.index.intersection(fixed.index)
    if len(common) == 0:
        return ""
    rows = []
    for cls in CLASSES:
        ids = [c for c in common if free.loc[c, "status"] == cls]
        if not ids:
            continue
        fx = fixed.loc[ids, "status"]
        rows.append(
            {
                "free-boundary class": cls,
                "cases": len(ids),
                "fixed-boundary converged": int((fx == "converged").sum()),
                "fixed-boundary bad_jacobian": int((fx == "bad_jacobian").sum()),
                "fixed-boundary other failure": int(
                    ((fx != "converged") & (fx != "bad_jacobian")).sum()
                ),
            }
        )
    return md_table(pd.DataFrame(rows), index=False)


def table_badjac_rerun(df: pd.DataFrame) -> str:
    """Outcome of the cases that failed the initial Jacobian at the base level
    when rerun at higher Fourier resolution or with a spectrally condensed
    initial boundary, free-boundary and (where run) fixed-boundary."""
    base = df[(df["level"] == "base") & (df["regime"] == "vacuum")]
    ids = set(base.loc[base["status"] == "bad_jacobian", "case_id"])
    levels = [lv for lv in ("m8", "m10", "base_cond") if lv in set(df["level"])]
    sub = df[
        df["case_id"].isin(ids) & df["level"].isin(levels) & (df["field"] == "f1280")
    ]
    if sub.empty:
        return ""
    rows = []
    for lv in levels:
        free = sub[(sub["level"] == lv) & (sub["regime"] == "vacuum")]
        fixed = sub[(sub["level"] == lv) & (sub["regime"] == "vacuum_fixed")]
        counts = free["status"].value_counts()
        row = {
            "level": lv,
            "what changes": LEVELS[lv].description or "",
            "cases rerun": len(free),
            "free-boundary converged": int(counts.get("converged", 0)),
            "still bad_jacobian": int(counts.get("bad_jacobian", 0)),
            "jacobian_75_times": int(counts.get("jacobian_75_times", 0)),
            "other": int(
                len(free)
                - counts.get("converged", 0)
                - counts.get("bad_jacobian", 0)
                - counts.get("jacobian_75_times", 0)
            ),
            "fixed-boundary converged": f"{int(fixed['converged'].sum())}/{len(fixed)}"
            if len(fixed)
            else "",
            "B.n RMS median of converged": float(
                free.loc[free["converged"], "bn_rms"].median()
            )
            if free["converged"].any()
            else np.nan,
        }
        if lv == "base_cond" and "condense.spectral_width_reduction" in free:
            row["spectral width reduction (median)"] = float(
                pd.to_numeric(
                    free["condense.spectral_width_reduction"], errors="coerce"
                ).median()
            )
        rows.append(row)
    return md_table(pd.DataFrame(rows), index=False)


def table_condensed_all(df: pd.DataFrame) -> str:
    """Base level against the condensed-boundary level on every case run at
    both (vacuum, free boundary): status transitions and accuracy ratio."""
    vac = df[(df["regime"] == "vacuum") & (df["field"] == "f1280")]
    base = vac[vac["level"] == "base"].set_index("case_id")
    cond = vac[vac["level"] == "base_cond"].set_index("case_id")
    common = base.index.intersection(cond.index)
    if len(common) < 10:
        return ""
    b = base.loc[common]
    c = cond.loc[common]
    lines = [
        f"Cases run at both levels: {len(common)}. Converged at base: {int(b['converged'].sum())}; "
        f"with the condensed boundary: {int(c['converged'].sum())}. "
        f"Newly converging: {int((~b['converged'] & c['converged']).sum())}; "
        f"newly failing: {int((b['converged'] & ~c['converged']).sum())}.",
        "",
    ]
    both = common[b["converged"].to_numpy() & c["converged"].to_numpy()]
    if len(both):
        ratio_bn = (c.loc[both, "bn_rms"] / b.loc[both, "bn_rms"]).astype(float)
        ratio_it = (c.loc[both, "itfsq_total"] / b.loc[both, "itfsq_total"]).astype(
            float
        )
        lines.append(
            f"On the {len(both)} cases converged at both: B.n RMS ratio condensed/base median "
            f"{ratio_bn.median():.2f} (p10 {ratio_bn.quantile(0.1):.2f}, p90 {ratio_bn.quantile(0.9):.2f}); "
            f"iterations ratio median {ratio_it.median():.2f}."
        )
        lines.append("")
    cross = pd.crosstab(b["status"], c["status"])
    lines += [
        "Status transitions (rows: base, columns: condensed):",
        "",
        md_table(cross),
        "",
    ]
    return "\n".join(lines)


def table_failures(
    df: pd.DataFrame,
    regime: str = "vacuum",
    per_class: int = 4,
    csv_path: Path | None = None,
) -> str:
    base = df[(df["level"] == "base") & (df["regime"] == regime) & ~df["converged"]]
    if base.empty:
        return "none"
    cols = {
        "case_id": "case",
        "nfp": "nfp",
        "coil_dist_over_a": "coil dist / a",
        "aspect_ratio": "aspect",
        "mean_iota": "iota",
        "max_elongation": "max elong.",
        "status": "class",
        "ns_reached": "ns reached",
        "itfsq_total": "iterations",
        "fsqr": "fsqr",
        "fsqz": "fsqz",
        "restart_counts.BAD_JACOBIAN": "Jacobian resets",
        "log.n_convergence_problem": "delt resets",
        "exception": "message",
    }
    have = [c for c in cols if c in base]
    full = base.sort_values(["status", "case_id"])[have].rename(columns=cols)
    if csv_path is not None:
        full.to_csv(csv_path, index=False)
    t = (
        base.sort_values(["status", "case_id"])
        .groupby("status", sort=False)
        .head(per_class)[have]
        .rename(columns=cols)
    )
    if "message" in t:
        t["message"] = t["message"].map(
            lambda v: textwrap.shorten(str(v), 80) if isinstance(v, str) else ""
        )
    return md_table(t, index=False)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def build_report(
    results_dirs: list[Path], manifest_path: Path | None, out_dir: Path
) -> Path:
    _style()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = load_ledgers(results_dirs)
    if df.empty:
        msg = "no results found"
        raise FileNotFoundError(msg)
    manifest = load_manifest(manifest_path or DEFAULT_MANIFEST)
    df = with_manifest(df, manifest)
    metas = load_meta(results_dirs)

    figures = {
        "outcomes": fig_outcomes(df, out_dir),
        "map": fig_outcome_map(df, out_dir),
        "residuals": fig_residual_gallery(df, out_dir),
        "accuracy": fig_accuracy(df, out_dir),
        "ladders": fig_ladders(df, out_dir),
        "field": fig_field_ladder(df, out_dir),
        "speed": fig_speed(df, out_dir),
        "flow": fig_flow_control(df, out_dir),
        "examples": fig_examples(df, out_dir),
    }

    base = df[df["level"] == "base"]
    n_cases = base["case_id"].nunique()
    meta0 = metas[0] if metas else {}
    lines: list[str] = []
    a = lines.append
    a("# VMEC++ free-boundary bench: results")
    a("")
    a(f"Generated from {', '.join(str(d) for d in results_dirs)}.")
    a("")
    a("## Setup")
    a("")
    a(
        f"* Solver: VMEC++ at `{meta0.get('repo_describe', '?')}` (Python package version {meta0.get('vmecpp_version')}), SIMSOPT {meta0.get('simsopt_version')}, Python {meta0.get('python')}, {meta0.get('platform')}."
    )
    a(
        f"* Cases: {n_cases} QUASR configurations from the manifest ({int((manifest['subset'] == 'reference').sum())} reference, {int((manifest['subset'] == 'main').sum())} stratified), {int(manifest['ladder'].sum())} on the ladders."
    )
    a(
        f"* Runs in the ledger: {len(df)} across regimes {sorted(set(df['regime']))}, levels {sorted(set(df['level']))}, field specs {sorted(set(df['field']))}."
    )
    for m in metas:
        started = min([*m.get("previous_runs", []), m.get("started")], key=str)
        a(
            f"* `{m.get('run_name')}`: started {started}, last finished {m.get('finished', 'running')}, workers x threads = {m.get('workers')} x {m.get('threads')}, tracing {m.get('trace_lines')} lines x {m.get('trace_turns')} turns."
        )
    a("")
    a(
        "Regimes: "
        + "; ".join(
            f"`{k}`: {v.description}"
            for k, v in REGIMES.items()
            if k in set(df["regime"])
        )
        + "."
    )
    a("")
    a(
        "Levels: "
        + "; ".join(
            f"`{k}`: mpol=ntor={v.mpol}, ns={list(v.ns_array)}, nzeta={v.nzeta}, ftol={v.ftol:.0e}, niter={v.niter}, nvacskip={v.nvacskip}, delt={v.delt}"
            for k, v in LEVELS.items()
            if k in set(df["level"])
        )
        + "."
    )
    a("")
    a(
        "Field specs: "
        + "; ".join(
            f"`{k}`: {v.coil_points}-point coil polygons, {v.points_rz}^2 R-Z grid, margin {v.margin}"
            for k, v in FIELDS.items()
            if k in set(df["field"])
        )
        + "."
    )
    a("")
    a("## Robustness")
    a("")
    a("Outcome classes per regime at the base level (count and share of cases):")
    a("")
    a(table_outcomes(df))
    a("")
    if figures["outcomes"]:
        a(f"![outcomes]({figures['outcomes']})")
        a("")
    if figures["map"]:
        a(f"![outcome map]({figures['map']})")
        a("")
    a(
        "Convergence rate by coil-to-surface distance (in minor radii, terciles of the QUASR index):"
    )
    a("")
    a(table_rate_by(df, "dist_bin", "coil distance"))
    a("")
    a("By number of field periods:")
    a("")
    a(table_rate_by(df, "nfp", "nfp"))
    a("")
    a("By mean rotational transform:")
    a("")
    a(table_rate_by(df, "iota_bin", "iota"))
    a("")
    if figures["residuals"]:
        a(f"![residual histories]({figures['residuals']})")
        a("")
    control = table_control(df)
    if control:
        a(
            "Fixed-boundary control: the same case with the QUASR surface imposed. A free-boundary failure whose fixed-boundary twin also fails is a boundary-representation or initial-guess problem, not a free-boundary one."
        )
        a("")
        a(control)
        a("")
    rerun = table_badjac_rerun(df)
    if rerun:
        a(
            "Initial-Jacobian failures rerun at higher Fourier resolution (m8, m10: "
            "ns = 51) and with a spectrally condensed initial boundary (base_cond):"
        )
        a("")
        a(rerun)
        a("")
    condensed = table_condensed_all(df)
    if condensed:
        a("Condensed initial boundary on every case run at both levels:")
        a("")
        a(condensed)
        a("")
    a(
        "Non-converged vacuum runs at the base level, up to four per class "
        "(the full list is in `failures_vacuum.csv`):"
    )
    a("")
    a(table_failures(df, "vacuum", csv_path=out_dir / "failures_vacuum.csv"))
    a("")
    a("## Accuracy (vacuum, converged, base level)")
    a("")
    a(
        "All errors are relative to the exact Biot-Savart field of the coils, evaluated on the VMEC++ boundary. Distances are in minor radii. The floor is the error of the mgrid table itself on the boundary; a solver cannot beat it."
    )
    a("")
    a(table_accuracy(df))
    a("")
    a(confinement_summary(df))
    a("")
    if figures["accuracy"]:
        a(f"![accuracy]({figures['accuracy']})")
        a("")
    if figures["examples"]:
        a(f"![examples]({figures['examples']})")
        a("")
    if figures["ladders"] or figures["field"]:
        a("## Convergence with resolution (ladder subset, vacuum)")
        a("")
        a(table_ladder(df))
        a("")
        paired = table_ladder_paired(df)
        if paired:
            a(
                "Same ladders restricted to the cases that converged at every level (like against like):"
            )
            a("")
            a(paired)
            a("")
        if figures["ladders"]:
            a(f"![ladders]({figures['ladders']})")
            a("")
        if figures["field"]:
            a(f"![field ladder]({figures['field']})")
            a("")
    a("## Speed (converged runs, base level)")
    a("")
    a(
        "Wall times are laptop timings with three threads per run and three runs "
        "at a time; ledgers were produced under different background load, so "
        "iteration counts are the portable number and wall times are indicative "
        "within one ledger only."
    )
    a("")
    a(table_speed(df))
    a("")
    if figures["speed"]:
        a(f"![speed]({figures['speed']})")
        a("")
    if figures["flow"]:
        a(f"![flow control]({figures['flow']})")
        a("")
    fb = table_finite_beta(df)
    if fb:
        a("## Finite beta")
        a("")
        a(
            "Beta regimes have no exact ground truth (plasma currents contribute to the field), so only convergence and sanity are reported: the achieved beta, the outboard shift of the magnetic axis against the vacuum solution of the same case, and for the current-carrying regime the net toroidal current."
        )
        a("")
        a(fb)
        a("")
    a("## Ledger columns")
    a("")
    a(
        "Every run is one JSON record in `results.jsonl`; see `README.md` for the definition of the metrics and `classify.py` for the classes."
    )
    report_path = out_dir / "report.md"
    report_path.write_text("\n".join(lines))
    df.drop(columns=[c for c in ("fsq_trace", "lcfs") if c in df]).to_csv(
        out_dir / "ledger.csv", index=False
    )
    return report_path


# ---------------------------------------------------------------------------
# Comparison of two ledgers
# ---------------------------------------------------------------------------


def build_comparison(baseline: Path, candidate: Path, out_dir: Path) -> Path:
    _style()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    b = load_ledgers([baseline])
    c = load_ledgers([candidate])
    if b.empty or c.empty:
        msg = "both ledgers must contain results"
        raise FileNotFoundError(msg)
    j = b.merge(c, on="key", suffixes=("_b", "_c"))
    lines = [
        "# Bench comparison",
        "",
        f"Baseline `{baseline}` against candidate `{candidate}`: {len(j)} runs in common.",
        "",
    ]
    cross = pd.crosstab(j["status_b"], j["status_c"])
    lines += [
        "Status transitions (rows: baseline, columns: candidate):",
        "",
        md_table(cross),
        "",
    ]
    newly_failing = j[j["converged_b"] & ~j["converged_c"]]
    newly_converging = j[~j["converged_b"] & j["converged_c"]]
    lines += [f"Newly failing: {len(newly_failing)}", ""]
    if len(newly_failing):
        lines += [
            md_table(newly_failing[["key", "status_b", "status_c"]], index=False),
            "",
        ]
    lines += [f"Newly converging: {len(newly_converging)}", ""]
    if len(newly_converging):
        lines += [
            md_table(newly_converging[["key", "status_b", "status_c"]], index=False),
            "",
        ]
    both = j[j["converged_b"] & j["converged_c"]]
    if len(both):
        rows = []
        for col, name in (
            ("itfsq_total", "iterations"),
            ("t_run", "wall time"),
            ("bn_rms", "B.n RMS"),
            ("fl_dev_rms", "field-line deviation"),
        ):
            if f"{col}_b" in both and f"{col}_c" in both:
                ratio = pd.to_numeric(
                    both[f"{col}_c"], errors="coerce"
                ) / pd.to_numeric(both[f"{col}_b"], errors="coerce")
                rows.append(
                    {
                        "quantity": f"{name} candidate / baseline",
                        "n": int(ratio.notna().sum()),
                        "p10": float(ratio.quantile(0.1)),
                        "median": float(ratio.median()),
                        "p90": float(ratio.quantile(0.9)),
                    }
                )
        lines += [
            "Ratios on runs converged in both (below 1 is better for all four):",
            "",
            md_table(pd.DataFrame(rows), index=False),
            "",
        ]
        fig, ax = plt.subplots(figsize=(4.5, 4.2))
        ax.scatter(
            both["itfsq_total_b"],
            both["itfsq_total_c"],
            s=16,
            color=SERIES[0],
            alpha=0.8,
            edgecolor=SURFACE,
            linewidth=0.5,
        )
        lim = [
            min(both["itfsq_total_b"].min(), both["itfsq_total_c"].min()) * 0.8,
            max(both["itfsq_total_b"].max(), both["itfsq_total_c"].max()) * 1.2,
        ]
        ax.plot(lim, lim, color=MUTED, lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("iterations, baseline")
        ax.set_ylabel("iterations, candidate")
        ax.set_title("Iterations to converge, per run")
        fig.tight_layout()
        fig.savefig(
            out_dir / "fig_compare_iterations.png", dpi=150, bbox_inches="tight"
        )
        plt.close(fig)
        lines += ["![iterations](fig_compare_iterations.png)", ""]
    path = out_dir / "compare.md"
    path.write_text("\n".join(lines))
    return path
