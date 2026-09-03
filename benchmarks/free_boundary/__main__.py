# SPDX-FileCopyrightText: 2024-present Proxima Fusion GmbH <info@proximafusion.com>
#
# SPDX-License-Identifier: MIT
"""Command line of the free-boundary bench.

Run from the repository root::

    python -m benchmarks.free_boundary select --n-main 300
    python -m benchmarks.free_boundary fetch
    python -m benchmarks.free_boundary run --out results/main --levels base
    python -m benchmarks.free_boundary report --results results/main --out report
"""

from __future__ import annotations

import argparse
import concurrent.futures
import sys
from pathlib import Path

from .cases import (
    DEFAULT_MANIFEST,
    load_database,
    load_manifest,
    select_cases,
    serial_path,
)
from .report import build_comparison, build_report
from .runner import remeasure, retrace, run_matrix, worker_main


def _add_manifest_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)


def _selected_manifest(args):
    manifest = load_manifest(args.manifest)
    if getattr(args, "subset", "all") != "all":
        if args.subset == "ladder":
            manifest = manifest[manifest["ladder"].astype(bool)]
        else:
            manifest = manifest[manifest["subset"] == args.subset]
    if getattr(args, "ids", None):
        wanted = {int(i) for i in args.ids.split(",")}
        manifest = manifest[manifest["ID"].isin(wanted)]
    if getattr(args, "limit", None):
        manifest = manifest.head(args.limit)
    return manifest


def cmd_select(args) -> int:
    df = load_database(args.cache_dir)
    manifest = select_cases(
        df, n_main=args.n_main, n_ladder=args.n_ladder, seed=args.seed
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.out, index=False)
    print(f"wrote {len(manifest)} cases to {args.out}")  # noqa: T201
    print(manifest.groupby(["subset"]).size().to_string())  # noqa: T201
    return 0


def cmd_fetch(args) -> int:
    manifest = _selected_manifest(args)
    ids = [int(i) for i in manifest["ID"]]

    def fetch(config_id: int) -> str:
        return str(serial_path(config_id, args.cache_dir))

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        paths = list(pool.map(fetch, ids))
    print(f"{len(paths)} serial files available")  # noqa: T201
    return 0


def cmd_run(args) -> int:
    manifest = _selected_manifest(args)
    run_matrix(
        manifest,
        levels=args.levels.split(","),
        fields=args.fields.split(","),
        regimes=args.regimes.split(","),
        out_dir=args.out,
        workers=args.workers,
        threads=args.threads,
        timeout=args.timeout,
        trace=not args.no_trace,
        trace_lines=args.trace_lines,
        trace_turns=args.trace_turns,
        trace_tol=args.trace_tol,
        cache_dir=args.cache_dir,
        run_name=args.name or "",
        redo=tuple(args.redo.split(",")) if args.redo else (),
    )
    return 0


def cmd_retrace(args) -> int:
    n = retrace(
        Path(args.results),
        n_lines=args.trace_lines,
        n_turns=args.trace_turns,
        tol=args.trace_tol,
        only_missing=args.only_missing,
        only_confinement=args.only_confinement,
        ids={int(i) for i in args.ids.split(",")} if args.ids else None,
        cache_dir=args.cache_dir,
    )
    print(f"retraced {n} records")  # noqa: T201
    return 0


def cmd_remeasure(args) -> int:
    n = remeasure(Path(args.results), cache_dir=args.cache_dir)
    print(f"remeasured {n} records")  # noqa: T201
    return 0


def cmd_worker(args) -> int:
    return worker_main(args.spec, args.out)


def cmd_report(args) -> int:
    try:
        build_report(
            [Path(p) for p in args.results],
            args.manifest,
            args.out,
            ledger_csv=args.ledger_csv,
        )
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)  # noqa: T201
        return 2
    return 0


def cmd_compare(args) -> int:
    build_comparison(Path(args.baseline), Path(args.candidate), args.out)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m benchmarks.free_boundary")
    parser.add_argument("--cache-dir", type=Path, default=None)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("select", help="build the case manifest from the QUASR index")
    p.add_argument("--n-main", type=int, default=300)
    p.add_argument("--n-ladder", type=int, default=12)
    p.add_argument("--seed", type=int, default=20260902)
    p.add_argument("--out", type=Path, default=None)
    p.set_defaults(func=cmd_select)

    p = sub.add_parser("fetch", help="download the QUASR serial files of the manifest")
    _add_manifest_arg(p)
    p.add_argument(
        "--subset", default="all", choices=["all", "reference", "main", "ladder"]
    )
    p.add_argument("--ids", default=None)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--workers", type=int, default=8)
    p.set_defaults(func=cmd_fetch)

    p = sub.add_parser("run", help="run the bench matrix")
    _add_manifest_arg(p)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--name", default=None)
    p.add_argument("--levels", default="base")
    p.add_argument("--fields", default="f1280")
    p.add_argument("--regimes", default="vacuum,beta1,beta2_current")
    p.add_argument(
        "--subset", default="all", choices=["all", "reference", "main", "ladder"]
    )
    p.add_argument("--ids", default=None, help="comma-separated QUASR IDs")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--threads", type=int, default=3)
    p.add_argument("--timeout", type=float, default=1800.0)
    p.add_argument("--no-trace", action="store_true")
    p.add_argument("--trace-lines", type=int, default=8)
    p.add_argument("--trace-turns", type=int, default=100)
    p.add_argument("--trace-tol", type=float, default=1e-9)
    p.add_argument(
        "--redo", default=None, help="comma-separated statuses to drop and re-run"
    )
    p.set_defaults(func=cmd_run)

    p = sub.add_parser(
        "retrace", help="recompute field-line metrics from stored boundaries"
    )
    p.add_argument("--results", required=True)
    p.add_argument("--trace-lines", type=int, default=8)
    p.add_argument("--trace-turns", type=int, default=100)
    p.add_argument("--trace-tol", type=float, default=1e-9)
    p.add_argument("--only-missing", action="store_true")
    p.add_argument(
        "--only-confinement",
        action="store_true",
        help="only compute fl_lost_quasr where it is missing",
    )
    p.add_argument("--ids", default=None, help="comma-separated QUASR IDs")
    p.set_defaults(func=cmd_retrace)

    p = sub.add_parser(
        "remeasure", help="recompute normal-field metrics from stored boundaries"
    )
    p.add_argument("--results", required=True)
    p.set_defaults(func=cmd_remeasure)

    p = sub.add_parser("worker", help=argparse.SUPPRESS)
    p.add_argument("spec")
    p.add_argument("out")
    p.set_defaults(func=cmd_worker)

    p = sub.add_parser("report", help="tables and plots from one or more result dirs")
    _add_manifest_arg(p)
    p.add_argument("--results", nargs="+", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument(
        "--ledger-csv", action="store_true", help="also write a flat ledger.csv"
    )
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("compare", help="diff two result ledgers of the same matrix")
    p.add_argument("--baseline", required=True)
    p.add_argument("--candidate", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_compare)

    args = parser.parse_args(argv)
    if args.command == "select" and args.out is None:
        args.out = DEFAULT_MANIFEST
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
