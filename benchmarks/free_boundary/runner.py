# SPDX-FileCopyrightText: 2024-present Proxima Fusion GmbH <info@proximafusion.com>
#
# SPDX-License-Identifier: MIT
"""Run the bench matrix: cases x regimes x levels x field specs.

Each (case, level, field) triple is one worker subprocess that builds the
response table once and solves every requested regime with it. Subprocesses
isolate the bench from solver crashes, enforce a wall-clock limit per job and
let the solver's own log (the legacy VMEC table on stdout) be captured per
run. Results are appended to ``results.jsonl`` as they arrive, so a run can be
resumed and the ledger diffed between solver versions.
"""

from __future__ import annotations

import concurrent.futures
import contextlib
import ctypes
import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import threading
import time
import traceback
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import simsopt

import vmecpp

from . import cases as case_lib
from . import metrics
from .classify import classify, residual_tail_summary
from .settings import FIELDS, LEVELS, REGIMES, presets_as_dict

REPO_ROOT = case_lib.REPO_ROOT


@dataclass
class JobSpec:
    case_id: int
    level: str
    field: str
    regimes: list[str]
    threads: int = 2
    trace: bool = True
    trace_lines: int = 8
    trace_turns: int = 100
    trace_tol: float = 1e-9
    cache_dir: str = ""
    work_dir: str = ""
    run_name: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def key(self) -> str:
        return f"{self.case_id:07d}_{self.level}_{self.field}"


def result_key(case_id: int, regime: str, level: str, field_name: str) -> str:
    return f"{case_id:07d}:{regime}:{level}:{field_name}"


# ---------------------------------------------------------------------------
# Capturing the solver's stdout (C++ side) per run
# ---------------------------------------------------------------------------


def _openmp_runtimes() -> list[str]:
    """Paths of the OpenMP runtimes linked by the VMEC++ and SIMSOPT extensions.

    Each wheel may bundle its own copy, so every copy has to be told about the
    thread count separately.
    """
    found: list[str] = []
    for module_name in ("vmecpp.cpp._vmecpp", "simsoptpp"):
        try:
            module = __import__(module_name, fromlist=["__file__"])
            ext = Path(module.__file__)
        except Exception:  # noqa: BLE001, S112
            continue
        if sys.platform == "darwin":
            cmd = ["otool", "-L", str(ext)]
        else:
            cmd = ["ldd", str(ext)]
        try:
            out = subprocess.run(
                cmd, capture_output=True, text=True, check=False
            ).stdout
        except OSError:
            continue
        for ln in out.splitlines():
            parts = ln.split()
            if not parts or "omp" not in ln.lower():
                continue
            path = (
                parts[0]
                if sys.platform == "darwin"
                else (parts[2] if "=>" in ln else parts[0])
            )
            if path.startswith(("@loader_path", "@rpath")):
                candidate = ext.parent / path.split("/", 1)[1] if "/" in path else None
                if candidate is None or not candidate.exists():
                    # simsopt bundles its runtime next to the extension in .dylibs/
                    matches = list(ext.parent.glob(".dylibs/*omp*"))
                    candidate = matches[0] if matches else None
                path = str(candidate) if candidate else ""
            if path and path not in found:
                found.append(path)
    return found


_OPENMP_RUNTIMES: list[str] | None = None


def set_omp_threads(n: int) -> bool:
    """Best-effort omp_set_num_threads on every OpenMP runtime in the process.

    The SIMSOPT field-line tracer evaluates one point per call, so an OpenMP
    team of several threads makes it many times slower; the response-table
    builder on the other hand scales with threads. Returns False if no runtime
    could be reached, in which case OMP_NUM_THREADS stays in charge.
    """
    global _OPENMP_RUNTIMES  # noqa: PLW0603
    if _OPENMP_RUNTIMES is None:
        _OPENMP_RUNTIMES = _openmp_runtimes()
    ok = False
    for path in [*_OPENMP_RUNTIMES, "libomp.dylib", "libgomp.so.1", "libomp.so"]:
        try:
            ctypes.CDLL(path).omp_set_num_threads(int(n))
            ok = True
        except (OSError, AttributeError):
            continue
    return ok


def _flush_c_stdio() -> None:
    sys.stdout.flush()
    with contextlib.suppress(Exception):
        ctypes.CDLL(None).fflush(None)


class FdCapture:
    """Redirect file descriptor 1 into a temporary file for the duration of a
    block; the text is available as ``.text`` afterwards."""

    text: str = ""

    def __enter__(self):
        _flush_c_stdio()
        self._tmp = tempfile.TemporaryFile(mode="w+b")
        self._saved = os.dup(1)
        os.dup2(self._tmp.fileno(), 1)
        return self

    def __exit__(self, *_: object) -> None:
        _flush_c_stdio()
        os.dup2(self._saved, 1)
        os.close(self._saved)
        self._tmp.seek(0)
        self.text = self._tmp.read().decode(errors="replace")
        self._tmp.close()


def stage_iterations(text: str) -> dict[str, int]:
    """Last printed iteration of every multigrid stage, keyed by ns."""
    stages: dict[str, int] = {}
    current: str | None = None
    for ln in text.splitlines():
        header = re.match(r"\s*NS =\s+(\d+)", ln)
        if header:
            current = header.group(1)
            stages.setdefault(current, 0)
            continue
        row = re.match(r"\s*(\d+)\s*\|", ln)
        if row and current is not None:
            stages[current] = max(stages[current], int(row.group(1)))
    return stages


def parse_log(text: str) -> dict:
    """Markers of the solver's legacy log that matter for the taxonomy."""
    lines = text.splitlines()
    return {
        "stage_iters": stage_iterations(text),
        "vacuum_on_iters": [
            int(m)
            for m in re.findall(
                r"VACUUM PRESSURE TURNED ON AT\s+(\d+) ITERATIONS", text
            )
        ],
        "n_convergence_problem": text.count("HAVING A CONVERGENCE PROBLEM"),
        "n_initial_jacobian_sign": text.count("INITIAL JACOBIAN CHANGED SIGN"),
        "theta_flipped": "need to flip theta" in text,
        "grid_exceeded": "exceeded the vacuum field grid" in text,
        "fatal": [ln.strip()[:200] for ln in lines if "FATAL" in ln or "WARNING" in ln][
            :3
        ],
    }


# ---------------------------------------------------------------------------
# Worker: one (case, level, field), all regimes
# ---------------------------------------------------------------------------


def _first_line(text: str) -> str:
    text = text.strip()
    return text.splitlines()[0][:400] if text else ""


def run_worker(spec: JobSpec) -> list[dict]:
    level = LEVELS[spec.level]
    field_spec = FIELDS[spec.field]
    cache_dir = Path(spec.cache_dir) if spec.cache_dir else None
    work_dir = Path(spec.work_dir) if spec.work_dir else Path(tempfile.mkdtemp())

    t0 = time.time()
    case = case_lib.load_case(spec.case_id, cache_dir, spec.meta)
    t_load = time.time() - t0

    t0 = time.time()
    set_omp_threads(spec.threads)
    table = case_lib.build_response_table(case, field_spec, level.nzeta, work_dir)
    t_mgrid = time.time() - t0
    set_omp_threads(1)

    t0 = time.time()
    bs, check, ref_points = metrics.reference_field(case)
    case_info = {
        "phiedge": case.phiedge,
        "b_char": case.b_char,
        "r_char": case.r_char,
        "extcur": float(case.extcur[0]),
        "n_coils": len(case.coils),
        "n_base_coils": len(case.base_coils),
        "exact_field_check": check,
        "reference_coil_points": ref_points,
    }
    try:
        case_info.update(metrics.mgrid_floor(case, table, level.nzeta, bs))
    except Exception:  # noqa: BLE001
        case_info["case_metrics_error"] = _first_line(
            traceback.format_exc().splitlines()[-1]
        )
    t_case_metrics = time.time() - t0

    results: list[dict] = []
    for regime_name in spec.regimes:
        regime = REGIMES[regime_name]
        rec: dict = {
            "key": result_key(spec.case_id, regime_name, spec.level, spec.field),
            "case_id": spec.case_id,
            "regime": regime_name,
            "level": spec.level,
            "field": spec.field,
            "run_name": spec.run_name,
            "nfp": case.nfp,
            "threads": spec.threads,
            "t_load": t_load,
            "t_mgrid": t_mgrid,
            "t_case_metrics": t_case_metrics,
            **case_info,
            "ftol": level.ftol,
            "ns_target": level.ns_final,
            "mpol": level.mpol,
            "ntor": level.ntor,
            "nzeta": level.nzeta,
            "niter_budget": level.niter,
        }
        t0 = time.time()
        vmec_input = case_lib.make_input(case, regime, level)
        rec["t_input"] = time.time() - t0
        if level.condense:
            data = case.condensed[level.condense_epsilon][1]
            rec["condense"] = {
                "spectral_width_reduction": data.get("spectral_width_reduction"),
                "max_RZ_error": data.get("max_RZ_error"),
                "initial_objective": data.get("initial_objective"),
                "final_objective": data.get("final_objective"),
            }
        rec["pres_scale"] = float(vmec_input.pres_scale)
        rec["curtor"] = float(vmec_input.curtor)
        rec["target_beta"] = case_lib.target_beta(case, regime)

        exception: str | None = None
        output = None
        t0 = time.time()
        set_omp_threads(spec.threads)
        with FdCapture() as cap:
            try:
                output = vmecpp.run(
                    vmec_input,
                    magnetic_field=table if regime.free_boundary else None,
                    verbose=1,
                    max_threads=spec.threads,
                )
            except Exception as exc:  # noqa: BLE001
                exception = f"{type(exc).__name__}: {_first_line(str(exc))}"
        rec["t_run"] = time.time() - t0
        rec["exception"] = exception
        rec["log"] = parse_log(cap.text)
        set_omp_threads(1)
        print(  # noqa: T201
            f"{rec['key']}: run {rec['t_run']:.1f} s, exception={exception!r}",
            file=sys.stderr,
            flush=True,
        )

        converged = False
        fsq_final: float | None = None
        tail: dict | None = None
        ier_flag: int | None = None
        if output is not None:
            wout = output.wout
            ier_flag = int(wout.ier_flag)
            fsqr, fsqz, fsql = float(wout.fsqr), float(wout.fsqz), float(wout.fsql)
            fsq_final = fsqr + fsqz + fsql
            converged = (
                ier_flag == 0
                and int(wout.ns) == level.ns_final
                and fsqr <= level.ftol
                and fsqz <= level.ftol
                and fsql <= level.ftol
            )
            fsqt = np.asarray(wout.fsqt, dtype=float)
            niter_last = int(wout.niter)
            tail = residual_tail_summary(fsqt[-niter_last:], level.ftol, niter_last)
            reasons = np.asarray(wout.restart_reason_timetrace, dtype=int)
            names = {int(r.value): r.name for r in vmecpp.RestartReason}
            restart_counts = {
                names.get(int(v), str(int(v))): int(c)
                for v, c in zip(*np.unique(reasons, return_counts=True), strict=True)
                if int(v) != 1  # 1 is "no restart"
            }
            step = max(1, len(fsqt) // 400)
            rec["fsq_trace"] = [
                round(float(np.log10(v)), 3) if np.isfinite(v) and v > 0 else None
                for v in fsqt[::step]
            ]
            rec["fsq_trace_step"] = step
            rec.update(
                {
                    "ier_flag": ier_flag,
                    "reason": wout.reason,
                    "ns_reached": int(wout.ns),
                    "niter_last_stage": niter_last,
                    "itfsq_total": int(wout.itfsq),
                    "fsqr": fsqr,
                    "fsqz": fsqz,
                    "fsql": fsql,
                    "fsq_min_last_stage": float(np.nanmin(fsqt[-niter_last:]))
                    if niter_last > 0 and len(fsqt)
                    else float("nan"),
                    "restart_counts": restart_counts,
                    "tail": tail,
                    "t_per_iter_ms": 1000.0 * rec["t_run"] / max(int(wout.itfsq), 1),
                    **metrics.physics_summary(wout),
                }
            )
            if converged or (ier_flag == 0 and int(wout.ns) == level.ns_final):
                t0 = time.time()
                rec["lcfs"] = {
                    "xm": np.asarray(wout.xm, dtype=int).tolist(),
                    "xn": np.asarray(wout.xn, dtype=int).tolist(),
                    "rmnc": np.asarray(wout.rmnc)[:, -1].tolist(),
                    "zmns": np.asarray(wout.zmns)[:, -1].tolist(),
                    "raxis_cc": np.asarray(wout.raxis_cc).tolist(),
                    "zaxis_cs": np.asarray(wout.zaxis_cs).tolist(),
                }
                try:
                    surface = metrics.lcfs_surface(wout)
                    rec.update(metrics.coil_distance(surface, case))
                    if (
                        case_lib.target_beta(case, regime) == 0.0
                        and regime.current_fraction == 0.0
                    ):
                        rec.update(metrics.normal_field_error(bs, surface))
                        rec.update(metrics.b2_mismatch(wout, bs))
                        rec.update(metrics.quasr_boundary_comparison(wout, case))
                        if spec.trace and converged and regime.free_boundary:
                            rec.update(
                                metrics.fieldline_deviation(
                                    bs,
                                    wout,
                                    n_lines=spec.trace_lines,
                                    n_turns=spec.trace_turns,
                                    tol=spec.trace_tol,
                                )
                            )
                            rec.update(
                                metrics.quasr_surface_confinement(
                                    bs,
                                    case,
                                    float(wout.Aminor_p),
                                    n_lines=spec.trace_lines,
                                    tol=spec.trace_tol,
                                )
                            )
                except Exception:  # noqa: BLE001
                    rec["metrics_error"] = _first_line(
                        traceback.format_exc().splitlines()[-1]
                    )
                rec["t_metrics"] = time.time() - t0
        rec["converged"] = bool(converged)
        rec["status"] = classify(
            converged=converged,
            exception=exception,
            log_text=cap.text,
            ier_flag=ier_flag,
            fsq_final=fsq_final,
            ftol=level.ftol,
            tail=tail,
            niter_budget=level.niter,
        )
        results.append(rec)
    return results


def worker_main(spec_path: str, out_path: str) -> int:
    spec = JobSpec(**json.loads(Path(spec_path).read_text()))
    results = run_worker(spec)
    Path(out_path).write_text(json.dumps(results, default=_json_default))
    return 0


def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    msg = f"not serialisable: {type(obj)}"
    raise TypeError(msg)


# ---------------------------------------------------------------------------
# Re-tracing field lines from stored boundaries
# ---------------------------------------------------------------------------


def retrace(
    results_dir: Path,
    *,
    n_lines: int = 8,
    n_turns: int = 100,
    tol: float = 1e-9,
    only_missing: bool = False,
    only_confinement: bool = False,
    ids: set[int] | None = None,
    cache_dir: Path | None = None,
    log=print,
) -> int:
    """Recompute the field-line metrics of every converged free-boundary vacuum
    record from its stored boundary, without re-running the solver.

    Useful after a change to the tracing (or to trace more turns), and to make
    the field-line columns of a ledger uniform. Rewrites ``results.jsonl``.
    """
    results_dir = Path(results_dir)
    path = results_dir / "results.jsonl"
    records = load_results(results_dir)
    set_omp_threads(1)
    fields: dict[int, tuple] = {}
    n_done = 0
    for rec in records:
        if (
            rec.get("regime") != "vacuum"
            or not rec.get("converged")
            or "lcfs" not in rec
        ):
            continue
        if only_missing and np.isfinite(rec.get("fl_dev_rms", np.nan)):
            continue
        case_id = int(rec["case_id"])
        if ids is not None and case_id not in ids:
            continue
        if only_confinement and np.isfinite(rec.get("fl_lost_quasr", np.nan)):
            continue
        if case_id not in fields:
            case = case_lib.load_case(case_id, cache_dir)
            fields[case_id] = (metrics.reference_field(case)[0], case)
        bs, case = fields[case_id]
        if not only_confinement:
            rec.update(
                metrics.fieldline_deviation(
                    bs,
                    metrics.lcfs_namespace(rec),
                    n_lines=n_lines,
                    n_turns=n_turns,
                    tol=tol,
                )
            )
        rec.update(
            metrics.quasr_surface_confinement(
                bs, case, float(rec["Aminor_p"]), n_lines=n_lines, tol=tol
            )
        )
        n_done += 1
        log(
            f"retraced {rec['key']}: dev rms {rec.get('fl_dev_rms', np.nan):.2e} a, "
            f"lost {rec.get('fl_lost_fraction', np.nan):.2f}, "
            f"lost from QUASR surface {rec.get('fl_lost_quasr', np.nan):.2f}"
        )
    tmp = path.with_suffix(".jsonl.tmp")
    with open(tmp, "w") as f:
        for rec in records:
            f.write(json.dumps(rec, default=_json_default) + "\n")
    tmp.replace(path)
    return n_done


# ---------------------------------------------------------------------------
# Matrix driver
# ---------------------------------------------------------------------------


def solver_metadata() -> dict:
    meta = {
        "vmecpp_version": getattr(vmecpp, "__version__", None),
        "vmecpp_file": getattr(vmecpp, "__file__", None),
        "simsopt_version": getattr(simsopt, "__version__", None),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
    }
    try:
        meta["repo_describe"] = subprocess.run(
            ["git", "describe", "--tags", "--always", "--dirty"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip()
    except OSError:
        meta["repo_describe"] = None
    return meta


def load_results(results_dir: Path) -> list[dict]:
    path = Path(results_dir) / "results.jsonl"
    if not path.exists():
        return []
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


def run_matrix(
    manifest,
    *,
    levels: list[str],
    fields: list[str],
    regimes: list[str],
    out_dir: Path,
    workers: int = 3,
    threads: int = 3,
    timeout: float = 1800.0,
    trace: bool = True,
    trace_lines: int = 8,
    trace_turns: int = 100,
    trace_tol: float = 1e-9,
    cache_dir: Path | None = None,
    run_name: str = "",
    redo: tuple[str, ...] = (),
    log=print,
) -> None:
    out_dir = Path(out_dir)
    jobs_dir = out_dir / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)
    results_path = out_dir / "results.jsonl"
    existing = load_results(out_dir)
    if redo and existing:
        # Drop the records of the classes to redo before they are re-run. Only
        # safe while no other run is appending to this ledger.
        kept = [r for r in existing if r.get("status") not in redo]
        tmp = results_path.with_suffix(".jsonl.tmp")
        with open(tmp, "w") as f:
            f.writelines(json.dumps(rec, default=_json_default) + "\n" for rec in kept)
        tmp.replace(results_path)
        log(f"dropped {len(existing) - len(kept)} records with status in {redo}")
        existing = kept
    done_keys = {r["key"] for r in existing}

    meta_path = out_dir / "meta.json"
    meta = {
        "run_name": run_name or out_dir.name,
        "started": time.strftime("%Y-%m-%d %H:%M:%S"),
        "levels": levels,
        "fields": fields,
        "regimes": regimes,
        "workers": workers,
        "threads": threads,
        "timeout": timeout,
        "trace": trace,
        "trace_lines": trace_lines,
        "trace_turns": trace_turns,
        "trace_tol": trace_tol,
        "n_cases": len(manifest),
        "presets": presets_as_dict(),
        **solver_metadata(),
    }
    if meta_path.exists():
        previous = json.loads(meta_path.read_text())
        meta["previous_runs"] = [
            *previous.get("previous_runs", []),
            previous.get("started"),
        ]
    meta_path.write_text(json.dumps(meta, indent=2))

    specs: list[JobSpec] = []
    meta_columns = [c for c in manifest.columns if c != "ID"]
    for _, row in manifest.iterrows():
        case_id = int(row["ID"])
        case_meta = {
            c: (row[c].item() if hasattr(row[c], "item") else row[c])
            for c in meta_columns
        }
        for level in levels:
            for field_name in fields:
                pending = [
                    r
                    for r in regimes
                    if result_key(case_id, r, level, field_name) not in done_keys
                ]
                if not pending:
                    continue
                specs.append(
                    JobSpec(
                        case_id=int(case_id),
                        level=level,
                        field=field_name,
                        regimes=pending,
                        threads=threads,
                        trace=trace,
                        trace_lines=trace_lines,
                        trace_turns=trace_turns,
                        trace_tol=trace_tol,
                        cache_dir=str(cache_dir or case_lib.default_cache_dir()),
                        work_dir=str(out_dir / "work"),
                        run_name=meta["run_name"],
                        meta=case_meta,
                    )
                )
    log(f"{len(specs)} jobs to run ({len(done_keys)} results already present)")

    lock = threading.Lock()
    counts: dict[str, int] = {}
    n_done = 0
    t_start = time.time()

    def record(records: list[dict]) -> None:
        nonlocal n_done
        with lock, open(results_path, "a") as f:
            for rec in records:
                f.write(json.dumps(rec, default=_json_default) + "\n")
                counts[rec["status"]] = counts.get(rec["status"], 0) + 1
            n_done += 1
        elapsed = time.time() - t_start
        summary = ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
        log(f"[{n_done}/{len(specs)} jobs, {elapsed / 60:.1f} min] {summary}")

    def run_one(spec: JobSpec) -> None:
        spec_path = jobs_dir / f"{spec.key}.json"
        out_path = jobs_dir / f"{spec.key}.out.json"
        log_path = jobs_dir / f"{spec.key}.log"
        spec_path.write_text(json.dumps(asdict(spec)))
        if out_path.exists():
            out_path.unlink()
        env = dict(os.environ)
        env["OMP_NUM_THREADS"] = str(spec.threads)
        cmd = [
            sys.executable,
            "-m",
            "benchmarks.free_boundary",
            "worker",
            str(spec_path),
            str(out_path),
        ]
        timed_out = False
        t0 = time.time()
        with open(log_path, "w") as logf:
            try:
                proc = subprocess.run(
                    cmd,
                    cwd=REPO_ROOT,
                    env=env,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    timeout=timeout,
                    check=False,
                )
                returncode = proc.returncode
            except subprocess.TimeoutExpired:
                timed_out = True
                returncode = -1
        wall = time.time() - t0
        if out_path.exists() and not timed_out:
            records = json.loads(out_path.read_text())
        else:
            tail = ""
            if log_path.exists():
                tail = "\n".join(log_path.read_text().splitlines()[-15:])[-2000:]
            records = [
                {
                    "key": result_key(spec.case_id, r, spec.level, spec.field),
                    "case_id": spec.case_id,
                    "regime": r,
                    "level": spec.level,
                    "field": spec.field,
                    "run_name": spec.run_name,
                    "threads": spec.threads,
                    "converged": False,
                    "status": "timeout" if timed_out else "error_other",
                    "exception": (
                        f"job timeout after {timeout:.0f} s"
                        if timed_out
                        else f"worker exit code {returncode}: {tail}"
                    ),
                    "t_run": wall,
                }
                for r in spec.regimes
            ]
        record(records)

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(run_one, specs))
    meta["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    meta["status_counts"] = counts
    meta_path.write_text(json.dumps(meta, indent=2))
