# SPDX-FileCopyrightText: 2024-present Proxima Fusion GmbH <info@proximafusion.com>
#
# SPDX-License-Identifier: MIT
"""Failure taxonomy for a free-boundary run.

Every run ends in exactly one of these classes. The first group comes from the
solver's own status; the last group splits "ran out of iterations" by what the
force residual was doing at the end, which is what a solver developer wants to
know (a run that is still converging is a budget problem, a stalled one is a
method problem).
"""

from __future__ import annotations

import numpy as np

CLASSES = (
    "converged",
    "vacuum_grid_exceeded",
    "jacobian_75_times",
    "bad_jacobian",
    "early_failure",
    "nan_residual",
    "max_iter_near",
    "max_iter_slow",
    "max_iter_oscillating",
    "max_iter_stalled",
    "timeout",
    "error_other",
)

DESCRIPTIONS = {
    "converged": "all three force residuals below ftol at the final ns",
    "vacuum_grid_exceeded": "the boundary left the mgrid box (field clamped, run aborted)",
    "jacobian_75_times": "flux surfaces kept self-intersecting (Jacobian reset 75 times)",
    "bad_jacobian": "initial Jacobian bad even after the ns = 3 retry",
    "early_failure": "solver failed in the first iterations of a stage",
    "nan_residual": "non-finite force residual",
    "max_iter_near": "budget exhausted within a decade of ftol",
    "max_iter_slow": "budget exhausted while still converging (would reach ftol within 3x budget)",
    "max_iter_oscillating": "budget exhausted with the residual oscillating",
    "max_iter_stalled": "budget exhausted on a residual plateau",
    "timeout": "wall-clock limit of the bench hit",
    "error_other": "any other exception (see the message)",
}


def residual_tail_summary(
    fsq: np.ndarray, ftol: float, n_last: int
) -> dict[str, float]:
    """Behaviour of the total residual over the last stage of a run.

    ``fsq`` is the per-iteration residual of the last multigrid stage. The slope
    is a least-squares fit of log10(fsq) over its last quarter, in decades per
    iteration; ``iters_to_ftol`` extrapolates that slope from the last value.
    """
    fsq = np.asarray(fsq, dtype=float)
    fsq = fsq[np.isfinite(fsq) & (fsq > 0)]
    if len(fsq) < 8:
        return {
            "tail_slope": float("nan"),
            "tail_std": float("nan"),
            "iters_to_ftol": float("nan"),
        }
    k = max(8, len(fsq) // 4)
    tail = np.log10(fsq[-k:])
    x = np.arange(k, dtype=float)
    slope = float(np.polyfit(x, tail, 1)[0])
    std = float(np.std(tail - np.polyval(np.polyfit(x, tail, 1), x)))
    decades_left = float(np.log10(fsq[-1]) - np.log10(ftol))
    iters_to_ftol = decades_left / (-slope) if slope < 0 else float("inf")
    return {
        "tail_slope": slope,
        "tail_std": std,
        "iters_to_ftol": float(iters_to_ftol),
        "n_last": int(n_last),
    }


def classify(
    *,
    converged: bool,
    exception: str | None,
    log_text: str,
    ier_flag: int | None,
    fsq_final: float | None,
    ftol: float,
    tail: dict[str, float] | None,
    niter_budget: int,
    timed_out: bool = False,
) -> str:
    if timed_out:
        return "timeout"
    if converged:
        return "converged"
    text = (exception or "") + "\n" + (log_text or "")
    if "exceeded the vacuum field grid" in text:
        return "vacuum_grid_exceeded"
    if "JACOBIAN_75_TIMES_BAD" in text or ier_flag == 4:
        return "jacobian_75_times"
    if "failed during the first iterations" in text:
        return "early_failure"
    if ier_flag == 1 or "BAD_JACOBIAN" in text:
        return "bad_jacobian"
    if fsq_final is not None and not np.isfinite(fsq_final):
        return "nan_residual"
    if exception and ier_flag not in (0, None):
        return "error_other"
    if fsq_final is None:
        return "error_other"
    if fsq_final <= 10.0 * ftol:
        return "max_iter_near"
    if tail:
        if (
            np.isfinite(tail.get("iters_to_ftol", np.inf))
            and tail["iters_to_ftol"] < 3.0 * niter_budget
        ):
            return "max_iter_slow"
        if tail.get("tail_std", 0.0) > 0.3:
            return "max_iter_oscillating"
    return "max_iter_stalled"
