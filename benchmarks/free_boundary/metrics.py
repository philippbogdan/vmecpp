# SPDX-FileCopyrightText: 2024-present Proxima Fusion GmbH <info@proximafusion.com>
#
# SPDX-License-Identifier: MIT
"""Accuracy metrics for a converged vacuum free-boundary solution.

In vacuum with zero net current the plasma carries no current, so the exact
solution is a flux surface of the coil field alone. That gives ground truth
independent of VMEC++: the coil field is evaluated by Biot-Savart on the
SIMSOPT coil curves (spectrally accurate quadrature), never through the mgrid
table. Three views of the same question "is the VMEC++ boundary a flux surface
of the coil field?" are measured:

* ``normal_field_error``: B . n / |B| on the boundary (local, differential),
* ``fieldline_deviation``: distance between field lines started on the
  boundary and the boundary itself after many toroidal turns (global,
  integrated), which is the "geometry against traced field lines" test,
* ``b2_mismatch``: |B|^2 on the plasma side (VMEC++) against the coil field.

``mgrid_floor`` measures how well the response table VMEC++ actually receives
represents the coil field on the boundary. No solver can be more accurate than
that floor, so every accuracy number should be read next to it.
"""

from __future__ import annotations

import time
import types
import typing

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree
from simsopt.field import BiotSavart, coils_via_symmetries
from simsopt.field.tracing import (
    MaxRStoppingCriterion,
    MaxZStoppingCriterion,
    MinRStoppingCriterion,
    MinZStoppingCriterion,
    ToroidalTransitStoppingCriterion,
    compute_fieldlines,
)
from simsopt.geo import SurfaceRZFourier

from .cases import QuasrCase, resampled_curve

# ---------------------------------------------------------------------------
# Geometry from the wout
# ---------------------------------------------------------------------------


def lcfs_rz(wout, theta: np.ndarray, zeta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """R and Z of the last closed flux surface at the given angles (flattened)."""
    xm = np.asarray(wout.xm, dtype=float)
    xn = np.asarray(wout.xn, dtype=float)
    angle = np.outer(theta, xm) - np.outer(zeta, xn)
    r = np.cos(angle) @ np.asarray(wout.rmnc)[:, -1]
    z = np.sin(angle) @ np.asarray(wout.zmns)[:, -1]
    return r, z


def lcfs_surface(
    wout, n_phi_per_period: int = 16, n_theta: int = 64
) -> SurfaceRZFourier:
    """The VMEC++ boundary as a SIMSOPT surface over the full torus."""
    nfp = int(wout.nfp)
    xm = np.asarray(wout.xm, dtype=int)
    xn = np.asarray(wout.xn, dtype=int)
    mpol = int(xm.max())
    ntor = int(np.abs(xn).max() // nfp) if len(xn) else 0
    surface = SurfaceRZFourier(
        nfp=nfp,
        stellsym=True,
        mpol=mpol,
        ntor=ntor,
        quadpoints_phi=np.linspace(0.0, 1.0, n_phi_per_period * nfp, endpoint=False),
        quadpoints_theta=np.linspace(0.0, 1.0, n_theta, endpoint=False),
    )
    rmnc = np.asarray(wout.rmnc)[:, -1]
    zmns = np.asarray(wout.zmns)[:, -1]
    for i in range(len(xm)):
        m = int(xm[i])
        n = int(xn[i] // nfp)
        surface.set_rc(m, n, float(rmnc[i]))
        if m > 0 or n > 0:
            surface.set_zs(m, n, float(zmns[i]))
    return surface


def lcfs_namespace(record: dict) -> types.SimpleNamespace:
    """A minimal wout-like object rebuilt from a ledger record's ``lcfs`` entry,
    enough for every geometry function in this module."""
    lcfs = record["lcfs"]
    return types.SimpleNamespace(
        nfp=int(record["nfp"]),
        Aminor_p=float(record["Aminor_p"]),
        Rmajor_p=float(record["Rmajor_p"]),
        xm=np.asarray(lcfs["xm"], dtype=int),
        xn=np.asarray(lcfs["xn"], dtype=int),
        rmnc=np.asarray(lcfs["rmnc"], dtype=float)[:, None],
        zmns=np.asarray(lcfs["zmns"], dtype=float)[:, None],
        raxis_cc=np.asarray(lcfs["raxis_cc"], dtype=float),
        zaxis_cs=np.asarray(lcfs["zaxis_cs"], dtype=float),
    )


def cross_section_rz(
    wout, phi: float, n_theta: int = 2048
) -> tuple[np.ndarray, np.ndarray]:
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    return lcfs_rz(wout, theta, np.full(n_theta, phi))


def section_centroid(
    wout, phi: np.ndarray, n_theta: int = 256
) -> tuple[np.ndarray, np.ndarray]:
    """Area centroid of the boundary cross-section at each toroidal angle: a
    reference point inside the section for winding counts that does not depend
    on the axis Fourier convention."""
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    r_c = np.empty(len(phi))
    z_c = np.empty(len(phi))
    for i, p in enumerate(phi):
        r, z = lcfs_rz(wout, theta, np.full(n_theta, p))
        r2, z2 = np.roll(r, -1), np.roll(z, -1)
        cross = r * z2 - r2 * z
        area = 0.5 * np.sum(cross)
        if abs(area) < 1e-14:
            r_c[i], z_c[i] = r.mean(), z.mean()
        else:
            r_c[i] = np.sum((r + r2) * cross) / (6.0 * area)
            z_c[i] = np.sum((z + z2) * cross) / (6.0 * area)
    return r_c, z_c


def magnetic_axis_rz(wout, phi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = np.arange(len(wout.raxis_cc))
    arg = np.outer(phi, n * wout.nfp)
    r = np.cos(arg) @ np.asarray(wout.raxis_cc)
    z = np.sin(arg) @ np.asarray(wout.zaxis_cs)
    return r, z


# ---------------------------------------------------------------------------
# Coil field
# ---------------------------------------------------------------------------


def exact_field(case: QuasrCase, n_points: int | None = None) -> BiotSavart:
    """Biot-Savart field of the coil set.

    With ``n_points`` the base curves are re-quadratured (the QUASR curves carry
    160 points; the trapezoidal rule on a periodic curve converges
    exponentially, so this only matters very close to a coil).
    """
    if n_points is None:
        return BiotSavart(case.coils)
    curves = [resampled_curve(c.curve, n_points) for c in case.base_coils]
    currents = [c.current for c in case.base_coils]
    return BiotSavart(coils_via_symmetries(curves, currents, case.nfp, True))


def reference_field(
    case: QuasrCase, n_points: int = 640, tolerance: float = 1e-8
) -> tuple[BiotSavart, float, int]:
    """The ground-truth field for the metrics, with its self-check.

    The stored coil curves carry 160 quadrature points, which is spectrally
    accurate except very close to a coil; the re-quadratured rebuild agrees
    with a further doubling to machine precision. When the two differ by more
    than ``tolerance`` on the QUASR boundary the rebuild is used. Returns the
    field, the maximum relative difference, and the point count used.
    """
    check = exact_field_check(case, n_points)
    if check > tolerance:
        return exact_field(case, n_points), check, n_points
    return exact_field(case), check, len(case.base_coils[0].curve.quadpoints)


def exact_field_check(case: QuasrCase, n_points: int = 640) -> float:
    """Max relative difference between the stored coil set and a re-quadratured
    symmetric rebuild, on the QUASR boundary: validates both the quadrature and
    the assumption that the coil list is the symmetric expansion of the base
    coils."""
    xyz = case.boundary.gamma().reshape(-1, 3)
    a = exact_field(case)
    b = exact_field(case, n_points)
    a.set_points(xyz)
    b.set_points(xyz)
    diff = np.linalg.norm(a.B() - b.B(), axis=1)
    return float(diff.max() / np.linalg.norm(a.B(), axis=1).max())


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def normal_field_error(bs: BiotSavart, surface) -> dict[str, float]:
    """Area-weighted RMS and maximum of B . n / |B| on a surface."""
    normal = surface.normal().reshape(-1, 3)
    area = np.linalg.norm(normal, axis=1)
    unit_normal = normal / area[:, None]
    xyz = surface.gamma().reshape(-1, 3)
    bs.set_points(xyz)
    b = bs.B()
    bn = np.sum(b * unit_normal, axis=1)
    bmag = np.linalg.norm(b, axis=1)
    w = area / area.sum()
    return {
        "bn_rms": float(np.sqrt(np.sum(w * bn**2) / np.sum(w * bmag**2))),
        "bn_max": float(np.max(np.abs(bn) / bmag)),
    }


def b2_mismatch(wout, bs: BiotSavart, n_theta: int = 64, n_zeta_per_period: int = 16):
    """|B|^2 on the plasma side of the boundary (VMEC++, extrapolated from the
    two outermost half-grid surfaces) against the coil field.

    The linear extrapolation over half a radial cell adds an error of order
    (1/ns)^2 that is not the solver's; treat this as a secondary metric.
    """
    nfp = int(wout.nfp)
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi / nfp, n_zeta_per_period, endpoint=False)
    tt, zz = np.meshgrid(theta, zeta, indexing="ij")
    tt = tt.ravel()
    zz = zz.ravel()
    r, z = lcfs_rz(wout, tt, zz)
    xm = np.asarray(wout.xm_nyq, dtype=float)
    xn = np.asarray(wout.xn_nyq, dtype=float)
    bmnc = np.asarray(wout.bmnc)
    b_edge = 1.5 * bmnc[:, -1] - 0.5 * bmnc[:, -2]
    b_vmec = np.cos(np.outer(tt, xm) - np.outer(zz, xn)) @ b_edge
    xyz = np.stack([r * np.cos(zz), r * np.sin(zz), z], axis=1)
    bs.set_points(xyz)
    b_coil = np.linalg.norm(bs.B(), axis=1)
    rel = (b_vmec**2 - b_coil**2) / b_coil**2
    return {
        "b2_rms": float(np.sqrt(np.mean(rel**2))),
        "b2_max": float(np.max(np.abs(rel))),
    }


def _point_to_polyline_distance(points: np.ndarray, polyline: np.ndarray) -> np.ndarray:
    """Distance from each 2-D point to a closed polyline (segment-exact)."""
    a = polyline
    b = np.roll(polyline, -1, axis=0)
    ab = b - a
    ab2 = np.maximum(np.sum(ab**2, axis=1), 1e-300)
    out = np.empty(len(points))
    for i, p in enumerate(points):
        t = np.clip(np.sum((p - a) * ab, axis=1) / ab2, 0.0, 1.0)
        proj = a + t[:, None] * ab
        out[i] = np.min(np.linalg.norm(p - proj, axis=1))
    return out


def fieldline_deviation(
    bs: BiotSavart,
    wout,
    n_lines: int = 8,
    n_turns: int = 100,
    tol: float = 1e-9,
) -> dict[str, float]:
    """Trace field lines of the coil field from the VMEC++ boundary and measure
    how far their phi = 0 crossings stray from the boundary cross-section.

    Distances are in units of the minor radius ``Aminor_p``. A line is "lost"
    when it leaves the tracing box or strays by more than half a minor radius;
    lost lines are counted but excluded from the RMS and maximum. The traced
    rotational transform is estimated from the same lines using the VMEC++
    magnetic axis as the poloidal reference.
    """
    nfp = int(wout.nfp)
    a_minor = float(wout.Aminor_p)
    r_major = float(wout.Rmajor_p)
    r_sec, z_sec = cross_section_rz(wout, 0.0)
    polyline = np.stack([r_sec, z_sec], axis=1)
    theta0 = 2.0 * np.pi * (np.arange(n_lines) + 0.5) / n_lines
    r0, z0 = lcfs_rz(wout, theta0, np.zeros(n_lines))

    # The tracing box is the extent of the whole surface (the axis of an nfp = 1
    # configuration can wander by more than the major radius) plus a margin.
    tt, zz = np.meshgrid(
        np.linspace(0.0, 2.0 * np.pi, 64, endpoint=False),
        np.linspace(0.0, 2.0 * np.pi / nfp, 32, endpoint=False),
        indexing="ij",
    )
    r_all, z_all = lcfs_rz(wout, tt.ravel(), zz.ravel())
    stopping = [
        ToroidalTransitStoppingCriterion(n_turns, False),
        MinRStoppingCriterion(max(float(r_all.min()) - 3.0 * a_minor, 1e-3)),
        MaxRStoppingCriterion(float(r_all.max()) + 3.0 * a_minor),
        MinZStoppingCriterion(float(z_all.min()) - 3.0 * a_minor),
        MaxZStoppingCriterion(float(z_all.max()) + 3.0 * a_minor),
    ]
    # SIMSOPT integrates dx/dt = B, so t is arc length divided by |B|; the
    # transit criterion ends the trace, tmax is only a safety net.
    bs.set_points(np.stack([r0, np.zeros(n_lines), z0], axis=1))
    b_min = max(float(np.linalg.norm(bs.B(), axis=1).min()), 1e-3)
    # Centroid of the section on a fine toroidal grid, interpolated at the hits.
    phi_grid = np.linspace(0.0, 2.0 * np.pi, 256, endpoint=False)
    rc_grid, zc_grid = section_centroid(wout, phi_grid)

    def centroid_at(phi_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        pg = np.append(phi_grid, 2.0 * np.pi)
        return (
            np.interp(phi_values, pg, np.append(rc_grid, rc_grid[0])),
            np.interp(phi_values, pg, np.append(zc_grid, zc_grid[0])),
        )

    t0 = time.time()
    # Crossings are recorded on n_planes toroidal planes: plane 0 gives the
    # geometric deviation, all planes give a sampling of the line fine enough
    # (less than pi of poloidal angle between samples at any iota met here)
    # to count its winding about the magnetic axis.
    n_planes = 64
    _, hits = compute_fieldlines(
        bs,
        list(r0),
        list(z0),
        tmax=8.0 * np.pi * r_major * (n_turns + 1) / b_min,
        tol=tol,
        phis=list(2.0 * np.pi * np.arange(n_planes) / n_planes),
        stopping_criteria=stopping,
    )
    seconds = time.time() - t0

    devs: list[np.ndarray] = []
    iotas: list[float] = []
    lost = 0
    for raw_hits in hits:
        hit = np.asarray(raw_hits)
        crossings = hit[hit[:, 1] == 0] if len(hit) else hit
        if len(crossings) < 0.9 * n_turns:
            lost += 1
            continue
        pts = np.stack(
            [np.hypot(crossings[:, 2], crossings[:, 3]), crossings[:, 4]], axis=1
        )
        d = _point_to_polyline_distance(pts, polyline) / a_minor
        if d.max() > 0.5:
            lost += 1
            continue
        devs.append(d)
        # Rotational transform from the winding about the section centroid,
        # sampled at every plane crossing in time order.
        xyz = hit[np.argsort(hit[:, 0])][:, 2:5]
        phi = np.unwrap(np.arctan2(xyz[:, 1], xyz[:, 0]))
        r_ax, z_ax = centroid_at(np.mod(phi, 2.0 * np.pi))
        theta = np.unwrap(
            np.arctan2(xyz[:, 2] - z_ax, np.hypot(xyz[:, 0], xyz[:, 1]) - r_ax)
        )
        if abs(phi[-1] - phi[0]) > 0:
            iotas.append(float((theta[-1] - theta[0]) / (phi[-1] - phi[0])))
    all_devs = np.concatenate(devs) if devs else np.array([np.nan])
    return {
        "fl_dev_rms": float(np.sqrt(np.mean(all_devs**2))),
        "fl_dev_max": float(np.max(all_devs)),
        "fl_lost_fraction": lost / n_lines,
        "fl_iota": float(np.mean(iotas)) if iotas else float("nan"),
        "fl_turns": n_turns,
        "fl_seconds": seconds,
        "fl_nfp": nfp,
    }


def _safe_cross_section(surface, phi_fraction: float, n_theta: int):
    """SIMSOPT's cross_section fails for surfaces whose cylindrical angle is not
    monotonic along the toroidal parameter (strongly non-axisymmetric nfp = 1
    shapes); return None in that case so callers can report NaN instead."""
    try:
        return surface.cross_section(phi_fraction, thetas=n_theta)
    except Exception:  # noqa: BLE001
        return None


def quasr_surface_confinement(
    bs: BiotSavart,
    case: QuasrCase,
    a_minor: float,
    n_lines: int = 8,
    n_turns: int = 30,
    tol: float = 1e-9,
) -> dict[str, float]:
    """Fraction of field lines started on the QUASR outermost surface that
    leave it by more than half a minor radius within ``n_turns`` turns.

    A configuration whose own surface does not confine the coil field's lines
    has its last closed flux surface inside the QUASR boundary (an optimisation
    residual of the database); "lost" lines from the VMEC++ boundary of such a
    case say nothing about the solver.
    """
    sec = _safe_cross_section(case.boundary, 0.0, 2048)
    if sec is None:
        return {"fl_lost_quasr": float("nan")}
    polyline = np.stack([np.hypot(sec[:, 0], sec[:, 1]), sec[:, 2]], axis=1)
    starts = _safe_cross_section(case.boundary, 0.0, n_lines)
    if starts is None:
        return {"fl_lost_quasr": float("nan")}
    r0 = np.hypot(starts[:, 0], starts[:, 1])
    z0 = starts[:, 2]
    r_lo, r_hi = float(polyline[:, 0].min()), float(polyline[:, 0].max())
    z_hi = float(np.abs(polyline[:, 1]).max())
    stopping = [
        ToroidalTransitStoppingCriterion(n_turns, False),
        MinRStoppingCriterion(max(r_lo - 3.0 * a_minor, 1e-3)),
        MaxRStoppingCriterion(r_hi + 3.0 * a_minor),
        MinZStoppingCriterion(-(z_hi + 3.0 * a_minor)),
        MaxZStoppingCriterion(z_hi + 3.0 * a_minor),
    ]
    bs.set_points(np.stack([r0, np.zeros(n_lines), z0], axis=1))
    b_min = max(float(np.linalg.norm(bs.B(), axis=1).min()), 1e-3)
    _, hits = compute_fieldlines(
        bs,
        list(r0),
        list(z0),
        tmax=8.0 * np.pi * case.r_char * (n_turns + 1) / b_min,
        tol=tol,
        phis=[0.0],
        stopping_criteria=stopping,
    )
    lost = 0
    for raw in hits:
        hit = np.asarray(raw)
        if len(hit) < 0.9 * n_turns:
            lost += 1
            continue
        pts = np.stack([np.hypot(hit[:, 2], hit[:, 3]), hit[:, 4]], axis=1)
        if _point_to_polyline_distance(pts, polyline).max() / a_minor > 0.5:
            lost += 1
    return {"fl_lost_quasr": lost / n_lines}


def quasr_boundary_comparison(
    wout, case: QuasrCase, n_phi: int = 6, n_theta: int = 1024
) -> dict[str, float]:
    """Geometric distance between the VMEC++ boundary and the QUASR surface it
    started from, and the ratio of enclosed volumes.

    Both should coincide if both are exact flux surfaces of the same field with
    the same enclosed flux. The distance is the mean over cross-sections of the
    symmetric RMS (and the overall maximum) point-to-curve distance, in minor
    radii, against the full-resolution QUASR surface.
    """
    a_minor = float(wout.Aminor_p)
    rms: list[float] = []
    mx = 0.0
    for k in range(n_phi):
        phi_fraction = 0.5 * k / (n_phi * case.nfp)
        sec = _safe_cross_section(case.boundary, phi_fraction, n_theta)
        if sec is None:
            continue
        quasr = np.stack([np.hypot(sec[:, 0], sec[:, 1]), sec[:, 2]], axis=1)
        r, z = cross_section_rz(wout, 2.0 * np.pi * phi_fraction, n_theta)
        vmec = np.stack([r, z], axis=1)
        d1 = _point_to_polyline_distance(vmec[::4], quasr) / a_minor
        d2 = _point_to_polyline_distance(quasr[::4], vmec) / a_minor
        d = np.concatenate([d1, d2])
        rms.append(float(np.sqrt(np.mean(d**2))))
        mx = max(mx, float(d.max()))
    quasr_volume = abs(float(case.boundary.volume()))
    return {
        "quasr_dist_rms": float(np.mean(rms)) if rms else float("nan"),
        "quasr_dist_max": mx if rms else float("nan"),
        "volume_ratio_quasr": float(abs(wout.volume) / quasr_volume),
    }


def coil_distance(surface, case: QuasrCase) -> dict[str, float]:
    """Minimum distance from the boundary to any coil filament point [m]."""
    coil_points = np.concatenate([c.curve.gamma() for c in case.coils], axis=0)
    tree = cKDTree(coil_points)
    d, _ = tree.query(surface.gamma().reshape(-1, 3))
    return {"lcfs_coil_min_dist": float(d.min())}


def mgrid_floor(
    case: QuasrCase, table, nphi: int, bs: BiotSavart | None = None
) -> dict[str, float]:
    """Relative error of the mgrid response table against the exact coil field on
    the QUASR boundary, evaluated the way VMEC++ uses the table: exactly on the
    toroidal planes, bilinear in R and Z.

    This bounds the accuracy of any free-boundary solution built on this table;
    it contains the coil polygon error and the R-Z interpolation error.
    """
    bs = bs or exact_field(case)
    p = table.parameters
    nr, nz = p.number_of_r_grid_points, p.number_of_z_grid_points
    r = np.linspace(p.r_grid_minimum, p.r_grid_maximum, nr)
    z = np.linspace(p.z_grid_minimum, p.z_grid_maximum, nz)
    scale = float(case.extcur[0])
    comps = [
        np.asarray(getattr(table, name)[0]).reshape(nphi, nz, nr) * scale
        for name in ("b_r", "b_p", "b_z")
    ]
    errs: list[np.ndarray] = []
    mags: list[np.ndarray] = []
    for k in range(nphi):
        phi = k * (2.0 * np.pi / case.nfp) / nphi
        sec = _safe_cross_section(case.boundary, phi / (2.0 * np.pi), 200)
        if sec is None:
            continue
        rr = np.hypot(sec[:, 0], sec[:, 1])
        zz = sec[:, 2]
        pp = np.arctan2(sec[:, 1], sec[:, 0])
        bs.set_points(sec)
        b = bs.B()
        exact = np.stack(
            [
                b[:, 0] * np.cos(pp) + b[:, 1] * np.sin(pp),
                -b[:, 0] * np.sin(pp) + b[:, 1] * np.cos(pp),
                b[:, 2],
            ],
            axis=1,
        )
        interp = np.stack(
            [
                RegularGridInterpolator((z, r), comps[c][k], method="linear")(
                    np.stack([zz, rr], axis=1)
                )
                for c in range(3)
            ],
            axis=1,
        )
        errs.append(np.linalg.norm(interp - exact, axis=1))
        mags.append(np.linalg.norm(exact, axis=1))
    if not errs:
        nan = float("nan")
        return {
            "mgrid_floor_rms": nan,
            "mgrid_floor_max": nan,
            "mgrid_cell_over_a": nan,
        }
    err = np.concatenate(errs)
    mag = np.concatenate(mags)
    return {
        "mgrid_floor_rms": float(np.sqrt(np.mean(err**2)) / np.sqrt(np.mean(mag**2))),
        "mgrid_floor_max": float(err.max() / mag.max()),
        "mgrid_cell_over_a": float(
            (r[1] - r[0]) / (0.5 * (case.extent[1] - case.extent[0]))
        ),
    }


def physics_summary(wout) -> dict[str, typing.Any]:
    """Global quantities of a solution worth keeping in the ledger."""
    iotaf = np.asarray(wout.iotaf)
    delbsq = np.asarray(wout.delbsq)
    return {
        "volume": float(wout.volume),
        "aspect": float(wout.aspect),
        "Rmajor_p": float(wout.Rmajor_p),
        "Aminor_p": float(wout.Aminor_p),
        "b0": float(wout.b0),
        "betatotal": float(wout.betatotal),
        "betapol": float(wout.betapol),
        "betator": float(wout.betator),
        "ctor": float(wout.ctor),
        "iota_axis": float(iotaf[0]),
        "iota_edge": float(iotaf[-1]),
        "raxis_phi0": float(np.sum(wout.raxis_cc)),
        "delbsq_last": float(delbsq[-1]) if len(delbsq) else float("nan"),
    }
