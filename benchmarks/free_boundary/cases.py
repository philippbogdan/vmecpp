# SPDX-FileCopyrightText: 2024-present Proxima Fusion GmbH <info@proximafusion.com>
#
# SPDX-License-Identifier: MIT
"""QUASR case set: database index, stratified selection, download, and the
conversion of a QUASR configuration into what a VMEC++ free-boundary run needs.

The QUASR database (https://quasr.flatironinstitute.org) ships stellarator
configurations as SIMSOPT serialisations that bundle nested flux surfaces with
the coil set that produces them. The outermost surface is only used as the
initial guess and as an independent reference; the coils define the problem.
"""

from __future__ import annotations

import dataclasses
import gzip
import json
import os
import time
import typing
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from simsopt._core import load as simsopt_load
from simsopt.field import BiotSavart, coils_to_makegrid
from simsopt.geo import CurveXYZFourier

import vmecpp

from .settings import FieldSpec, Level, Regime

VACUUM_PERMEABILITY = 4.0e-7 * np.pi  # mu0 [T m / A]

QUASR_BASE_URL = "https://quasr.flatironinstitute.org"
PACKAGE_DIR = Path(__file__).parent
REPO_ROOT = PACKAGE_DIR.parent.parent
REPO_QUASR_DIR = REPO_ROOT / "tests" / "data" / "quasr"
DEFAULT_MANIFEST = PACKAGE_DIR / "cases" / "manifest.csv"

# The twelve configurations checked into tests/data/quasr: the reference subset.
REFERENCE_IDS = (
    954,
    9914,
    19493,
    19609,
    19940,
    29346,
    50136,
    65579,
    112718,
    165868,
    167381,
    336902,
)

MANIFEST_COLUMNS = [
    "ID",
    "subset",
    "ladder",
    "nfp",
    "nc_per_hp",
    "helicity",
    "aspect_ratio",
    "minor_radius",
    "volume",
    "mean_iota",
    "iota_axis",
    "iota_edge",
    "iota_shear",
    "mean_elongation",
    "max_elongation",
    "max_kappa",
    "max_msc",
    "min_coil2surface_dist",
    "coil_dist_over_a",
    "min_coil2coil_dist",
    "qs_error",
    "Nsurfaces",
    "message",
]


def default_cache_dir() -> Path:
    env = os.environ.get("VMECPP_FB_BENCH_CACHE")
    if env:
        return Path(env)
    return Path.home() / ".cache" / "vmecpp_free_boundary_bench"


def quasr_serial_url(config_id: int) -> str:
    id_str = f"{config_id:07d}"
    return f"{QUASR_BASE_URL}/simsopt_serials/{id_str[0:4]}/serial{id_str}.json"


def _download(url: str, path: Path, attempts: int = 3) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "vmecpp-free-boundary-bench"}
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read()
            path.write_bytes(data)
            return
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(2.0 * (attempt + 1))
    msg = f"could not download {url}: {last_error}"
    raise RuntimeError(msg)


# ---------------------------------------------------------------------------
# Database index and case selection
# ---------------------------------------------------------------------------


def load_database(cache_dir: Path | None = None) -> pd.DataFrame:
    """The QUASR index (one row per configuration), downloaded once and cached."""
    cache_dir = cache_dir or default_cache_dir()
    path = cache_dir / "database.json.gz"
    if not path.exists():
        _download(f"{QUASR_BASE_URL}/database.json.gz", path)
    with gzip.open(path) as f:
        raw = json.load(f)
    df = pd.DataFrame(raw["data"], columns=raw["columns"]).reset_index(drop=True)
    return add_derived_columns(df)


def add_derived_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Derived quantities used for stratification and for the report."""
    df = df.copy()
    df["coil_dist_over_a"] = df["min_coil2surface_dist"] / df["minor_radius"]
    profiles = df["iota_profile"].map(lambda p: [float(v) for v in p])
    tf = df["tf_profile"].map(lambda p: [float(v) for v in p])
    df["iota_axis"] = profiles.map(lambda p: p[0])
    df["iota_edge"] = profiles.map(lambda p: p[-1])

    def shear(args: tuple[list[float], list[float]]) -> float:
        iota, flux = args
        if len(iota) < 2 or flux[-1] == flux[0]:
            return 0.0
        return (iota[-1] - iota[0]) / (flux[-1] - flux[0])

    df["iota_shear"] = list(map(shear, zip(profiles, tf, strict=True)))
    df["helicity"] = df["helicity"].astype(int)
    return df


def _bin_labels(values: pd.Series, edges: list[float], labels: list[str]) -> pd.Series:
    return pd.cut(values, bins=edges, labels=labels, include_lowest=True).astype(str)


def select_cases(
    df: pd.DataFrame,
    n_main: int = 300,
    n_ladder: int = 12,
    seed: int = 20260902,
    reference_ids: typing.Sequence[int] = REFERENCE_IDS,
) -> pd.DataFrame:
    """Stratified sample of the QUASR index.

    Strata are the product of the number of field periods, terciles of the
    coil-to-surface distance in minor radii (the "close to coils" axis) and
    terciles of the aspect ratio, each tercile taken over the whole database.
    Cases are drawn round-robin over the strata so every stratum is represented
    even where the database is sparse. The reference subset (the configurations
    checked into ``tests/data/quasr``) is always included, and ``n_ladder``
    further cases are flagged for the resolution ladder.
    """
    rng = np.random.default_rng(seed)
    pool = df[~df["ID"].isin(reference_ids)].reset_index(drop=True)

    dist_edges = [-np.inf, *pool["coil_dist_over_a"].quantile([1 / 3, 2 / 3]), np.inf]
    aspect_edges = [-np.inf, *pool["aspect_ratio"].quantile([1 / 3, 2 / 3]), np.inf]
    pool["dist_bin"] = _bin_labels(
        pool["coil_dist_over_a"], dist_edges, ["near", "mid", "far"]
    )
    pool["aspect_bin"] = _bin_labels(
        pool["aspect_ratio"], aspect_edges, ["low", "mid", "high"]
    )
    pool["nfp_bin"] = pool["nfp"].clip(upper=6).astype(int).astype(str)
    pool["stratum"] = (
        pool["nfp_bin"] + "/" + pool["dist_bin"] + "/" + pool["aspect_bin"]
    )

    # Shuffle within each stratum, then take round-robin until n_main reached.
    order = rng.permutation(len(pool))
    pool = pool.iloc[order].reset_index(drop=True)
    queues = {name: list(group.index) for name, group in pool.groupby("stratum")}
    names = sorted(queues)
    chosen: list = []
    while len(chosen) < n_main and any(queues.values()):
        for name in names:
            if queues[name] and len(chosen) < n_main:
                chosen.append(queues[name].pop())
    main = pool.loc[chosen].copy()
    main["subset"] = "main"

    reference = df[df["ID"].isin(reference_ids)].copy()
    reference["subset"] = "reference"
    reference["stratum"] = ""

    selected = pd.concat([reference, main], ignore_index=True)
    selected["ladder"] = selected["subset"] == "reference"
    # Ladder extras: round-robin over the main strata with a fresh permutation,
    # so the ladder subset spans distance and aspect-ratio bins as well.
    extras = main.iloc[rng.permutation(len(main))]
    seen: set = set()
    picked = 0
    for _, row in extras.iterrows():
        key = (row["dist_bin"], row["aspect_bin"])
        if key in seen and picked < n_ladder:
            continue
        seen.add(key)
        selected.loc[selected["ID"] == row["ID"], "ladder"] = True
        picked += 1
        if picked >= n_ladder:
            break
    selected = selected.sort_values(["subset", "ID"], ascending=[False, True])
    return selected[MANIFEST_COLUMNS].reset_index(drop=True)


def load_manifest(path: Path | str = DEFAULT_MANIFEST) -> pd.DataFrame:
    return pd.read_csv(path)


# ---------------------------------------------------------------------------
# Serial files and case loading
# ---------------------------------------------------------------------------


def serial_path(config_id: int, cache_dir: Path | None = None) -> Path:
    """Local path of a QUASR serial file, downloading it if needed.

    The copies checked into the repository (Git LFS) are used when present and
    actually fetched; otherwise the file is downloaded into the cache.
    """
    name = f"serial{config_id:07d}.json"
    repo_copy = REPO_QUASR_DIR / name
    if repo_copy.exists() and repo_copy.stat().st_size > 4096:
        return repo_copy
    cache_dir = cache_dir or default_cache_dir()
    local = cache_dir / "quasr" / name
    if not (local.exists() and local.stat().st_size > 4096):
        _download(quasr_serial_url(config_id), local)
    return local


@dataclass
class QuasrCase:
    """A QUASR configuration reduced to what a free-boundary run needs."""

    config_id: int
    nfp: int
    boundary: typing.Any  # simsopt surface at full torus (outermost QUASR surface)
    coils: list  # full coil set (all symmetry copies)
    base_coils: list  # the unique coils (first half period)
    extcur: np.ndarray  # one circuit: the first coil's current [A]
    phiedge: float  # enclosed vacuum toroidal flux [Wb]
    r_char: float  # characteristic major radius, R00 of the boundary
    b_char: float  # characteristic |B| on the boundary cross-section [T]
    extent: tuple[float, float, float, float]  # R_min, R_max, Z_min, Z_max
    meta: dict  # the manifest row (QUASR index values), possibly empty
    condensed: dict = dataclasses.field(default_factory=dict)  # epsilon -> surface

    @property
    def stellsym(self) -> bool:
        return bool(self.boundary.stellsym)


def boundary_extent(surface) -> tuple[float, float, float, float]:
    xyz = surface.gamma()
    r = np.hypot(xyz[..., 0], xyz[..., 1])
    z = xyz[..., 2]
    return float(r.min()), float(r.max()), float(z.min()), float(z.max())


def enclosed_toroidal_flux(coils, surface, n_theta: int = 1000) -> tuple[float, float]:
    """Enclosed vacuum toroidal flux (phiedge) and a characteristic |B|.

    The flux is the line integral of the coil vector potential around the
    boundary cross-section at phi = 0 (exact for the Biot-Savart field); the
    sign follows the toroidal field direction at the cross-section centre.
    """
    biot_savart = BiotSavart(coils)
    loop = surface.cross_section(0.0, thetas=n_theta)
    biot_savart.set_points(loop)
    a_field = biot_savart.A()
    b_char = float(np.mean(np.linalg.norm(biot_savart.B(), axis=1)))
    dl = np.roll(loop, -1, axis=0) - loop
    flux = float(np.sum(0.5 * (a_field + np.roll(a_field, -1, axis=0)) * dl))
    r_centre = float(np.hypot(*loop[:, :2].mean(axis=0)))
    biot_savart.set_points(np.array([[r_centre, 0.0, 0.0]]))
    sign = float(np.sign(biot_savart.B()[0, 1])) or 1.0
    return abs(flux) * sign, b_char


def load_case(
    config_id: int, cache_dir: Path | None = None, meta: dict | None = None
) -> QuasrCase:
    surfaces, coils = simsopt_load(str(serial_path(config_id, cache_dir)))
    boundary = surfaces[-1]
    nfp = int(boundary.nfp)
    if not boundary.stellsym:
        msg = f"QUASR {config_id}: only stellarator-symmetric cases are supported"
        raise ValueError(msg)
    n_base = len(coils) // (2 * nfp)
    phiedge, b_char = enclosed_toroidal_flux(coils, boundary)
    rz = boundary.to_RZFourier()
    return QuasrCase(
        config_id=config_id,
        nfp=nfp,
        boundary=boundary,
        coils=list(coils),
        base_coils=list(coils[:n_base]),
        # VMEC++ normalises a multi-coil circuit's response by its first coil's
        # current, so extcur must restore that reference current.
        extcur=np.array([float(coils[0].current.get_value())]),
        phiedge=phiedge,
        r_char=float(rz.get_rc(0, 0)),
        b_char=b_char,
        extent=boundary_extent(boundary),
        meta=dict(meta or {}),
    )


def resampled_curve(curve, n_points: int) -> CurveXYZFourier:
    """The same Fourier curve evaluated on ``n_points`` uniform quadrature points."""
    new = CurveXYZFourier(np.linspace(0.0, 1.0, n_points, endpoint=False), curve.order)
    new.x = curve.x
    return new


def write_coils_file(case: QuasrCase, path: Path, coil_points: int) -> Path:
    """A MAKEGRID coils file with all coils in one current group.

    Every coil is written as a closed polygon with ``coil_points`` vertices; the
    straight-segment error of that polygon is part of the vacuum-field floor
    measured by ``metrics.mgrid_floor``.
    """
    curves = [resampled_curve(c.curve, coil_points) for c in case.base_coils]
    currents = [c.current for c in case.base_coils]
    coils_to_makegrid(
        str(path),
        curves,
        currents,
        groups=[1] * len(case.coils),
        nfp=case.nfp,
        stellsym=True,
    )
    return path


def makegrid_parameters(
    case: QuasrCase, field: FieldSpec, nphi: int
) -> vmecpp.MakegridParameters:
    r_min, r_max, z_min, z_max = case.extent
    dr = r_max - r_min
    dz = z_max - z_min
    return vmecpp.MakegridParameters(
        normalize_by_currents=True,
        assume_stellarator_symmetry=True,
        number_of_field_periods=case.nfp,
        r_grid_minimum=r_min - field.margin * dr,
        r_grid_maximum=r_max + field.margin * dr,
        number_of_r_grid_points=field.points_rz,
        z_grid_minimum=z_min - field.margin * dz,
        z_grid_maximum=z_max + field.margin * dz,
        number_of_z_grid_points=field.points_rz,
        number_of_phi_grid_points=nphi,
    )


def build_response_table(
    case: QuasrCase, field: FieldSpec, nphi: int, work_dir: Path
) -> vmecpp.MagneticFieldResponseTable:
    work_dir.mkdir(parents=True, exist_ok=True)
    # One file per process: concurrent jobs of the same case must not share it.
    coils_file = write_coils_file(
        case,
        work_dir / f"coils.quasr{case.config_id:07d}_{field.name}_{os.getpid()}",
        field.coil_points,
    )
    params = makegrid_parameters(case, field, nphi)
    try:
        return vmecpp.MagneticFieldResponseTable.from_coils_file(coils_file, params)
    finally:
        coils_file.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# VmecInput assembly
# ---------------------------------------------------------------------------


def condensed_boundary(case: QuasrCase, epsilon: float):
    """The QUASR boundary with its poloidal angle reparametrised by SIMSOPT's
    spectral condensation (shape preserved to ``epsilon`` minor radii),
    computed once per case and cached on it."""
    if epsilon not in case.condensed:
        rz = case.boundary.to_RZFourier()
        surface, data = rz.condense_spectrum(epsilon=epsilon, verbose=False)
        case.condensed[epsilon] = (surface, data)
    return case.condensed[epsilon][0]


def boundary_coefficients(
    surface, mpol: int, ntor: int
) -> tuple[np.ndarray, np.ndarray, float]:
    """VMEC++ rbc/zbs arrays of shape (mpol, 2*ntor+1) from a SIMSOPT surface."""
    rz = surface.to_RZFourier()
    # SurfaceRZFourier uses m up to mpol inclusive, unlike VMEC++; resizing to
    # (mpol - 1, ntor) makes rz.rc/rz.zs exactly (mpol, 2*ntor+1).
    resized = rz.change_resolution(mpol - 1, ntor)
    rz = rz if resized is None else resized
    return rz.rc.copy(), rz.zs.copy(), float(rz.get_rc(0, 0))


def target_beta(case: QuasrCase, regime: Regime) -> float:
    """The volume-averaged beta a regime asks for on this case."""
    if regime.shift_target is None:
        return regime.target_beta
    iota = abs(float(case.meta.get("iota_edge", float("nan"))))
    aspect = float(case.meta.get("aspect_ratio", float("nan")))
    if not (np.isfinite(iota) and np.isfinite(aspect)) or aspect <= 0.0:
        msg = f"QUASR {case.config_id}: iota_edge and aspect_ratio needed for {regime.name}"
        raise ValueError(msg)
    return float(min(regime.beta_cap, regime.shift_target * 2.0 * iota**2 / aspect))


def pressure_scale(case: QuasrCase, regime: Regime) -> float:
    beta = target_beta(case, regime)
    if beta == 0.0:
        return 0.0
    return beta * case.b_char**2 / VACUUM_PERMEABILITY


def characteristic_current(case: QuasrCase) -> float:
    return 2.0 * np.pi * case.r_char * case.b_char / VACUUM_PERMEABILITY


def make_input(
    case: QuasrCase,
    regime: Regime,
    level: Level,
    *,
    free_boundary: bool | None = None,
    return_outputs_even_if_not_converged: bool = True,
) -> vmecpp.VmecInput:
    """Assemble the VmecInput for a case, regime and level.

    Mirrors ``tests/test_free_boundary_quasr.py``: the QUASR boundary is the
    initial guess, ``phiedge`` is the exact enclosed vacuum flux, pressure and
    current profiles are power series scaled to the requested beta and current.
    """
    if free_boundary is None:
        free_boundary = regime.free_boundary
    surface = (
        condensed_boundary(case, level.condense_epsilon)
        if level.condense
        else case.boundary
    )
    rbc, zbs, r_axis_guess = boundary_coefficients(surface, level.mpol, level.ntor)
    ns_array = np.asarray(level.ns_array, dtype=np.int64)

    vmec_input = vmecpp.VmecInput.default()
    vmec_input.lasym = False
    vmec_input.nfp = case.nfp
    vmec_input.mpol = level.mpol
    vmec_input.ntor = level.ntor
    vmec_input.ntheta = 0
    vmec_input.nzeta = level.nzeta
    vmec_input.ns_array = ns_array
    vmec_input.ftol_array = np.full(len(ns_array), level.ftol, dtype=float)
    vmec_input.niter_array = np.full(len(ns_array), level.niter, dtype=np.int64)
    vmec_input.nstep = 200
    vmec_input.delt = level.delt
    vmec_input.tcon0 = level.tcon0
    vmec_input.phiedge = case.phiedge
    vmec_input.gamma = 0.0
    vmec_input.bloat = 1.0

    vmec_input.pmass_type = "power_series"
    vmec_input.am = np.array(regime.pressure_shape, dtype=float)
    vmec_input.pres_scale = pressure_scale(case, regime)

    vmec_input.ncurr = 1
    vmec_input.pcurr_type = "power_series"
    vmec_input.ac = np.array(regime.current_shape, dtype=float)
    vmec_input.curtor = regime.current_fraction * characteristic_current(case)

    if free_boundary:
        vmec_input.lfreeb = True
        vmec_input.extcur = case.extcur
        vmec_input.nvacskip = level.nvacskip
    else:
        vmec_input.lfreeb = False
        vmec_input.extcur = np.array([])

    vmec_input.raxis_c = np.zeros(level.ntor + 1)
    vmec_input.zaxis_s = np.zeros(level.ntor + 1)
    vmec_input.raxis_c[0] = r_axis_guess
    vmec_input.rbc = rbc
    vmec_input.zbs = zbs
    vmec_input.return_outputs_even_if_not_converged = (
        return_outputs_even_if_not_converged
    )
    return vmec_input
