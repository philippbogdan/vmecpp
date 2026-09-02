# SPDX-FileCopyrightText: 2024-present Proxima Fusion GmbH <info@proximafusion.com>
#
# SPDX-License-Identifier: MIT
"""Presets for the free-boundary bench.

Three independent axes are swept:

* ``Regime``: the plasma physics imposed on top of a QUASR coil set (vacuum,
  finite beta without net current, finite beta with net toroidal current).
* ``Level``: the VMEC++ discretisation and flow-control settings.
* ``FieldSpec``: how faithfully the coil field is represented in the mgrid
  response table that VMEC++ actually sees (coil polygon resolution and
  R-Z grid resolution). This is a separate axis on purpose: the error of the
  vacuum field handed to the solver bounds the accuracy any solver can reach.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class Regime:
    """A physics regime to impose on a QUASR configuration."""

    name: str
    target_beta: float
    # Net toroidal current as a fraction of I_char = 2*pi*<R>*<B_phi>/mu0.
    current_fraction: float = 0.0
    # power_series coefficients of the (unscaled) pressure shape p(s) ~ (1-s).
    pressure_shape: tuple[float, ...] = (1.0, -1.0)
    # power_series coefficients of the toroidal current density shape.
    current_shape: tuple[float, ...] = (1.0, -1.0)
    # If set, target_beta is ignored and beta is scaled per case to a nominal
    # Shafranov shift of this many minor radii: beta = shift * 2 * eps * iota^2
    # with eps = 1 / aspect ratio and iota the QUASR edge rotational transform,
    # capped at beta_cap. Fixed-beta regimes ignore the equilibrium limit
    # (which scales as eps * iota^2), so a 1 percent target is far above it
    # for low-iota configurations.
    shift_target: float | None = None
    beta_cap: float = 0.04
    # False runs the same case fixed-boundary (the QUASR surface imposed): a
    # control that separates failures of the free-boundary coupling from
    # failures of the boundary representation or the initial axis guess.
    free_boundary: bool = True
    description: str = ""


REGIMES: dict[str, Regime] = {
    "vacuum": Regime(
        "vacuum", 0.0, description="zero pressure, zero net toroidal current"
    ),
    "beta1": Regime(
        "beta1", 0.01, description="about 1 percent volume-averaged beta, no current"
    ),
    "beta2_current": Regime(
        "beta2_current",
        0.02,
        current_fraction=0.02,
        description="about 2 percent beta plus a net toroidal current",
    ),
    "vacuum_fixed": Regime(
        "vacuum_fixed",
        0.0,
        free_boundary=False,
        description="vacuum, fixed boundary on the truncated QUASR surface (control)",
    ),
    "beta_scaled": Regime(
        "beta_scaled",
        0.0,
        shift_target=0.2,
        description="beta scaled per case to a nominal Shafranov shift of 0.2 a",
    ),
}


@dataclass(frozen=True)
class Level:
    """VMEC++ resolution and flow-control settings for one run."""

    name: str
    mpol: int
    ntor: int
    ns_array: tuple[int, ...]
    nzeta: int
    ftol: float = 1e-9
    niter: int = 2000
    nvacskip: int = 6
    delt: float = 1.0
    tcon0: float = 1.0
    description: str = ""

    @property
    def ns_final(self) -> int:
        return self.ns_array[-1]


LEVELS: dict[str, Level] = {
    # The settings of tests/test_free_boundary_quasr.py: the reference point.
    "base": Level(
        "base", 6, 6, (8, 16, 31), 24, description="tests/test_free_boundary_quasr.py"
    ),
    # Fourier ladder at a fixed radial sequence.
    "m6": Level("m6", 6, 6, (8, 16, 31, 51), 24, description="Fourier ladder"),
    "m8": Level("m8", 8, 8, (8, 16, 31, 51), 32, description="Fourier ladder"),
    "m10": Level("m10", 10, 10, (8, 16, 31, 51), 40, description="Fourier ladder"),
    "m12": Level("m12", 12, 12, (8, 16, 31, 51), 48, description="Fourier ladder"),
    # Radial ladder at fixed mpol = ntor = 8 (ns = 51 is shared with m8).
    "ns31": Level("ns31", 8, 8, (8, 16, 31), 32, description="radial ladder"),
    "ns71": Level("ns71", 8, 8, (8, 16, 31, 51, 71), 32, description="radial ladder"),
    "ns101": Level(
        "ns101", 8, 8, (8, 16, 31, 51, 71, 101), 32, description="radial ladder"
    ),
    # Tolerance ladder at the m8 resolution.
    "ftol8": Level("ftol8", 8, 8, (8, 16, 31, 51), 32, ftol=1e-8),
    "ftol11": Level("ftol11", 8, 8, (8, 16, 31, 51), 32, ftol=1e-11, niter=4000),
    "ftol13": Level("ftol13", 8, 8, (8, 16, 31, 51), 32, ftol=1e-13, niter=8000),
    # Flow-control variants of the base level.
    "nvac1": Level("nvac1", 6, 6, (8, 16, 31), 24, nvacskip=1),
    "delt05": Level("delt05", 6, 6, (8, 16, 31), 24, delt=0.5),
}


@dataclass(frozen=True)
class FieldSpec:
    """How the coil field is turned into the mgrid response table."""

    name: str
    # R and Z grid points of the mgrid box (the toroidal count equals nzeta).
    points_rz: int
    # Fraction of the boundary extent added on each side of the mgrid box.
    margin: float
    # Straight-segment polygon points per coil written to the coils file.
    coil_points: int
    description: str = ""


FIELDS: dict[str, FieldSpec] = {
    "f160": FieldSpec(
        "f160", 101, 0.4, 160, description="coil curves as stored by QUASR"
    ),
    "f1280": FieldSpec("f1280", 101, 0.4, 1280, description="bench default"),
    "f1280w": FieldSpec(
        "f1280w", 151, 1.0, 1280, description="wide box for finite-beta regimes"
    ),
    "f2560": FieldSpec("f2560", 201, 0.4, 2560, description="high fidelity"),
}


def presets_as_dict() -> dict:
    """All presets, JSON-serialisable, for the run metadata."""
    return {
        "regimes": {k: asdict(v) for k, v in REGIMES.items()},
        "levels": {k: asdict(v) for k, v in LEVELS.items()},
        "fields": {k: asdict(v) for k, v in FIELDS.items()},
    }
