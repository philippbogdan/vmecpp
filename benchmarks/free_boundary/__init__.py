# SPDX-FileCopyrightText: 2024-present Proxima Fusion GmbH <info@proximafusion.com>
#
# SPDX-License-Identifier: MIT
"""Free-boundary test bench for VMEC++.

A reproducible benchmark of the free-boundary solver on QUASR coil sets: which
configurations converge, how the others fail, how accurate the converged
solutions are against the exact coil field, and how long each solve takes.

See ``README.md`` in this directory for the design and ``python -m
benchmarks.free_boundary --help`` for the command line.
"""
