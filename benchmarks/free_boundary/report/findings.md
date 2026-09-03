# VMEC++ free-boundary bench: findings

A reproducible benchmark of the VMEC++ free-boundary solver on 312 QUASR
coil sets, built so that a change to the solver can be judged on numbers.
`report.md` next to this file is the machine-written report (tables and
figures generated from the ledgers in `../results/`); this file is the reading
of it. The code is `benchmarks/free_boundary/` (see its README). Solver:
VMEC++ built from commit 91ab6b3 of main (2 Sep 2026, package version
0.7.4.dev33+g91ab6b351); SIMSOPT 1.11.1.

```
QUASR index (371k) --stratified sample--> 312 cases (12 already in tests/data/quasr)
        |
   coils --> mgrid table --> vmecpp.run (free boundary)  <-- QUASR surface as initial guess
        |                          |
   exact Biot-Savart field         wout + solver log
        |                          |
        +---- accuracy metrics ----+---- outcome class ----> ledger --> report
```

Three axes are swept: the physics regime (vacuum; three finite-beta regimes;
a fixed-boundary control), the VMEC++ level (the resolution and flow-control
settings of `tests/test_free_boundary_quasr.py` as the base, plus Fourier,
radial and tolerance ladders and two flow-control variants) and the fidelity
of the coil field handed to the solver. Every run ends in one of twelve
outcome classes; every converged vacuum run is measured against the exact
coil field, which is the ground truth in vacuum because the plasma carries
no current.

## Headline findings

1. **At the settings of the existing test module, 57 percent of QUASR
   configurations converge in vacuum (178 of 312).** The failures split into
   two kinds. 53 (17 percent) fail the initial Jacobian check before the
   first iteration; the same 53 fail identically with the boundary held fixed,
   so this is the initial boundary, not the free-boundary coupling. The other
   81 (26 percent) fail only when the boundary is free: 49 by repeated
   Jacobian resets, 31 by running out of iterations (14 on a plateau, 13
   oscillating, 4 still converging), 1 by leaving the mgrid box. 77 of the 81
   converge with the boundary held fixed; the remaining four (all repeated
   Jacobian resets) fail both ways.

2. **The converged vacuum boundaries are flux surfaces of the coil field to
   about 0.2 percent, not to numerical precision.** B.n/|B| on the VMEC++
   boundary has a median of 1.7e-3 (p90 5.8e-3, worst 2.1e-2). Field lines of the coil
   field started on the boundary stray from it by a median of 0.035 minor
   radii over 50 turns (p90 0.13). The mgrid table VMEC++ was given is accurate
   to 7.6e-5 (median), 25 times better than the solution, so this is the
   solver's error, not the input's.

3. **The error falls slowly with Fourier resolution and hardly at all with
   radial resolution.** On the cases that converge at every level, B.n/|B|
   goes 7.7e-4 to 5.8e-4 to 4.0e-4 to 3.3e-4 for mpol = ntor = 6, 8, 10, 12;
   raising ns from 31 to 101 at mpol 8 moves it from 7.7e-4 to 5.7e-4 with
   nothing gained past ns = 51; tightening ftol from 1e-8 to 1e-13 moves it
   from 6.6e-4 to 5.0e-4 at five times the iterations (the field-line
   deviation responds more, 0.032 to 0.0085 minor radii). Nothing in the sweep
   approaches the 1e-12 that the product owner asked about. Raising the
   resolution also changes which cases converge: 10 of the 20 ladder cases
   converge at the base level, 8 to 9 at mpol 8 to 12, 6 at ftol = 1e-13.

4. **The existing test's coils file caps its own accuracy.** With the coil
   curves as QUASR stores them (160 points per coil) the mgrid table is wrong
   by 2.5e-4 (median over the ladder subset), within a factor of three of the
   solver error at that resolution; with 1280-point polygons the floor drops
   to 4e-5 and with 2560 points on a 201 x 201 grid to 1e-5, while the solver
   error stays at 6e-4 to 9e-4. The bench uses 1280-point polygons by default.

5. **The two finite-beta regimes sketched in the test module are above the
   equilibrium beta limit for most low-iota cases.** The Shafranov shift
   scales as beta / (2 eps iota^2); at 1 percent beta none of the 66 cases
   with iota below 0.25 converge, and 16 runs leave the mgrid box. A regime
   that scales beta per case to a nominal shift of 0.2 minor radii converges
   in 96 cases and produces exactly that shift (median 0.196 a, outboard in
   93 of 93), with the target beta met to 3 percent. The 2 percent plus net
   current regime converges in 9 cases.

6. **Free boundary costs about twice the iterations of fixed boundary and
   about 2.7 times the wall time.** Converged vacuum runs take a median of
   1370 force iterations (p90 2660) against 761 fixed-boundary; the vacuum
   pressure is switched on at iteration 64 (median). Neither flow-control
   variant changes much: `nvacskip = 1` converges the same 10 of 20 ladder
   cases in the same number of iterations at a higher cost per iteration
   (NESTOR every iteration), `delt = 0.5` converges 11 of 20.

7. **Spectral condensation of the initial boundary recovers part of the
   initial-Jacobian class.** Reparametrising the poloidal angle of the QUASR
   surface with SIMSOPT's `condense_spectrum` (shape kept to 1e-3 minor
   radii) before the run turns 10 of the 53 initial-Jacobian cases into
   converged free-boundary solutions (14 of 53 with the boundary held fixed);
   39 still fail the initial Jacobian. The constrained condensation only
   reduced the spectral width by a median 10 percent; on the 21 cases where
   it achieved more than 20 percent, 10 converged. Raising the Fourier
   resolution to mpol 8 or 10 recovers none of the 53. The lever is the
   parametrisation of the input, and a less constrained condensation is the
   obvious next test.

8. **Where traced field lines are lost, it is the configuration, not the
   solver, in 26 of 28 cases.** 28 of the 178 converged vacuum cases lose at
   least one of eight traced lines (19 lose half or more). In 26 of the 28
   the same tracing started on the QUASR surface itself loses lines too: the
   coil field's last closed flux surface lies inside the QUASR boundary
   (an optimisation residual of the database, concentrated at iota above
   1.5 and aspect ratio 12 and above). Field-line deviation is therefore
   reported only over lines that stay confined, and the two remaining cases
   (QUASR 258805 and 1094744) are the ones where the VMEC++ boundary sits
   outside a confined region that the QUASR surface is inside.

## Robustness

`report.md`, section Robustness. The outcome classes are the solver's own
terminal states plus a split of budget exhaustion by what the residual was
doing (near ftol, still falling, oscillating, on a plateau).

Convergence is not strongly organised by coil proximity: 47 percent of cases
with coils closer than 1.2 minor radii converge, 63 percent between 1.2 and
1.9, 58 percent beyond. It is organised by field-period number (nfp = 1: 35
percent; nfp 2 to 8: 44 to 80 percent) and by the initial boundary. The
initial-Jacobian failures have a median spectral width (SIMSOPT definition)
of 39 against 11.5 for the converged cases; no case with a spectral width
below 6 failed the initial Jacobian.

The fixed-boundary control separates the classes cleanly:

```
free-boundary class        cases   fixed-boundary twin
converged                  178     176 converged, 2 fail
bad_jacobian                53     53 bad_jacobian          -> initial boundary
jacobian_75_times           49     45 converged, 4 fail     -> mostly free-boundary coupling
max_iter (any)              31     31 converged             -> free-boundary coupling
vacuum_grid_exceeded         1      1 converged
```

The residual histories (figure `fig_residuals_vacuum.png`) show the three
free-boundary failure modes: Jacobian resets that never settle, a plateau of
the residual at 1e-2 to 1e-4 after the vacuum pressure is switched on, and a
comb-like oscillation at the same level.

## Accuracy

`report.md`, section Accuracy, and figures `fig_accuracy.png`,
`fig_examples.png`, `fig_ladders.png`, `fig_field_ladder.png`.

The metrics are local (B.n/|B| on the boundary), integrated (field-line
deviation after 50 turns) and geometric (distance to the QUASR surface the
run started from). They agree with each other: the deviation of traced lines
is the integrated normal-field error over the turns, the traced rotational
transform agrees with the VMEC++ edge iota to 0.8 percent (median), and the
VMEC++ boundary sits a median 0.022 minor radii from the QUASR surface, which
is itself an optimised approximation to a flux surface.

The floor of the mgrid table is measured for every run, evaluated the way
VMEC++ uses the table (on the toroidal planes, bilinear in R and Z), on the
QUASR surface the run starts from; the converged boundary lies within 0.02
minor radii of that surface (median), so the floor at the final boundary is
the same to within the smoothness of the table. Every accuracy figure carries
it as a dashed line. The solver error is above the floor by a factor of 25
(median) and does not move when the floor is lowered by a further factor of
25 on the same cases (`fig_field_ladder.png`), which is the strict version of
the test: the plateau is the solver's.

## Speed

`report.md`, section Speed, and figures `fig_speed.png`,
`fig_flow_control.png`. Iteration counts are portable; wall times are laptop
timings (Apple Silicon, 10 cores, three threads per run, three runs at a
time) and only comparable within one ledger.

## Finite beta

`report.md`, section Finite beta. There is no exact ground truth at
finite beta (the plasma current contributes to the field outside), so the
bench reports convergence, the achieved beta, the outboard shift of the axis
against the vacuum solution of the same case, and for the current-carrying
regime the net current (achieved to 1 percent of the prescribed value in the
nine converged cases). The scaled regime hits its 4 percent cap for the
high-iota cases, and those are the ones that fail: 65 percent of cases with
iota below 0.25 converge at scaled beta but 13 percent of those above 0.75.

## What the bench does not claim

* The 57 percent is a score on a stratified stress suite that spans the whole
  QUASR index (field periods, coil distance, aspect ratio, iota up to 4.5),
  not a failure rate on a typical workload. The twelve configurations already
  in `tests/data/quasr` converge in 5 of 12.
* The finite-beta regimes measure robustness and sanity, not correctness; the
  2 percent plus current regime, with a median axis shift of 0.84 minor radii
  in the nine converged cases, is a stress result.
* No accuracy statement at finite beta. That needs a virtual-casing
  evaluation or an analytic equilibrium; the runner has the hooks.
* No comparison with VMEC2000, GVEC or DESC. The bench measures VMEC++
  against the coil field, which any of them can be measured against with the
  same code.
* No tuning. The two flow-control variants are data points.
* "Not converged" means the tolerance was not reached within the budget at
  the base settings, nothing more; the classes say what the run was doing.
* One machine, one solver build, one manifest seed. The manifest and the
  ledgers are checked in so the numbers can be reproduced or diffed.

## How to use it

From the repository root, with VMEC++ and SIMSOPT installed (the shipped
ledgers under `benchmarks/free_boundary/results/` are Git LFS files):

```bash
python -m benchmarks.free_boundary fetch
python -m benchmarks.free_boundary run --out results/after --levels base \
    --fields f1280 --regimes vacuum,vacuum_fixed
python -m benchmarks.free_boundary report --results results/after --out report/after
python -m benchmarks.free_boundary compare \
    --baseline benchmarks/free_boundary/results/main \
    --candidate results/after --out report/compare
```

`compare` lists the status transitions between two ledgers of the same
matrix, the newly failing and newly converging cases, and the ratio of
iterations, wall time, B.n and field-line deviation on the runs that converged
in both. That is the validation the product owner described for flow-control
changes ("converges significantly faster on average and doesn't stop
converging on any equilibria it previously worked on"), made mechanical.
