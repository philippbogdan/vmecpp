# Free-boundary test bench

A reproducible benchmark of the VMEC++ free-boundary solver on QUASR coil sets,
so that a change to the solver can be judged on numbers: what converges, how the
rest fails, how accurate the converged solutions are against the exact coil
field, and how long each solve takes.

```
QUASR index (371k configs) --select--> manifest.csv (312 cases, stratified)
                                            |
                serial<ID>.json (surfaces + coils, one file per case)
                                            |
        +-----------------------------------+-----------------------------------+
        |                                   |                                   |
   coils file (polygons)           outermost surface                    Biot-Savart on the
        |                          = initial guess + phiedge            SIMSOPT coil curves
   mgrid response table                     |                           = ground truth
   (VMEC++ makegrid)  ------------>  vmecpp.run(lfreeb)  <---------------------+
                                            |                                   |
                             wout + solver log + timings                        |
                                            |                                   |
                         classify (taxonomy)   metrics (B.n, field lines, |B|^2)
                                            |
                                     results.jsonl  --report-->  tables + plots
```

## The three axes

* **Regime** (`settings.REGIMES`): `vacuum` (zero pressure, zero current; the
  only regime with an exact ground truth), `vacuum_fixed` (the same case with
  the QUASR surface imposed as a fixed boundary: a control that separates
  failures of the free-boundary coupling from failures of the initial
  boundary), `beta_scaled` (beta scaled per case to a nominal Shafranov shift
  of 0.2 minor radii, `beta = 0.2 * 2 * eps * iota^2`, capped at 4 percent),
  and the two regimes sketched in `tests/test_free_boundary_quasr.py`:
  `beta1` (1 percent) and `beta2_current` (2 percent plus a net toroidal
  current of 2 percent of `2 pi R B / mu0`). The fixed-beta regimes sit above
  the equilibrium beta limit (about `eps * iota^2`) for low-iota configurations
  and are kept as a stress test.
* **Level** (`settings.LEVELS`): resolution and flow control. `base` is the
  setting of the existing test module; `m6..m12` is a Fourier ladder at a fixed
  radial sequence, `ns31..ns101` a radial ladder at `mpol = ntor = 8`,
  `ftol8..ftol13` a tolerance ladder, `nvac1` and `delt05` flow-control
  variants, `base_cond` the base setting with the initial boundary's poloidal
  angle reparametrised by SIMSOPT's `condense_spectrum` (shape kept to 1e-3
  minor radii) before truncation.
* **Field** (`settings.FIELDS`): how faithfully the coil field reaches the
  solver. `f160` uses the coil curves as QUASR stores them (160 points) on a
  101 x 101 grid, `f1280` (the default) 1280-point polygons, `f2560` 2560-point
  polygons on 201 x 201, `f1280w` a wide box for finite-beta runs. The error of
  the table on the boundary is measured for every run (`mgrid_floor_rms`) and is
  the floor below which no solver result can be trusted.

## Case set

`cases/manifest.csv` is a stratified sample of the QUASR index: the product of
the number of field periods, terciles of the coil-to-surface distance in minor
radii and terciles of the aspect ratio, drawn round-robin with a fixed seed, plus
the twelve configurations checked into `tests/data/quasr` (subset `reference`).
Twenty cases carry `ladder = True` and are used for the resolution ladders.
The manifest carries the QUASR index values (iota profile, elongation, coil
distances, quasi-symmetry error) so results can be sliced by them.

## Outcome taxonomy

Every run ends in exactly one class (`classify.CLASSES`): `converged`,
`vacuum_grid_exceeded`, `jacobian_75_times`, `bad_jacobian`, `early_failure`,
`nan_residual`, `max_iter_near` (within a decade of ftol), `max_iter_slow` (still
converging, would reach ftol within three budgets), `max_iter_oscillating`,
`max_iter_stalled`, `timeout`, `error_other`. The split of the budget-exhausted
runs uses the force-residual history of the last stage.

## Accuracy metrics (vacuum, converged)

* `bn_rms`, `bn_max`: area-weighted RMS and maximum of the pointwise ratio
  B . n / |B| of the exact coil field on the VMEC++ boundary, on a 128 x 48
  (per period) grid, which resolves the residual at the boundary's truncation
  modes. Independent of the mgrid table. An independent SIMSOPT evaluation
  reproduces the ledger values to four digits.
* `fl_dev_rms`, `fl_dev_max`, `fl_lost_fraction`, `fl_iota`: field lines of the
  coil field traced from the boundary; distance of their phi = 0 crossings to
  the boundary in minor radii (lines that stray by more than half a minor
  radius are counted as lost and excluded); the traced rotational transform
  from the winding about the section centroid.
* `fl_lost_quasr`: the same tracing started on the QUASR surface itself. A
  configuration whose own surface loses lines has its last closed flux surface
  inside the QUASR boundary, so lost lines from the VMEC++ boundary of that
  case say nothing about the solver.
* `b2_rms`, `b2_max`: |B|^2 on the plasma side (extrapolated from the two
  outermost half-grid surfaces) against the coil field. Secondary: the
  extrapolation adds an error of order (1/ns)^2.
* `quasr_dist_rms`, `quasr_dist_max`, `volume_ratio_quasr`: geometric distance
  to and volume ratio against the QUASR surface the run started from.
* `mgrid_floor_rms`, `mgrid_floor_max`: error of the response table on the
  boundary, evaluated as VMEC++ uses it (on the toroidal planes, bilinear in R
  and Z).
* `exact_field_check`: the Biot-Savart reference against a re-quadratured
  rebuild of the coil set from its base coils; validates the ground truth.
  When the two differ by more than 1e-8 the rebuilt (640-point) coils are
  used as the reference (`reference_coil_points` records which).
* Convergence is decided from `ier_flag`, the final `ns` and the three force
  residuals against `ftol`, never from the `reason` string of the wout (which
  reads "convergence was not reached" for code 0, converged or not).

## Running

From the repository root, with VMEC++ and SIMSOPT installed:

```bash
python -m benchmarks.free_boundary select --n-main 300      # rebuild the manifest
python -m benchmarks.free_boundary fetch                    # download serial files
python -m benchmarks.free_boundary run --out results/main --levels base \
    --fields f1280 --regimes vacuum --workers 3 --threads 3
python -m benchmarks.free_boundary run --out results/main --levels base \
    --fields f1280w --regimes beta_scaled,beta1,beta2_current
python -m benchmarks.free_boundary run --out results/ladder --subset ladder \
    --levels m6,m8,m10,m12,ns31,ns71,ns101 --fields f1280 --regimes vacuum
python -m benchmarks.free_boundary report --results results/main results/ladder \
    --out report
python -m benchmarks.free_boundary compare --baseline results/main \
    --candidate results/main_after_change --out report/compare
python -m benchmarks.free_boundary retrace --results results/main   # field-line metrics only
```

`run --redo error_other,timeout` drops the records of the named classes and
runs them again. `retrace` recomputes the field-line metrics of a ledger from
the stored boundaries without re-running the solver.

Serial files and the QUASR index are cached under
`~/.cache/vmecpp_free_boundary_bench` (override with `VMECPP_FB_BENCH_CACHE`).
A `run` is resumable: results already in `results.jsonl` are skipped. Do not
point two concurrent `run` commands at the same output directory. Each
(case, level, field) job is a subprocess with a wall-clock limit; the solver's
own log is captured per run and its markers (vacuum activation iteration,
convergence-problem resets, grid overruns, per-stage iterations) are kept in
the ledger.

## Reading the numbers

Every accuracy number should be read next to `mgrid_floor_rms` of the same run:
a solver cannot beat the field it is given. The `f160` field spec reproduces the
setting of the existing test module and has a floor of a few 1e-4; `f1280`
brings it to a few 1e-5 and `f2560` below 1e-5.
