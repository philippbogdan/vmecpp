# Shipped ledgers

One directory per run of the bench matrix, each with the ledger
(`results.jsonl.gz`, one JSON record per run, stored with Git LFS) and the run
metadata (`meta.json`: settings, presets, environment).

| ledger | what | records |
|---|---|---|
| `main` | all 312 cases at the base level: vacuum, fixed-boundary control, three finite-beta regimes | 1560 |
| `ladder` | the 20 ladder cases in vacuum at the Fourier, radial, tolerance and flow-control levels and at three coil-field fidelities | 280 |
| `badjac` | the 53 cases that failed the initial Jacobian, rerun at mpol = ntor = 8 and 10 | 106 |
| `condensed` | the same 53 cases with a spectrally condensed initial boundary, free and fixed boundary | 106 |

`python -m benchmarks.free_boundary report --results benchmarks/free_boundary/results/main ...`
reads the compressed ledgers directly; `compare` diffs a new run against one
of them. The vacuum field-line metrics in `main` were traced over 50 toroidal
turns, those in `ladder` over 100 (the `fl_turns` field of each record).
