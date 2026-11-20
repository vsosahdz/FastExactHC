# Reproducibility

This document explains how the experiments of the paper were run, what each file in
`results/` contains, which table of the paper it feeds, and how to rebuild those tables
from the raw data.

## What is included

| Included | Where |
|---|---|
| C code of the two methods (bitmatrix, wfg-bitmatrix) | `src/`, built with `make` |
| Correctness checks on the examples of the paper and on 14 fronts generated as in the experiments | `make test`, `make test-fronts` |
| Raw data of every experiment reported in the paper (one exception below) | `results/` |
| Script that rebuilds the tables of the paper from the raw data | `scripts/make_tables.py` |

Not included:
- The Python driver that generated the instances and timed all the methods.
- The code of the baselines. These are third-party implementations, obtained from their authors as described in the paper: WFG 1.02 in C, QHV, QHV-II in C++, and WFG from pygmo 2.19.8.
- The raw data of the table for the DTLZ2 population of 1000 points (*Up to 20 Objectives*). That table cannot be rebuilt from this repository.

## Environment

All timings were measured on an Apple M5 with 16 GB of memory, macOS 26.6, using Apple clang 21. The methods were compiled with `-O3 -march=native`.
- Instances: generated with pymoo 0.6.1.5 (DTLZ, WFG). The MaF-type geometries (inverted linear and scaled) were generated analytically.
- Reference point: `r = z_nad + 0.1 (z_nad - z_id)`.
- Time limit: 10 s per computation.

## Experiments and raw data

| Folder | Experiment | Paper table / figure |
|---|---|---|
| `results/01_geometrias_13/` | *Breadth*: 13 instances, 200 points, d in {3, 5, 8, 10, 12, 14, 16, 18, 20}, one seed. A method that exceeds the limit is not run in higher dimensions of that instance | Running time with d = 10 (`tab:tiempos`) and the count of fastest methods in *Running Time by Geometry* |
| `results/02_tamano_maximo/` | Largest number of points (10, 20, 30, 50, 75, 100, 150) solved within 10 s on full fronts, d = 15 and 20 | Largest solvable size (`tab:tamano`) |
| `results/03_sensibilidad_componentes/` | Ablation: the same function with components switched on or off, 200 points. It uses its own sample of points, independent of `01` | Ablation (`tab:ablacion`) |
| `results/04_estadistico_10_ejecuciones/` | Earlier statistical study with 10 runs per combination | Mentioned in the text; its largest solved d matches the 30-run study |
| `results/05_estadistico_30_ejecuciones/` | *Statistical study*: concave and linear fronts, n in {10, 20, 50, 100, 200, 500, 1000}, d = 2 to 20, 30 seeds (7980 blocks) | `tab:maxd`, `tab:wilcoxon`, `tab:dificultad`, `tab:celdas`, `tab:mediana`, and the figures of the statistical study |

### Names used in the raw data

| In `results/` | In the paper |
|---|---|
| `imatrix` | bitmatrix |
| `wfg_imatrix` | wfg-bitmatrix |
| `wfg_c` | WFG (original C code, version 1.02) |
| `wfg` (folders 01 and 02) | WFG from pygmo (WFG_py) |
| `qhv`, `qhv2` | QHV, QHV-II |
| `qmc`, `qmc_depth` (folder 01) | quasi-Monte Carlo approximations; not reported in the paper |
| `concavo`, `lineal`, `convexo`, `lineal_invertido`, `desconectado`, `degenerado`, `concavo_escalado` | concave, linear, convex, inverted linear, disconnected, degenerate, scaled concave |
| `*_poblacion` | population instances (random decision vectors of DTLZ) |

### Status codes in `statistics_results_30.csv`

| `status` | Meaning |
|---|---|
| `ok` | finished within 10 s; `seconds` and `hv` are valid |
| `timeout` | exceeded 10 s |
| `timeout_skipped` | not run: the first seed of the same combination exceeded the limit, so the remaining seeds count as failures |
| `dominated_skip` | not run: the method failed on an easier combination (fewer points and fewer objectives), so it counts as failed here |
| `n/a` or empty (`qhv`) | not applicable: the original QHV code supports at most 15 objectives |

Any status other than `ok` counts as a run that did not finish. In rank-based comparisons it ranks last in its block (time = 10 s).

## Rebuilding the tables

```bash
pip install -r scripts/requirements.txt    # NumPy and SciPy
python3 scripts/make_tables.py             # prints the tables as Markdown
python3 scripts/make_tables.py --out tables.md
```

The script rebuilds, from the raw data and with the definitions given in the paper:
- `tab:tiempos`, `tab:tamano`, `tab:ablacion`;
- `tab:maxd`, the largest d solved in at least half of the seeds;
- `tab:wilcoxon`: Wilcoxon signed-rank test over all blocks with Holm correction, and the Friedman test without QHV;
- `tab:dificultad`, the geometric mean of the time ratio by difficulty;
- `tab:celdas`: per-combination Wilcoxon test on the 30 paired runs, α = 0.05, and the Vargha–Delaney Â12;
- `tab:mediana`, the median and interquartile range.

On the data in this repository, every value matches the paper except one rounding: convex, WFG_py, d = 10 is 2.2148 s.

## Checking the implementation

```bash
make test          # the 2D and 3D examples of the paper (hypervolume 6 and 14.125)
make test-fronts   # both methods on 14 fronts; relative tolerance 1e-9 against the recorded hypervolume
```

## Notes on the timings

- `./bitmatrix -t` measures the whole call with `clock_gettime`, including the construction of the comparison matrix W. The code is single-threaded.
- For WFG, QHV, and QHV-II, the paper uses the time reported by the programs themselves, excluding data reading and process startup.
- For d ≤ 12, wfg-bitmatrix calls bitmatrix directly. Their difference therefore measures the timing noise of the setup, which the paper reports (2.6% in a typical combination, up to about 10%).
