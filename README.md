# FastExactHC

Exact hypervolume computation by pivot partitioning with an inherited comparison
matrix and pruned inclusion–exclusion. This repository contains the C code of the
two exact methods of the paper and the raw data of the experiments.

- **bitmatrix**: pivot partitioning (as in QHV-II) with three exact improvements:
  a pruning rule for the inclusion–exclusion base case, a bitwise comparison
  matrix `W` that is built once and inherited throughout the recursion, and an
  ordering of the subregions by the number of points they receive.
- **wfg-bitmatrix**: WFG recursion for high dimension, with limit sets filtered by
  the same matrix `W`, switching to bitmatrix when 12 or fewer objectives remain.

All objectives are minimized. The hypervolume is computed with respect to a
reference point `r`, and only points strictly below `r` contribute.

## Requirements

A C99 compiler (`gcc` or `clang`) and `make`. Tested on macOS (Apple clang) and
valid for Linux. No external libraries are needed.

## Build and test

```bash
make          # builds the program ./bitmatrix and the shared library
make test     # checks both methods on the examples of the paper
make run      # runs the examples and prints the computation time
make clean    # removes the build products
```

`make` also builds `libbitmatrix.so` (Linux) or `libbitmatrix.dylib` (macOS).
To build without `-march=native`, use `make CFLAGS=-O3`.

## Usage

```bash
./bitmatrix [-m bitmatrix|wfg-bitmatrix] [-r r1,...,rd] [-k leaf] [-s switch] [-t] file
```

| Option | Meaning | Default |
|---|---|---|
| `-m` | method | `bitmatrix` |
| `-r` | reference point, comma separated | maximum of each objective + 1 |
| `-k` | maximum leaf size `k_max` | 10 |
| `-s` | number of objectives at or below which wfg-bitmatrix switches to bitmatrix | 12 |
| `-t` | also print the computation time | off |

The input file has one point per line, with the objective values separated by
spaces, tabs, or commas. Empty lines and lines starting with `#` are ignored.
The program prints the hypervolume on the first line.

```bash
$ ./bitmatrix -r 4,4,4 examples/example_3d.txt
14.125
$ ./bitmatrix -t -m wfg-bitmatrix -r 4,4,4 examples/example_3d.txt
14.125
time 0.000004 s (n=5, d=3, wfg-bitmatrix)
```

For fronts with many objectives (about 13 or more), `wfg-bitmatrix` is usually
the faster option; below that, both methods run the same bitmatrix engine.

## Using the library

The functions take the points as a row-major array (`points[i * d + o]` is
objective `o` of point `i`):

```c
double hypervolume_c_imatrix(const double *points, size_t n, size_t d,
                             const double *reference, size_t exact_limit);
double hypervolume_c_wfg_imatrix(const double *points, size_t n, size_t d,
                                 const double *reference, size_t exact_limit,
                                 size_t switch_dimension);
```

`exact_limit = 0` and `switch_dimension = 0` select the defaults (10 and 12).
The functions return `-1.0` only if memory allocation fails. From Python, the
shared library can be loaded with `ctypes`.

## Repository layout

| Path | Content |
|---|---|
| `src/hypervolume_c.c` | base pivot partitioning and shared helpers |
| `src/hypervolume_c_ie_pruned.c` | pruned inclusion–exclusion in the leaves |
| `src/hypervolume_c_imatrix.c` | bitmatrix (`hypervolume_c_imatrix`) |
| `src/hypervolume_c_wfg_imatrix.c` | wfg-bitmatrix (`hypervolume_c_wfg_imatrix`) |
| `src/main.c` | command-line program |
| `examples/` | the 2D and 3D examples of the paper (hypervolume 6 and 14.125 with `r = (4,4)` and `r = (4,4,4)`) |
| `results/` | raw data of the experiments of the paper |

Each source file includes the previous one, so a single translation unit
(`src/main.c` or `src/hypervolume_c_wfg_imatrix.c`) builds everything.

## Experimental results

The folder `results/` contains the raw data, without processing, of all the runs
reported in the paper; `results/README.md` (in Spanish) describes each experiment,
the environment, and the command used.

| Folder | Experiment |
|---|---|
| `01_geometrias_13` | 13 front geometries, 200 points, up to 20 objectives |
| `02_tamano_maximo` | largest front solved in less than 10 s with 15 and 20 objectives |
| `03_sensibilidad_componentes` | component sensitivity (ablation) study |
| `04_estadistico_10_ejecuciones` | statistical study, 10 runs per combination |
| `05_estadistico_30_ejecuciones` | statistical study, 30 runs per combination (main results) |

