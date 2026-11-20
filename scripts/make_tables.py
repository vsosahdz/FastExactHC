#!/usr/bin/env python3
"""Rebuild the tables of the paper from the raw data in results/.

Usage:
    python3 scripts/make_tables.py [--results results] [--out tables.md]

Requires Python 3.8+, NumPy, and SciPy (only for the statistical tests).
Each function below rebuilds one table of the paper and returns it as Markdown.
A run that did not finish within the 10 s limit counts as 10 s, so it ranks last
in its block, as stated in the paper.
"""
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

LIMIT = 10.0
METHODS = {"wfg_c": "WFG", "qhv": "QHV", "qhv2": "QHV-II", "imatrix": "bitmatrix", "wfg_imatrix": "wfg-bitmatrix"}
FRONTS = {"concavo": "concave", "lineal": "linear", "convexo": "convex", "lineal_invertido": "inverted linear",
          "desconectado": "disconnected", "wfg1": "WFG1", "wfg2": "WFG2", "wfg3": "WFG3",
          "dtlz7_poblacion": "DTLZ7 (pop.)", "dtlz2_poblacion": "DTLZ2 pop.", "dtlz1_poblacion": "DTLZ1 (pop.)",
          "degenerado": "degenerate", "concavo_escalado": "scaled concave"}
QHV_MAX_D = 15  # maximum number of objectives supported by the original QHV code


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


def g3(x):
    """Three significant digits, as in the paper."""
    return f"{x:.3g}" if x >= 0.001 else f"{x:.2e}"


# ---------------------------------------------------------------- statistical study (results/05)

def load_statistics(results):
    """times[(front, n, d, seed)][method] = seconds, or LIMIT if the run did not finish."""
    times = defaultdict(dict)
    for r in read_csv(results / "05_estadistico_30_ejecuciones/statistics_results_30.csv"):
        method, d = r["method"], int(r["objectives"])
        if method == "qhv" and d > QHV_MAX_D:
            continue
        ok = r["status"] == "ok" and r["seconds"] not in ("", "nan")
        times[(r["instance"], int(r["points"]), d, int(r["seed"]))][method] = float(r["seconds"]) if ok else LIMIT
    return times


def finished(t):
    return t < LIMIT


def table_maxd(times):
    """Largest d solved within 10 s in at least half of the seeds (Table 'tab:maxd')."""
    solved = defaultdict(int)
    for (front, n, d, seed), by_method in times.items():
        for method, t in by_method.items():
            solved[(front, n, method, d)] += finished(t)
    rows = []
    for front in ("concavo", "lineal"):
        for n in (10, 20, 50, 100, 200, 500, 1000):
            row = [FRONTS[front], n]
            for method in METHODS:
                ds = [d for (f, nn, m, d), c in solved.items() if (f, nn, m) == (front, n, method) and c >= 15]
                row.append(max(ds) if ds else "-")
            rows.append(row)
    return table(["Front", "n"] + list(METHODS.values()), rows)


def table_median(times):
    """Median [interquartile range] on representative combinations of the concave front ('tab:mediana')."""
    combos = [(50, 10), (50, 15), (100, 8), (100, 12), (200, 6), (200, 10), (500, 5), (500, 8), (1000, 4), (1000, 7)]
    rows = []
    for n, d in combos:
        row = [n, d]
        for method in METHODS:
            v = [times[k][method] for k in times if k[:3] == ("concavo", n, d) and method in times[k]]
            if not v:
                row.append("n/a")
            elif sum(not finished(t) for t in v) >= len(v) / 2:
                row.append(">10")
            else:
                q1, med, q3 = np.percentile(v, [25, 50, 75])
                row.append(f"{g3(med)} [{g3(q1)}, {g3(q3)}]")
        rows.append(row)
    return table(["n", "d"] + list(METHODS.values()), rows)


def paired_blocks(times, other):
    return [(b["imatrix"], b[other]) for b in times.values() if "imatrix" in b and other in b]


def table_difficulty(times):
    """Other method's time / bitmatrix time by difficulty, blocks solved by both ('tab:dificultad')."""
    bins = [(0, 0.01, "<0.01"), (0.01, 0.1, "0.01-0.1"), (0.1, 1, "0.1-1"), (1, 10, "1-10")]
    rows = []
    for other in ("wfg_c", "qhv", "qhv2", "wfg_imatrix"):
        row = [METHODS[other]]
        pairs = [(a, b) for a, b in paired_blocks(times, other) if finished(a) and finished(b)]
        for lo, hi, _ in bins:
            sel = [(a, b) for a, b in pairs if lo <= max(a, b) < hi]
            if not sel:
                row.append("-")
                continue
            gm = math.exp(np.mean([math.log(b / a) for a, b in sel]))
            faster = 100 * np.mean([a < b for a, b in sel])
            row.append(f"{gm:.2f} ({faster:.0f}%)")
        rows.append(row)
    return table(["Method"] + [label for *_, label in bins], rows)


def holm(pvalues):
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    adjusted, running = [0.0] * len(pvalues), 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(pvalues) - rank) * pvalues[i]))
        adjusted[i] = running
    return adjusted


def table_wilcoxon(times):
    """Pairwise comparison of bitmatrix against each method over all blocks ('tab:wilcoxon')."""
    others, rows, pvals = ("wfg_c", "qhv", "qhv2", "wfg_imatrix"), [], []
    for other in others:
        pairs = [(a, b) for a, b in paired_blocks(times, other) if finished(a) or finished(b)]
        both = [(a, b) for a, b in pairs if finished(a) and finished(b)]
        faster = 100 * np.mean([a < b for a, b in pairs])
        speedup = math.exp(np.mean([math.log(b / a) for a, b in both]))
        only_bm = sum(finished(a) and not finished(b) for a, b in pairs)
        only_other = sum(finished(b) and not finished(a) for a, b in pairs)
        pvals.append(stats.wilcoxon([a for a, _ in pairs], [b for _, b in pairs]).pvalue)
        rows.append([METHODS[other], len(pairs), f"{faster:.1f}%", f"{speedup:.2f}x", only_bm, only_other])
    for row, p in zip(rows, holm(pvals)):
        row.append(f"{p:.1e}" if p >= 1e-10 else "<1e-10")
    four = ("wfg_c", "qhv2", "imatrix", "wfg_imatrix")
    blocks = [[b[m] for m in four] for b in times.values() if all(m in b for m in four) and any(finished(b[m]) for m in four)]
    chi2, p = stats.friedmanchisquare(*np.array(blocks).T)
    ranks = np.mean([stats.rankdata(row) for row in blocks], axis=0)
    footer = (f"\n\nFriedman (without QHV): {len(blocks)} blocks, chi2 = {chi2:.1f}, p = {p:.1e}. "
              "Average rank: " + ", ".join(f"{METHODS[m]} {r:.2f}" for m, r in zip(four, ranks)))
    return table(["Method", "blocks", "bitmatrix faster", "speed-up (geo. mean)", "only bitmatrix", "only other", "p (Holm)"], rows) + footer


def a12(x, y):
    """Vargha-Delaney A12: probability that a run of x is faster (smaller) than a run of y."""
    x, y = np.asarray(x), np.asarray(y)
    return float(np.mean((x[:, None] < y[None, :]) + 0.5 * (x[:, None] == y[None, :])))


def table_cells(times):
    """Per-combination Wilcoxon test (30 paired runs) and Vargha-Delaney effect size ('tab:celdas')."""
    combos = defaultdict(lambda: defaultdict(dict))
    for (front, n, d, seed), b in times.items():
        combos[(front, n, d)][seed] = b
    rows = []
    for other in ("wfg_c", "qhv", "qhv2", "wfg_imatrix"):
        counts = {"all": [0, 0, 0, []], "nontrivial": [0, 0, 0, []]}
        for runs in combos.values():
            pairs = [(b["imatrix"], b[other]) for b in runs.values() if "imatrix" in b and other in b]
            if not pairs or not any(finished(a) or finished(b) for a, b in pairs):
                continue
            x, y = [a for a, _ in pairs], [b for _, b in pairs]
            p = stats.wilcoxon(x, y).pvalue if any(a != b for a, b in pairs) else 1.0
            outcome = 1 if p >= 0.05 else (0 if np.median(x) < np.median(y) else 2)
            groups = ["all"] + (["nontrivial"] if max(np.median(x), np.median(y)) > 0.01 else [])
            for g in groups:
                counts[g][outcome] += 1
                counts[g][3].append(a12(x, y))
        row = [METHODS[other]]
        for g in ("all", "nontrivial"):
            plus, eq, minus, effects = counts[g]
            row += [plus, eq, minus, f"{np.median(effects):.2f}" if effects else "-"]
        rows.append(row)
    return table(["Method", "+", "=", "-", "A12", "+ (nontrivial)", "= (nontrivial)", "- (nontrivial)", "A12 (nontrivial)"], rows)


# ---------------------------------------------------------------- front geometries (results/01)

def table_times(results, d=10):
    """Running time (s) with d = 10 and 200 points ('tab:tiempos')."""
    rows01 = read_csv(results / "01_geometrias_13/fronts_comparison.csv")
    rows_c = read_csv(results / "01_geometrias_13/fronts_wfg_c.csv")
    value = {}
    for r in rows01 + rows_c:
        if int(r["objectives"]) == d:
            value[(r["instance"], r["method"])] = (float(r["seconds"]) if r["status"] == "ok" else None)
    columns = [("wfg_c", "WFG"), ("wfg", "WFG_py"), ("qhv", "QHV"), ("qhv2", "QHV-II"), ("imatrix", "bitmatrix"), ("wfg_imatrix", "wfg-bm.")]
    instances = ["lineal", "concavo", "convexo", "lineal_invertido", "desconectado", "wfg1", "wfg2", "dtlz7_poblacion"]
    rows = []
    for inst in instances:
        row = [FRONTS[inst]]
        for key, _ in columns:
            t = value.get((inst, key))
            row.append(">10" if t is None else (f"{t:.3f}" if t < 1 else f"{t:.2f}"))
        rows.append(row)
    return table(["Instance"] + [label for _, label in columns], rows)


# ---------------------------------------------------------------- largest solvable size (results/02)

def parse_size_table(path):
    """Rows 'instance d | n | front | method1 | ...' -> {(instance, d, method): [n solved]}."""
    solved, header = defaultdict(list), None
    for line in path.read_text(encoding="utf-8").splitlines():
        cells = [c.strip() for c in line.split("|")]
        if line.startswith("instancia"):
            header = cells[3:]
            continue
        m = re.match(r"(\S+)\s+(\d+)$", cells[0]) if header else None
        if not m:
            continue
        inst, d, n = m.group(1), int(m.group(2)), int(cells[1])
        for method, cell in zip(header, cells[3:]):
            if cell.endswith("s") and cell not in ("10+",):
                solved[(inst, d, method)].append(n)
    return solved


def table_size(results):
    """Largest number of points solved within 10 s ('tab:tamano')."""
    full = parse_size_table(results / "02_tamano_maximo/tamano_maximo_frentes_llenos.txt")
    wfg_c = parse_size_table(results / "02_tamano_maximo/tamano_maximo_wfg_original.txt")
    rows = []
    for d in (15, 20):
        for inst in ("lineal", "concavo", "convexo", "lineal_invertido", "desconectado", "wfg1"):
            row = [FRONTS[inst], d]
            for source, method in ((wfg_c, "wfg_c"), (full, "qhv2"), (full, "imatrix"), (full, "wfg_imatrix")):
                ns = source.get((inst, d, method), [])
                row.append(max(ns) if ns else "-")
            rows.append(row)
    return table(["Front", "d", "WFG", "QHV-II", "bitmatrix", "wfg-bm."], rows)


# ---------------------------------------------------------------- ablation (results/03)

def table_ablation(results):
    """Ablation: time (s) / tree nodes (thousands) at the hardest d solved by bitmatrix ('tab:ablacion')."""
    rows = read_csv(results / "03_sensibilidad_componentes/ablation_results.csv")
    data = {(r["instance"], int(r["objectives"]), r["variant"]): r for r in rows}
    variants = ["base", "+poda", "+filtro", "+orden", "+matriz=imatrix",
                "imatrix-poda", "imatrix-filtro", "imatrix-orden", "imatrix-matriz", "imatrix_ciclico"]
    labels = {"base": "base", "+poda": "+ pruning", "+filtro": "+ filtering", "+orden": "+ ordering",
              "+matriz=imatrix": "+ matrix (bitmatrix)", "imatrix-poda": "bitmatrix w/o pruning",
              "imatrix-filtro": "bitmatrix w/o filtering", "imatrix-orden": "bitmatrix w/o ordering",
              "imatrix-matriz": "bitmatrix w/o matrix", "imatrix_ciclico": "bitmatrix, cyclic ordering"}
    fronts = ["concavo", "lineal", "lineal_invertido", "desconectado", "wfg1", "dtlz2_poblacion"]
    hardest = {}
    for f in fronts:
        ds = [d for (inst, d, v), r in data.items() if inst == f and v == "+matriz=imatrix" and r["status"] == "ok"]
        hardest[f] = max(ds)
    out = []
    for v in variants:
        row = [labels[v]]
        for f in fronts:
            r = data.get((f, hardest[f], v))
            row.append(f"{float(r['seconds']):.3f} / {float(r['nodes']) / 1000:,.0f}" if r and r["status"] == "ok" else ">10")
        ratios = []
        for f in fronts:
            for d in (5, 8, 10, 12):
                ref, r = data.get((f, d, "+matriz=imatrix")), data.get((f, d, v))
                if ref and r and ref["status"] == r["status"] == "ok" and float(ref["seconds"]) > 0.01:
                    ratios.append(float(r["seconds"]) / float(ref["seconds"]))
        row.append(f"{math.exp(np.mean(np.log(ratios))):.2f}" if ratios else "-")
        out.append(row)
    header = ["Variant"] + [f"{FRONTS[f]}, d={hardest[f]}" for f in fronts] + ["rel. time"]
    return table(header, out)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", default=Path(__file__).resolve().parent.parent / "results", type=Path)
    parser.add_argument("--out", type=Path, help="write the tables to this Markdown file instead of stdout")
    args = parser.parse_args()
    times = load_statistics(args.results)
    sections = [
        ("Running time with d = 10 and 200 points (tab:tiempos)", table_times(args.results)),
        ("Largest number of points solved within 10 s (tab:tamano)", table_size(args.results)),
        ("Ablation study (tab:ablacion)", table_ablation(args.results)),
        ("Largest d solved in at least half of the seeds (tab:maxd)", table_maxd(times)),
        ("Pairwise comparison over all blocks (tab:wilcoxon)", table_wilcoxon(times)),
        ("Speed-up by difficulty (tab:dificultad)", table_difficulty(times)),
        ("Per-combination Wilcoxon test and effect size (tab:celdas)", table_cells(times)),
        ("Median [IQR] on representative combinations (tab:mediana)", table_median(times)),
    ]
    text = "\n\n".join(f"## {title}\n\n{body}" for title, body in sections) + "\n"
    if args.out:
        args.out.write_text(text, encoding="utf-8")
    else:
        print(text)


if __name__ == "__main__":
    main()
