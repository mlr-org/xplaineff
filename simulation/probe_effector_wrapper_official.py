#!/usr/bin/env python3
"""Compare effector's official Regional*.fit() runtime with the benchmark wrapper.

This probe checks whether the timing decomposition in
simulation/benchmark_regional_runtime_effector.py materially changes effector's total runtime.
"""

import argparse
import csv
import inspect
import os
import statistics
import time
from collections import defaultdict

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np

from benchmark_regional_runtime_effector import (
    fit_rf,
    load_data,
    make_predict,
    make_space_partitioner,
    measure_regional_ale,
    measure_regional_pdp,
)

import effector


CELL_REGISTRY = {
    "small_toy_pdp": ("toy", "pdp", 1000, 10, 20, 2),
    "small_toy_ale": ("toy", "ale", 1000, 10, 20, 2),
    "mid_rf_pdp": ("rf", "pdp", 5000, 20, 20, 2),
    "mid_rf_ale": ("rf", "ale", 5000, 20, 20, 2),
    "high_rf_pdp_N20000": ("rf", "pdp", 20000, 20, 20, 2),
    "large_rf_ale_N20000": ("rf", "ale", 20000, 20, 20, 2),
    "high_toy_pdp_D100": ("toy", "pdp", 10000, 100, 20, 2),
    "high_toy_ale_split10": ("toy", "ale", 10000, 20, 20, 10),
    "high_rf_ale_K50": ("rf", "ale", 10000, 20, 50, 2),
    "high_rf_pdp_D50": ("rf", "pdp", 10000, 50, 20, 2),
    "high_rf_pdp_D100": ("rf", "pdp", 10000, 100, 20, 2),
}

DEFAULT_CELLS = [
    "small_toy_pdp",
    "small_toy_ale",
    "mid_rf_pdp",
    "mid_rf_ale",
    "large_rf_ale_N20000",
    "high_toy_pdp_D100",
    "high_toy_ale_split10",
    "high_rf_ale_K50",
]


def measure_official(effect, x, predict, d, resolution, n_split, min_node_size, numerical_features_grid_size):
    try:
        import effector.helpers as helpers
        helpers.NOF_INTERNAL_POINTS = int(resolution)
    except Exception:
        pass
    axis_limits = np.array([[-1.0] * d, [1.0] * d])
    space_partitioner = make_space_partitioner(n_split, min_node_size, numerical_features_grid_size)
    features = list(range(d))
    tic = time.time()
    if effect == "pdp":
        regional = effector.RegionalPDP(data=x, model=predict, axis_limits=axis_limits, nof_instances="all")
        params = inspect.signature(regional.fit).parameters
        kwargs = {"features": features, "space_partitioner": space_partitioner}
        if "centering" in params:
            kwargs["centering"] = True
        if "points_for_centering" in params:
            kwargs["points_for_centering"] = resolution
        if "points_for_mean_heterogeneity" in params:
            kwargs["points_for_mean_heterogeneity"] = resolution
        if "use_vectorized" in params:
            kwargs["use_vectorized"] = True
        regional.fit(**kwargs)
    else:
        binning = effector.axis_partitioning.Fixed(nof_bins=resolution, min_points_per_bin=0)
        regional = effector.RegionalALE(data=x, model=predict, axis_limits=axis_limits, nof_instances="all")
        params = inspect.signature(regional.fit).parameters
        kwargs = {"features": features, "space_partitioner": space_partitioner, "binning_method": binning}
        if "points_for_mean_heterogeneity" in params:
            kwargs["points_for_mean_heterogeneity"] = resolution
        regional.fit(**kwargs)
    return time.time() - tic


def summarize(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["cell"]].append(row)

    summary = []
    for cell, cell_rows in grouped.items():
        official = [float(x["total_time_sec"]) for x in cell_rows if x["method"] == "official_fit"]
        wrapper = [float(x["total_time_sec"]) for x in cell_rows if x["method"] == "wrapper_precompute_plus_split"]
        if not official or not wrapper:
            continue
        first = cell_rows[0]
        official_median = statistics.median(official)
        wrapper_median = statistics.median(wrapper)
        summary.append({
            "cell": cell,
            "model_type": first["model_type"],
            "effect": first["effect"],
            "N": first["N"],
            "D": first["D"],
            "resolution": first["resolution"],
            "n_split": first["n_split"],
            "official_fit_median_sec": official_median,
            "wrapper_total_median_sec": wrapper_median,
            "wrapper_over_official": wrapper_median / official_median if official_median else float("nan"),
        })
    return summary


def write_csv(path, rows):
    if not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def parse_cells(raw):
    if raw == "default":
        return DEFAULT_CELLS
    if raw == "all":
        return list(CELL_REGISTRY)
    cells = [x.strip() for x in raw.split(",") if x.strip()]
    unknown = sorted(set(cells) - set(CELL_REGISTRY))
    if unknown:
        raise ValueError("Unknown cell(s): {}".format(", ".join(unknown)))
    return cells


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", default="simulation/data/global_r_runtime")
    parser.add_argument("--outdir", default="simulation/results/runtime_runs/local_effector_wrapper_official_probe")
    parser.add_argument("--reps", type=int, default=2)
    parser.add_argument("--cells", default="default")
    parser.add_argument("--min-node-size", type=int, default=50)
    parser.add_argument("--numerical-features-grid-size", type=int, default=20)
    args = parser.parse_args()

    raw_path = os.path.join(args.outdir, "effector_wrapper_vs_official_raw.csv")
    summary_path = os.path.join(args.outdir, "effector_wrapper_vs_official_summary.csv")
    rows = read_csv(raw_path)
    data_cache = {}
    model_cache = {}
    selected_cells = parse_cells(args.cells)
    for cell in selected_cells:
        model_type, effect, n, d, resolution, n_split = CELL_REGISTRY[cell]
        key = (n, d)
        if key not in data_cache:
            x, y = load_data(args.datadir, n, d)
            if x is None:
                print(f"[skip] missing data for {cell}: N={n} D={d}", flush=True)
                continue
            data_cache[key] = (x, y)
        x, y = data_cache[key]

        model_key = (model_type, n, d)
        if model_key not in model_cache:
            print(f"[fit] {cell} model={model_type} N={n} D={d}", flush=True)
            model_cache[model_key] = fit_rf(x, y) if model_type == "rf" else None
        predict = make_predict(model_cache[model_key], model_type)
        wrapper = measure_regional_pdp if effect == "pdp" else measure_regional_ale

        print(f"[warmup] {cell}", flush=True)
        wrapper(x, predict, d, resolution, n_split, args.min_node_size, args.numerical_features_grid_size)
        measure_official(effect, x, predict, d, resolution, n_split, args.min_node_size, args.numerical_features_grid_size)

        for rep in range(1, args.reps + 1):
            print(f"[run] {cell} rep={rep} wrapper", flush=True)
            timing = wrapper(
                x,
                predict,
                d,
                resolution,
                n_split,
                args.min_node_size,
                args.numerical_features_grid_size,
            )
            rows.append({
                "cell": cell,
                "model_type": model_type,
                "effect": effect,
                "N": n,
                "D": d,
                "resolution": resolution,
                "n_split": n_split,
                "method": "wrapper_precompute_plus_split",
                "repetition": rep,
                "precompute_time_sec": timing["precompute"],
                "split_time_sec": timing["split"],
                "total_time_sec": timing["total"],
            })

            print(f"[run] {cell} rep={rep} official", flush=True)
            total = measure_official(
                effect,
                x,
                predict,
                d,
                resolution,
                n_split,
                args.min_node_size,
                args.numerical_features_grid_size,
            )
            rows.append({
                "cell": cell,
                "model_type": model_type,
                "effect": effect,
                "N": n,
                "D": d,
                "resolution": resolution,
                "n_split": n_split,
                "method": "official_fit",
                "repetition": rep,
                "precompute_time_sec": "",
                "split_time_sec": "",
                "total_time_sec": total,
            })
        write_csv(raw_path, rows)
        write_csv(summary_path, summarize(rows))

    summary = summarize(rows)
    write_csv(raw_path, rows)
    write_csv(summary_path, summary)

    print("\ncell,official_sec,wrapper_sec,wrapper/official", flush=True)
    for row in summary:
        print(
            f"{row['cell']},{row['official_fit_median_sec']:.6f},"
            f"{row['wrapper_total_median_sec']:.6f},{row['wrapper_over_official']:.3f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
