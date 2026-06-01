#!/usr/bin/env python3
"""effector split-search benchmark. Global precomputation is excluded from time_sec."""

import argparse
import csv
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OMP_THREAD_LIMIT", "1")
os.environ.setdefault("OMP_PROC_BIND", "FALSE")
os.environ.setdefault("KMP_INIT_AT_FORK", "FALSE")
os.environ.setdefault("KMP_AFFINITY", "disabled")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
import pandas as pd

try:
    import effector
except ImportError:
    sys.stderr.write("effector not installed. Run: pip install effector\n")
    sys.exit(1)


def parse_int_vec(raw):
    return [int(x) for x in raw.split(",") if x.strip()]


def load_data(datadir, n, d, seed=21):
    path = os.path.join(datadir, "benchmark_N{}_D{}_seed{}.csv".format(n, d, seed))
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    return df[["x{}".format(i + 1) for i in range(d)]].values


def toy_pred_fun(x):
    x1, x2, x3 = x[:, 0], x[:, 1], x[:, 2]
    return 5.0 * x1 + 5.0 * x2 + np.where(x3 > 0, 10.0 * x1, 0.0) - np.where(x3 > 0, 10.0 * x2, 0.0)


def precompute_regional_pdp(x, d, n_grid):
    axis_limits = np.array([[-1.0] * d, [1.0] * d])
    model = toy_pred_fun
    regional = effector.RegionalPDP(data=x, model=model, axis_limits=axis_limits, nof_instances="all")
    regional.y_ice = {}
    for feat in range(d):
        pdp = effector.PDP(data=x, model=model, axis_limits=axis_limits, nof_instances="all")
        try:
            pdp.fit(features=feat, centering=True, points_for_centering=n_grid)
        except TypeError:
            pdp.fit(features=feat)
        xs = np.linspace(axis_limits[0, feat], axis_limits[1, feat], n_grid)
        y_ice = pdp.eval(feature=feat, xs=xs, heterogeneity=True, return_all=True)
        regional.y_ice["feature_" + str(feat)] = y_ice.T
    return regional


def precompute_regional_ale(x, d, n_intervals):
    axis_limits = np.array([[-1.0] * d, [1.0] * d])
    model = toy_pred_fun
    binning = effector.axis_partitioning.Fixed(nof_bins=n_intervals, min_points_per_bin=0)
    regional = effector.RegionalALE(data=x, model=model, axis_limits=axis_limits, nof_instances="all")
    regional.global_data_effect = {}
    regional.global_bin_limits = {}
    for feat in range(d):
        global_ale = effector.ALE(data=x, model=model, axis_limits=axis_limits, nof_instances="all")
        global_ale.fit(features=feat, binning_method=binning, centering=False)
        key = "feature_" + str(feat)
        regional.global_data_effect[key] = global_ale.data_effect_ale[key]
        regional.global_bin_limits[key] = global_ale.bin_limits[key]
    return regional


def time_split(regional, method, d, n_intervals):
    space_partitioner = effector.space_partitioning.Best(max_depth=2, min_samples_leaf=50)
    tic = time.time()
    for feat in range(d):
        if method == "regional_pdp":
            heter = regional._create_heterogeneity_function(
                foi=feat,
                min_points=space_partitioner.min_points_per_subregion,
            )
        else:
            heter = regional._create_heterogeneity_function(
                feat,
                space_partitioner.min_points_per_subregion,
                n_intervals,
            )
        regional._fit_feature(feat, heter, space_partitioner, "all")
    return time.time() - tic


def record(method, n, d, n_grid, n_intervals, repetition, time_sec, status="ok", error_message=""):
    return {
        "module": "split_search",
        "package": "effector",
        "impl": "split",
        "method": method,
        "model_type": "toy",
        "N": n,
        "D": d,
        "n_grid": n_grid if "pdp" in method else "",
        "n_intervals": n_intervals if "ale" in method else "",
        "repetition": repetition,
        "time_sec": time_sec,
        "status": status,
        "error_message": error_message,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", default="simulation/data/global_r_runtime")
    parser.add_argument("--outdir", default="simulation/results/split_search_runtime")
    parser.add_argument("--reps", type=int, default=20)
    parser.add_argument("--N-vec", default="500,1000,2500,5000")
    parser.add_argument("--D", type=int, default=10)
    parser.add_argument("--n-grid", type=int, default=20)
    parser.add_argument("--n-intervals", type=int, default=20)
    parser.add_argument("--fail-fast", default="false")
    args = parser.parse_args()

    fail_fast = args.fail_fast.lower() in {"true", "1", "yes", "y", "on"}
    n_vec = parse_int_vec(args.N_vec)
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    for n in n_vec:
        x = load_data(args.datadir, n, args.D)
        if x is None:
            continue
        print("=== effector split: N={}, D={} ===".format(n, args.D), flush=True)

        for method in ["regional_pdp", "regional_ale"]:
            print("  {}".format(method), flush=True)
            try:
                if method == "regional_pdp":
                    regional = precompute_regional_pdp(x, args.D, args.n_grid)
                else:
                    regional = precompute_regional_ale(x, args.D, args.n_intervals)
                time_split(regional, method, args.D, args.n_intervals)
            except Exception as exc:
                if fail_fast:
                    raise
                print("    warmup skipped/failed: {}".format(exc), flush=True)

            for rep in range(1, args.reps + 1):
                try:
                    if method == "regional_pdp":
                        regional = precompute_regional_pdp(x, args.D, args.n_grid)
                    else:
                        regional = precompute_regional_ale(x, args.D, args.n_intervals)
                    elapsed = time_split(regional, method, args.D, args.n_intervals)
                    rows.append(record(method, n, args.D, args.n_grid, args.n_intervals, rep, elapsed))
                except Exception as exc:
                    if fail_fast:
                        raise
                    rows.append(record(
                        method,
                        n,
                        args.D,
                        args.n_grid,
                        args.n_intervals,
                        rep,
                        float("nan"),
                        status="error",
                        error_message=str(exc),
                    ))

    out = os.path.join(args.outdir, "split_search_effector_toy.csv")
    cols = [
        "module", "package", "impl", "method", "model_type", "N", "D",
        "n_grid", "n_intervals", "repetition", "time_sec", "status", "error_message",
    ]
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols)
        writer.writeheader()
        writer.writerows(rows)
    print("Written: {}".format(out), flush=True)


if __name__ == "__main__":
    main()
