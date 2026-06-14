#!/usr/bin/env python3
"""effector regional PDP/ALE runtime benchmark.

Measures regional precompute, split-search, and total timings separately.
"""

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
os.environ.setdefault("PYTHONHASHSEED", "21")

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

try:
    import effector
except ImportError:
    sys.stderr.write("effector not installed. Run: pip install effector\n")
    sys.exit(1)


RF_CONFIG = {
    "n_estimators": 100,
    "criterion": "squared_error",
    "max_features": 1.0,
    "max_depth": None,
    "min_samples_split": 2,
    "min_samples_leaf": 1,
    "min_weight_fraction_leaf": 0.0,
    "max_leaf_nodes": None,
    "min_impurity_decrease": 0.0,
    "bootstrap": True,
    "max_samples": None,
    "ccp_alpha": 0.0,
    "random_state": 21,
    "n_jobs": 1,
}


def parse_int_vec(raw):
    return [int(x) for x in raw.split(",") if x.strip()]


def parse_chr_vec(raw):
    return [x.strip() for x in raw.split(",") if x.strip()]


def parse_flag(raw):
    return str(raw).strip().lower() in {"true", "1", "yes", "y", "on"}


def load_data(datadir, n, d, seed=21):
    path = os.path.join(datadir, "benchmark_N{}_D{}_seed{}.csv".format(n, d, seed))
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path)
    x = df[["x{}".format(i + 1) for i in range(d)]].values
    y = df["y"].values
    return x, y


def fit_rf(x, y):
    return RandomForestRegressor(**RF_CONFIG).fit(x, y)


def rf_pred_fun(model, x):
    return model.predict(x)


def toy_pred_fun(x):
    x1, x2, x3 = x[:, 0], x[:, 1], x[:, 2]
    return 5.0 * x1 + 5.0 * x2 + np.where(x3 > 0, 10.0 * x1, 0.0) - np.where(x3 > 0, 10.0 * x2, 0.0)


def quantile_points(x, n_points):
    probs = np.linspace(0.0, 1.0, n_points)
    return np.quantile(x, probs)


def make_predict(model, model_type):
    if model_type == "rf":
        return lambda x: rf_pred_fun(model, x)
    return toy_pred_fun


def make_space_partitioner(n_split, min_node_size):
    return effector.space_partitioning.Best(max_depth=n_split, min_samples_leaf=min_node_size)


def measure_regional_pdp(x, predict, d, resolution, n_split, min_node_size):
    axis_limits = np.array([[-1.0] * d, [1.0] * d])
    space_partitioner = make_space_partitioner(n_split, min_node_size)
    features = list(range(d))
    regional = effector.RegionalPDP(data=x, model=predict, axis_limits=axis_limits, nof_instances="all")
    regional.y_ice = {}

    tic = time.time()
    for feat in features:
        pdp = effector.PDP(data=x, model=predict, axis_limits=axis_limits, nof_instances="all")
        try:
            pdp.fit(features=feat, centering=True, points_for_centering=resolution)
        except TypeError:
            pdp.fit(features=feat)
        xs = quantile_points(x[:, feat], resolution)
        y_ice = pdp.eval(feature=feat, xs=xs, heterogeneity=True, return_all=True)
        regional.y_ice["feature_" + str(feat)] = y_ice.T
    precompute = time.time() - tic

    tic = time.time()
    for feat in features:
        heterogeneity = regional._create_heterogeneity_function(
            foi=feat,
            min_points=space_partitioner.min_points_per_subregion,
        )
        regional._fit_feature(feat, heterogeneity, space_partitioner, "all")
    split = time.time() - tic
    return {"precompute": precompute, "split": split, "total": precompute + split}


def measure_regional_ale(x, predict, d, resolution, n_split, min_node_size):
    axis_limits = np.array([[-1.0] * d, [1.0] * d])
    space_partitioner = make_space_partitioner(n_split, min_node_size)
    binning = effector.axis_partitioning.Fixed(nof_bins=resolution, min_points_per_bin=0)
    features = list(range(d))
    regional = effector.RegionalALE(data=x, model=predict, axis_limits=axis_limits, nof_instances="all")
    regional.global_data_effect = {}
    regional.global_bin_limits = {}

    tic = time.time()
    for feat in features:
        global_ale = effector.ALE(data=x, model=predict, axis_limits=axis_limits, nof_instances="all")
        global_ale.fit(features=feat, binning_method=binning, centering=False)
        key = "feature_" + str(feat)
        regional.global_data_effect[key] = global_ale.data_effect_ale[key]
        regional.global_bin_limits[key] = global_ale.bin_limits[key]
    precompute = time.time() - tic

    tic = time.time()
    for feat in features:
        heterogeneity = regional._create_heterogeneity_function(
            feat,
            space_partitioner.min_points_per_subregion,
            resolution,
        )
        regional._fit_feature(feat, heterogeneity, space_partitioner, "all")
    split = time.time() - tic
    return {"precompute": precompute, "split": split, "total": precompute + split}


def make_cells(args):
    rows = []
    if "vs_N" in args.sub_experiments:
        for n in args.N_vec:
            rows.append(("vs_N", n, args.fixed_D, args.resolution, args.n_split))
    if "vs_D" in args.sub_experiments:
        for d in args.D_vec:
            rows.append(("vs_D", args.fixed_N, d, args.resolution, args.n_split))
    if "vs_res" in args.sub_experiments:
        for resolution in args.resolution_vec:
            rows.append(("vs_res", args.fixed_N, args.fixed_D, resolution, args.n_split))
    if "vs_split" in args.sub_experiments:
        for n_split in args.n_split_vec:
            rows.append(("vs_split", args.fixed_N, args.fixed_D, args.resolution, n_split))
    return list(dict.fromkeys(rows))


def record(args, model_type, effect, cell, repetition, timing=None, status="ok", error_message=""):
    sub_experiment, n, d, resolution, n_split = cell
    return {
        "module": "regional_runtime",
        "package": "effector",
        "impl": "default",
        "effect": effect,
        "method": "regional_{}".format(effect),
        "model_type": model_type,
        "sub_experiment": sub_experiment,
        "N": n,
        "D": d,
        "resolution": resolution,
        "n_grid": resolution if effect == "pdp" else "",
        "n_intervals": resolution if effect == "ale" else "",
        "n_split": n_split,
        "repetition": repetition,
        "precompute_time_sec": "" if timing is None else timing["precompute"],
        "split_time_sec": "" if timing is None else timing["split"],
        "total_time_sec": "" if timing is None else timing["total"],
        "status": status,
        "error_message": error_message,
    }


def run_model(args, model_type):
    cells = make_cells(args)
    data_cache = {}
    model_cache = {}
    for _, n, d, _, _ in cells:
        key = (n, d)
        if key in data_cache:
            continue
        x, y = load_data(args.datadir, n, d)
        if x is None:
            continue
        data_cache[key] = (x, y)
        model_cache[key] = fit_rf(x, y) if model_type == "rf" else None

    rows = []
    for effect in ["pdp", "ale"]:
        runner = measure_regional_pdp if effect == "pdp" else measure_regional_ale
        for cell in cells:
            sub_experiment, n, d, resolution, n_split = cell
            key = (n, d)
            if key not in data_cache:
                continue
            x, _ = data_cache[key]
            predict = make_predict(model_cache[key], model_type)
            msg = (
                "[{}] effector regional {} {} N={} D={} res={} n_split={}".format(
                    model_type, effect, sub_experiment, n, d, resolution, n_split
                )
            )
            print(msg + " | start", flush=True)
            try:
                runner(x, predict, d, resolution, n_split, args.min_node_size)
            except Exception as exc:
                if args.fail_fast:
                    raise
                print(msg + " | warmup skipped/failed: {}".format(exc), flush=True)

            for repetition in range(1, args.reps + 1):
                try:
                    timing = runner(x, predict, d, resolution, n_split, args.min_node_size)
                    rows.append(record(args, model_type, effect, cell, repetition, timing=timing))
                except Exception as exc:
                    if args.fail_fast:
                        raise
                    rows.append(record(
                        args,
                        model_type,
                        effect,
                        cell,
                        repetition,
                        status="error",
                        error_message=str(exc),
                    ))
            print(msg + " | done", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", default="simulation/data/global_r_runtime")
    parser.add_argument("--outdir", default="simulation/results/regional_runtime")
    parser.add_argument("--reps", type=int, default=30)
    parser.add_argument("--N-vec", default="5000,10000,25000,50000")
    parser.add_argument("--D-vec", default="10,20,50,100")
    parser.add_argument("--fixed-N", type=int, default=10000)
    parser.add_argument("--fixed-D", type=int, default=20)
    parser.add_argument("--resolution", type=int, default=20)
    parser.add_argument("--resolution-vec", default="10,20,50")
    parser.add_argument("--n-split", type=int, default=2)
    parser.add_argument("--n-split-vec", default="2,5,8,10")
    parser.add_argument("--min-node-size", type=int, default=50)
    parser.add_argument("--models", default="rf,toy")
    parser.add_argument("--sub-experiments", default="vs_N,vs_D,vs_res,vs_split")
    parser.add_argument("--output-suffix", default="")
    parser.add_argument("--fail-fast", default="false")
    args = parser.parse_args()
    args.N_vec = parse_int_vec(args.N_vec)
    args.D_vec = parse_int_vec(args.D_vec)
    args.resolution_vec = parse_int_vec(args.resolution_vec)
    args.n_split_vec = parse_int_vec(args.n_split_vec)
    args.models = parse_chr_vec(args.models)
    args.sub_experiments = parse_chr_vec(args.sub_experiments)
    args.fail_fast = parse_flag(args.fail_fast)
    invalid_models = sorted(set(args.models) - {"rf", "toy"})
    if invalid_models:
        raise ValueError("Unsupported regional model type(s): {}".format(", ".join(invalid_models)))
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    for model_type in args.models:
        rows.extend(run_model(args, model_type))

    filename = "regional_runtime_effector"
    if args.output_suffix:
        filename += "_" + args.output_suffix
    filename += ".csv"
    out = os.path.join(args.outdir, filename)
    columns = [
        "module", "package", "impl", "effect", "method", "model_type", "sub_experiment",
        "N", "D", "resolution", "n_grid", "n_intervals", "n_split", "repetition",
        "precompute_time_sec", "split_time_sec", "total_time_sec", "status", "error_message",
    ]
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print("Written: {}".format(out), flush=True)


if __name__ == "__main__":
    main()
