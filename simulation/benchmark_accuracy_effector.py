#!/usr/bin/env python3
"""Accuracy benchmark (effector): recover the oracle first split for the x2 regional effect."""

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd

try:
    import effector
    from effector.space_partitioning import Best
except ImportError:
    sys.stderr.write("effector not installed. Run: pip install effector\n")
    sys.exit(1)


N_VEC = [200, 500, 1000, 5000]
D_VEC = [5, 10, 20]
VARIANTS = ["num_0", "num_04", "cat"]
N_SEEDS = 30

FOI = 1
ORACLE_SPLIT_FEATURE = 2
N_INTERVALS = 20
MIN_NODE_SIZE = 40
NUMERICAL_GRID_SIZE = 40
MIN_HETEROGENEITY_DECREASE = 0.0


def parse_int_vec(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def parse_str_vec(s: str) -> list[str]:
    return [x.strip() for x in s.split(",") if x.strip()]


def variant_threshold(variant: str) -> float | None:
    if variant == "num_0":
        return 0.0
    if variant == "num_04":
        return 0.4
    return None


def load_csv(datadir: str, N: int, D: int, variant: str, seed: int) -> tuple[pd.DataFrame, np.ndarray]:
    fn = os.path.join(datadir, f"acc_N{N}_D{D}_{variant}_seed{seed}.csv")
    if not os.path.isfile(fn):
        sys.stderr.write(f"Missing dataset {fn}. Run: Rscript simulation/generate_accuracy_data.R\n")
        sys.exit(1)

    df = pd.read_csv(fn)
    if variant == "cat":
        df_num = df.copy()
        df_num["x3"] = np.where(df_num["x3"].astype(str).values == "0", 0.0, 1.0)
        X = df_num[[f"x{i + 1}" for i in range(D)]].values.astype(np.float64)
    else:
        X = df[[f"x{i + 1}" for i in range(D)]].values.astype(np.float64)
    return df, X


def fit_oracle_predictor(df: pd.DataFrame, X: np.ndarray, variant: str):
    x1 = X[:, 0]
    x2 = X[:, 1]
    x3 = X[:, 2]
    moderator_x1 = (x1 > 0).astype(np.float64)
    if variant == "cat":
        moderator_x3 = (x3 < 0.5).astype(np.float64)
    else:
        moderator_x3 = (x3 <= float(variant_threshold(variant))).astype(np.float64)

    design = np.column_stack(
        [
            np.ones(X.shape[0], dtype=np.float64),
            x1,
            x2,
            x2 * moderator_x3,
            x2 * moderator_x1,
        ]
    )
    coef, _, _, _ = np.linalg.lstsq(design, df["y"].to_numpy(dtype=np.float64), rcond=None)

    def pred(Xm: np.ndarray) -> np.ndarray:
        x1_new = Xm[:, 0]
        x2_new = Xm[:, 1]
        x3_new = Xm[:, 2]
        moderator_x1_new = (x1_new > 0).astype(np.float64)
        if variant == "cat":
            moderator_x3_new = (x3_new < 0.5).astype(np.float64)
        else:
            moderator_x3_new = (x3_new <= float(variant_threshold(variant))).astype(np.float64)
        return (
            coef[0]
            + coef[1] * x1_new
            + coef[2] * x2_new
            + coef[3] * (x2_new * moderator_x3_new)
            + coef[4] * (x2_new * moderator_x1_new)
        )

    return pred


def axis_limits_from_X(X: np.ndarray) -> np.ndarray:
    lo = X.min(axis=0)
    hi = X.max(axis=0)
    pad = 1e-6 * (hi - lo + 1e-9)
    return np.vstack([lo - pad, hi + pad])


def true_left_mask(X: np.ndarray, variant: str) -> np.ndarray:
    x3 = X[:, 2]
    if variant == "cat":
        return x3 < 0.5
    return x3 <= float(variant_threshold(variant))


def first_split_props(part) -> tuple[int | None, object | None, str | None]:
    if part is None:
        return None, None, None

    important_splits = getattr(part, "important_splits", None)
    if isinstance(important_splits, list) and len(important_splits) > 0:
        split_0 = important_splits[0]
        return split_0["foc_index"], split_0["foc_split_position"], split_0["foc_type"]

    tree = getattr(part, "splits_tree", None)
    if tree is not None and len(getattr(tree, "nodes", [])) > 1:
        for node in tree.nodes:
            info = getattr(node, "info", None) or {}
            if info.get("level") == 1 and "foc_index" in info:
                return info["foc_index"], info["foc_split_position"], info.get("foc_type")
    return None, None, None


def split_point_mae(foc_index: int | None, position, variant: str) -> float | None:
    if variant == "cat" or foc_index != ORACLE_SPLIT_FEATURE or position is None:
        return None
    return abs(float(position) - float(variant_threshold(variant)))


def node_acc(X: np.ndarray, variant: str, foc_index: int | None, position, foc_type: str | None) -> float:
    if foc_index is None or position is None:
        return float("nan")

    truth_left = true_left_mask(X, variant)
    if foc_type == "cat":
        pred_left = X[:, foc_index] == position
    else:
        pred_left = X[:, foc_index] < float(position)

    acc_left = np.mean(pred_left == truth_left)
    acc_right = np.mean((~pred_left) == truth_left)
    return float(max(acc_left, acc_right))


def run_regional_pdp(
    X: np.ndarray,
    predict,
    axis_limits: np.ndarray,
    feature_names: list[str],
    space_partitioner: Best,
):
    regional = effector.RegionalPDP(
        data=X,
        model=predict,
        axis_limits=axis_limits,
        nof_instances="all",
        feature_names=feature_names,
    )
    regional.fit(
        features=[FOI],
        candidate_conditioning_features="all",
        space_partitioner=space_partitioner,
    )
    return regional.partitioners[f"feature_{FOI}"]


def run_regional_ale(
    X: np.ndarray,
    predict,
    axis_limits: np.ndarray,
    feature_names: list[str],
    space_partitioner: Best,
):
    binning = effector.axis_partitioning.Fixed(nof_bins=N_INTERVALS, min_points_per_bin=0)
    regional = effector.RegionalALE(
        data=X,
        model=predict,
        axis_limits=axis_limits,
        nof_instances="all",
        feature_names=feature_names,
    )
    regional.fit(
        features=[FOI],
        candidate_conditioning_features="all",
        space_partitioner=space_partitioner,
        binning_method=binning,
        points_for_mean_heterogeneity=N_INTERVALS,
    )
    return regional.partitioners[f"feature_{FOI}"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", default="simulation/data/accuracy")
    parser.add_argument("--outdir", default="simulation/results/accuracy")
    parser.add_argument("--n-seeds", type=int, default=N_SEEDS)
    parser.add_argument("--N-vec", default=None)
    parser.add_argument("--D-vec", default=None)
    parser.add_argument("--variants", default=None)
    args, _ = parser.parse_known_args()

    n_vec = parse_int_vec(args.N_vec) if args.N_vec else N_VEC
    d_vec = parse_int_vec(args.D_vec) if args.D_vec else D_VEC
    variants = parse_str_vec(args.variants) if args.variants else VARIANTS
    os.makedirs(args.outdir, exist_ok=True)

    space_partitioner = Best(
        max_depth=1,
        min_samples_leaf=MIN_NODE_SIZE,
        min_heterogeneity_decrease_pcg=MIN_HETEROGENEITY_DECREASE,
        numerical_features_grid_size=NUMERICAL_GRID_SIZE,
    )

    rows: list[dict[str, object]] = []
    for variant in variants:
        for N in n_vec:
            for D in d_vec:
                feature_names = [f"x{i + 1}" for i in range(D)]
                for s in range(1, args.n_seeds + 1):
                    seed = 1000 + s
                    df, X = load_csv(args.datadir, N, D, variant, seed)
                    axis_limits = axis_limits_from_X(X)
                    predict = fit_oracle_predictor(df, X, variant)

                    part_pdp = run_regional_pdp(X, predict, axis_limits, feature_names, space_partitioner)
                    foc_index, position, foc_type = first_split_props(part_pdp)
                    rows.append(
                        {
                            "package": "effector",
                            "method": "effector_rpdp",
                            "variant": variant,
                            "N": N,
                            "D": D,
                            "seed": seed,
                            "foi_feature": "x2",
                            "oracle_split_feature": "x3",
                            "selected_split_feature": feature_names[foc_index] if foc_index is not None else "",
                            "selected_split_value": "" if position is None else str(position),
                            "split_feat_correct": foc_index == ORACLE_SPLIT_FEATURE,
                            "split_pt_error": (
                                "" if split_point_mae(foc_index, position, variant) is None
                                else split_point_mae(foc_index, position, variant)
                            ),
                            "node_acc": node_acc(X, variant, foc_index, position, foc_type),
                        }
                    )

                    part_ale = run_regional_ale(X, predict, axis_limits, feature_names, space_partitioner)
                    foc_index, position, foc_type = first_split_props(part_ale)
                    rows.append(
                        {
                            "package": "effector",
                            "method": "effector_rale",
                            "variant": variant,
                            "N": N,
                            "D": D,
                            "seed": seed,
                            "foi_feature": "x2",
                            "oracle_split_feature": "x3",
                            "selected_split_feature": feature_names[foc_index] if foc_index is not None else "",
                            "selected_split_value": "" if position is None else str(position),
                            "split_feat_correct": foc_index == ORACLE_SPLIT_FEATURE,
                            "split_pt_error": (
                                "" if split_point_mae(foc_index, position, variant) is None
                                else split_point_mae(foc_index, position, variant)
                            ),
                            "node_acc": node_acc(X, variant, foc_index, position, foc_type),
                        }
                    )

    outp = os.path.join(args.outdir, "accuracy_effector.csv")
    with open(outp, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Written: {outp}")


if __name__ == "__main__":
    main()
