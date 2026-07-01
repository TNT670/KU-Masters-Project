"""
Shared utilities for the cancer-tissue classifier notebooks
(breast_cancer.ipynb, lung_cancer.ipynb, colon_cancer.ipynb).

Centralizes the three things that used to be copy-pasted and hand-edited
across all 11 original notebooks:

  1. Reading the per-cancer JSON label config.
  2. Building the sample label array from that config (complete or
     leave-one-out mode), with a hard length check against the declared
     per-GSE sample counts.
  3. Resolving the preprocessed `.dat` filenames for each run mode,
     following the naming convention emitted by `preprocess.R`.

Label config schema (see breast_labels.json / lung_labels.json /
colon_labels.json). Each dataset's "labels" block is exactly one of:

  uniform     : {"type": "uniform", "value": 0|1}
  blocks      : {"type": "blocks", "blocks": [{"value": v, "count": n}, ...]}
  alternating : {"type": "alternating", "first_value": 0|1}
                -> label[i] = (i + first_value) % 2
  text        : {"type": "text", "file": "<name>.txt",
                 "match": "contains"|"endswith",
                 "patterns": ["IDC", ...], "match_label": 0|1}
                -> a line matching ANY pattern gets match_label,
                   every other line gets (1 - match_label).
"""

import json
import os
import numpy as np
import pandas as pd


def load_config(path):
    """Load a cancer label config JSON file."""
    with open(path) as f:
        return json.load(f)


def _build_block(ds, base_dir):
    """Build the label array for a single GSE dataset entry."""
    spec = ds["labels"]
    t = spec["type"]

    if t == "uniform":
        arr = [spec["value"]] * ds["n_samples"]

    elif t == "blocks":
        arr = []
        for b in spec["blocks"]:
            arr += [b["value"]] * b["count"]

    elif t == "alternating":
        first = spec["first_value"]
        arr = [(i + first) % 2 for i in range(ds["n_samples"])]

    elif t == "text":
        path = os.path.join(base_dir, spec["file"])
        match_label = spec["match_label"]
        other_label = 1 - match_label
        patterns = spec["patterns"]
        mode = spec["match"]
        arr = []
        with open(path) as f:
            for line in f:
                s = line.rstrip("\n")
                if mode == "contains":
                    hit = any(p in s for p in patterns)
                elif mode == "endswith":
                    hit = any(s.endswith(p) for p in patterns)
                else:
                    raise ValueError(f"unknown text match mode: {mode!r}")
                arr.append(match_label if hit else other_label)

    else:
        raise ValueError(f"unknown label type: {t!r}")

    if len(arr) != ds["n_samples"]:
        raise AssertionError(
            f'GSE{ds["gse"]}: built {len(arr)} labels but n_samples={ds["n_samples"]}'
        )
    return np.array(arr)


def build_labels(config, base_dir, held_out_gse=None):
    """
    Build the label array(s) for a run.

    Returns (train_labels, test_labels).
      - Complete mode (held_out_gse is None): test_labels is None and
        train_labels covers every GSE in config order.
      - LOO mode: the held-out GSE's block becomes test_labels and is
        omitted from train_labels.

    base_dir is the directory containing the text-label files (the cancer dir).
    """
    held = str(held_out_gse) if held_out_gse is not None else None
    gse_ids = [ds["gse"] for ds in config["datasets"]]
    if held is not None and held not in gse_ids:
        raise ValueError(f"HELD_OUT_GSE {held!r} not in config GSEs {gse_ids}")

    train_blocks = []
    test_labels = None
    for ds in config["datasets"]:
        block = _build_block(ds, base_dir)
        if held is not None and ds["gse"] == held:
            test_labels = block
        else:
            train_blocks.append(block)

    train_labels = np.concatenate(train_blocks) if train_blocks else np.array([])
    return train_labels, test_labels


def dat_paths(config, base_dir, held_out_gse=None):
    """
    Resolve preprocessed .dat file paths, following the preprocess.R naming:
      complete : {prefix}.dat
      LOO      : {prefix}_loo{GSE}_train.dat  and  {prefix}_loo{GSE}_test.dat
    """
    prefix = config["output_prefix"]
    if held_out_gse is None:
        return {"complete": os.path.join(base_dir, f"{prefix}.dat")}
    g = str(held_out_gse)
    return {
        "train": os.path.join(base_dir, f"{prefix}_loo{g}_train.dat"),
        "test": os.path.join(base_dir, f"{prefix}_loo{g}_test.dat"),
    }


def select_uncorrelated(df, labels, threshold=0.85):
    """
    Remove redundant highly-correlated features.

    Builds the absolute Pearson correlation matrix via the fast vectorized
    form (centered, L2-normalized columns; X.T @ X), treats each connected
    component of the "|corr| > threshold" graph as a cluster of mutually
    redundant features, and from each cluster keeps the single member most
    correlated with `labels` (point-biserial correlation), dropping the rest.
    Features not correlated with anything above the threshold are all kept.

    This replaces the original column-order one-liner, which (a) kept an
    arbitrary cluster member determined by gene-ID order and (b) could
    over-drop transitive chains (dropping a feature for correlating with one
    that was itself already dropped). Computed on `df` alone, so it stays
    leak-free when `df` is the training set.

    Returns (keep_cols, drop_cols) as lists of column labels.
    """
    cols = list(df.columns)
    n = len(cols)
    X = df.to_numpy(dtype=float)

    # Centered, L2-normalized columns -> Xs.T @ Xs is the Pearson matrix.
    Xc = X - X.mean(axis=0)
    norms = np.sqrt((Xc ** 2).sum(axis=0))
    norms[norms == 0] = 1.0  # guard zero-variance columns (their corr -> 0)
    Xs = Xc / norms

    corr = np.abs(Xs.T @ Xs)
    np.fill_diagonal(corr, 0.0)

    # |correlation of each feature with the label| -> representative score.
    y = np.asarray(labels, dtype=float)
    yc = y - y.mean()
    yn = np.sqrt((yc ** 2).sum())
    rel = np.abs(Xs.T @ (yc / yn)) if yn > 0 else np.zeros(n)

    # Union-find over the connected components of the > threshold graph.
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i, j in np.argwhere(np.triu(corr > threshold, k=1)):
        union(int(i), int(j))

    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)

    keep_idx = set()
    for members in groups.values():
        if len(members) == 1:
            keep_idx.add(members[0])
        else:
            keep_idx.add(max(members, key=lambda m: rel[m]))

    keep_cols = [cols[i] for i in range(n) if i in keep_idx]
    drop_cols = [cols[i] for i in range(n) if i not in keep_idx]
    return keep_cols, drop_cols


def subsample(X, y, p, stratify=True):
    """
    Return a random fraction p of (X, y) for a training-size sweep.

    p >= 1 returns the full set unchanged. With stratify=True the class balance
    is preserved (so iteration-to-iteration variance reflects model instability,
    not class-ratio luck). With stratify=False, draws are re-tried until both
    classes are present (matching the original notebooks' guard).
    """
    from sklearn.model_selection import train_test_split
    if p >= 1.0:
        return X, y
    if stratify:
        Xs, _, ys, _ = train_test_split(X, y, train_size=p, stratify=y)
        return Xs, ys
    while True:
        Xs, _, ys, _ = train_test_split(X, y, train_size=p)
        if len(np.unique(ys)) > 1:
            return Xs, ys


def score(y_true, preds):
    """Return [accuracy, precision, recall, f1] for binary predictions."""
    from sklearn.metrics import precision_recall_fscore_support as prfs, accuracy_score
    p, r, f, _ = prfs(y_true, preds, average="binary", zero_division=0.0)
    return [accuracy_score(y_true, preds), p, r, f]


def plot_depth_sweep(means, percents, hp_values, hp_label, title):
    """
    Grid of 4 subplots (accuracy/precision/recall/F1), each plotting the metric
    against a hyperparameter, with one line per training-size fraction.
    `means` is shaped [len(percents)][len(hp_values)][4].
    """
    import matplotlib.pyplot as plt
    means = np.array(means)
    metrics = ["Accuracy", "Precision", "Recall", "F-Score"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for m, ax in enumerate(axes.flat):
        for i, p in enumerate(percents):
            ax.plot(hp_values, means[i][:, m], label=f"{p*100:.0f}%")
        ax.set_xlabel(hp_label)
        ax.set_ylabel(f"Avg {metrics[m]}")
        ax.set_title(metrics[m])
    fig.suptitle(title)
    fig.tight_layout()
    handles, labels = axes.flat[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", title="Train %", bbox_to_anchor=(1.0, 0.5))
    plt.show()


def plot_increments(means, percents, title):
    """
    Plot accuracy/precision/recall/F1 against training-size fraction.
    `means` is shaped [len(percents)][4]. X-axis runs 100% (left) -> 10% (right).
    """
    import matplotlib.pyplot as plt
    means = np.array(means)
    metrics = ["Accuracy", "Precision", "Recall", "F-Score"]
    x = [p * 100 for p in percents]
    for m in range(4):
        plt.plot(x, means[:, m], marker="o", label=metrics[m])
    plt.gca().invert_xaxis()  # 100% on the left, like the original plots
    plt.xlabel("Percent of training data used")
    plt.ylabel("Average score")
    plt.title(title)
    plt.legend()
    plt.show()


def load_dat(path):
    """
    Load a space-delimited preprocessed table and remove the Batch column.
    Returns (dataframe_of_expression, batch_series).

    The complete-run and LOO-train tables carry a trailing 'Batch' column; the
    LOO-test table (a single held-out GSE, written by preprocess.R) does not, so
    the column is popped only when present and `batch` is None otherwise.
    Row order matches the GSE order in the label config.
    """
    df = pd.read_csv(path, sep=" ")
    batch = df.pop("Batch") if "Batch" in df.columns else None
    return df, batch
