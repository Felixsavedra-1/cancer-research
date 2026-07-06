"""Ranking-quality metrics for the benchmark harness.

Pure, network-free, dependency-light (NumPy only) so the benchmark's evaluation
is itself unit-testable. Everything here scores a ranking of residues against a
binary gold standard (1 = known driver residue, 0 = other candidate).

AUROC and average precision (AUPRC) are the two headline numbers; AUPRC is the
honest one for driver discovery because positives are rare, so a high AUROC can
still hide poor precision. See docs/METHODS.md → Benchmark.
"""

from __future__ import annotations

import numpy as np


def _as_arrays(labels, scores) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(labels, dtype=float)
    s = np.asarray(scores, dtype=float)
    if y.shape != s.shape:
        raise ValueError("labels and scores must have the same length")
    return y, s


def roc_auc(labels, scores) -> float:
    """Area under the ROC curve via the Mann-Whitney U (rank-sum) identity.

    Ties are handled with average ranks. Undefined (returns ``nan``) when one
    class is absent. Equivalent to the probability a random positive outranks a
    random negative.
    """
    y, s = _as_arrays(labels, scores)
    n_pos = float((y == 1).sum())
    n_neg = float((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    sorted_s = s[order]
    i = 0
    n = len(s)
    while i < n:
        j = i
        while j + 1 < n and sorted_s[j + 1] == sorted_s[i]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        ranks[order[i : j + 1]] = avg_rank
        i = j + 1
    sum_pos_ranks = ranks[y == 1].sum()
    auc = (sum_pos_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def average_precision(labels, scores) -> float:
    """Average precision (area under the precision-recall curve).

    Computed as the recall-weighted mean of precision at each threshold where a
    positive is retrieved: ``AP = Σ (R_k − R_{k−1}) · P_k``. Ties in score are
    broken pessimistically (positives placed after negatives of equal score) so
    the number cannot be inflated by tie ordering.
    """
    y, s = _as_arrays(labels, scores)
    n_pos = float((y == 1).sum())
    if n_pos == 0:
        return float("nan")
    # Pessimistic tie-break: sort by descending score, then ascending label, so
    # equal-scored negatives are ranked ahead of positives.
    order = np.lexsort((y, -s))
    y_sorted = y[order]
    tp = np.cumsum(y_sorted)
    fp = np.cumsum(1.0 - y_sorted)
    precision = tp / np.maximum(tp + fp, 1.0)
    recall = tp / n_pos
    # Sum precision at the ranks where a new positive is retrieved.
    is_pos = y_sorted == 1
    recall_prev = np.concatenate(([0.0], recall[:-1]))
    delta_recall = recall - recall_prev
    return float(np.sum(precision[is_pos] * delta_recall[is_pos]))


def precision_at_k(labels, scores, k: int) -> float:
    """Fraction of the top-``k`` ranked items that are positives."""
    y, s = _as_arrays(labels, scores)
    if k <= 0 or len(y) == 0:
        return float("nan")
    k = min(k, len(y))
    top = np.argsort(-s, kind="mergesort")[:k]
    return float(y[top].sum() / k)


def recall_at_k(labels, scores, k: int) -> float:
    """Fraction of all positives captured within the top-``k`` ranked items."""
    y, s = _as_arrays(labels, scores)
    n_pos = float((y == 1).sum())
    if n_pos == 0 or len(y) == 0:
        return float("nan")
    k = min(k, len(y))
    top = np.argsort(-s, kind="mergesort")[:k]
    return float(y[top].sum() / n_pos)


def prevalence(labels) -> float:
    """Positive base rate — the AUPRC a random ranker would achieve."""
    y = np.asarray(labels, dtype=float)
    return float((y == 1).mean()) if len(y) else float("nan")
