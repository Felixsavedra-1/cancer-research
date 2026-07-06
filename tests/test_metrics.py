"""Known-answer tests for the benchmark ranking metrics (network-free)."""

import math

from cancer_tool import metrics


def test_roc_auc_perfect_separation():
    labels = [0, 0, 1, 1]
    scores = [0.1, 0.2, 0.8, 0.9]
    assert metrics.roc_auc(labels, scores) == 1.0


def test_roc_auc_inverted_is_zero():
    labels = [1, 1, 0, 0]
    scores = [0.1, 0.2, 0.8, 0.9]
    assert metrics.roc_auc(labels, scores) == 0.0


def test_roc_auc_random_is_half_with_ties():
    # All equal scores → every positive/negative pair is a tie → AUC 0.5.
    labels = [0, 1, 0, 1]
    scores = [0.5, 0.5, 0.5, 0.5]
    assert metrics.roc_auc(labels, scores) == 0.5


def test_roc_auc_known_value():
    # One positive outranks one of two negatives → AUC = (1 pair won)/(1*2) ... worked example.
    labels = [1, 0, 0]
    scores = [0.6, 0.9, 0.3]
    # positive (0.6) beats the 0.3 negative but loses to the 0.9 negative → 1/2.
    assert metrics.roc_auc(labels, scores) == 0.5


def test_roc_auc_single_class_is_nan():
    assert math.isnan(metrics.roc_auc([1, 1, 1], [0.1, 0.2, 0.3]))


def test_average_precision_perfect():
    labels = [1, 1, 0, 0]
    scores = [0.9, 0.8, 0.2, 0.1]
    assert metrics.average_precision(labels, scores) == 1.0


def test_average_precision_worked_example():
    # Ranking (desc score): pos, neg, pos, neg.
    # AP = (1/1)*1  at first pos (recall 0->0.5, precision 1/1)
    #    + (2/3)*1  at second pos (recall 0.5->1.0, precision 2/3)  ... averaged by delta-recall
    # AP = 0.5*(1.0) + 0.5*(2/3) = 0.8333...
    labels = [1, 0, 1, 0]
    scores = [0.9, 0.7, 0.5, 0.3]
    assert abs(metrics.average_precision(labels, scores) - (0.5 * 1.0 + 0.5 * (2 / 3))) < 1e-9


def test_average_precision_pessimistic_ties():
    # Positive and negative share the top score → positive ranked behind → precision 1/2.
    labels = [1, 0]
    scores = [0.5, 0.5]
    assert abs(metrics.average_precision(labels, scores) - 0.5) < 1e-9


def test_precision_and_recall_at_k():
    labels = [1, 0, 1, 0, 0]
    scores = [0.9, 0.8, 0.7, 0.6, 0.5]
    assert metrics.precision_at_k(labels, scores, 2) == 0.5  # top2 = {pos, neg}
    assert metrics.recall_at_k(labels, scores, 2) == 0.5     # 1 of 2 positives
    assert metrics.recall_at_k(labels, scores, 3) == 1.0     # both positives in top3


def test_prevalence():
    assert metrics.prevalence([1, 0, 0, 0]) == 0.25
