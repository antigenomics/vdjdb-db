"""Cluster metrics, vendored from TCREMP with exactly one deliberate divergence.

⚠ **The one divergence is the tie-break in :func:`binominal_test`, and it is a correctness fix.**
Upstream assigns each cluster its majority label by sorting on ``fraction_matched`` alone and keeping
the first row, so a cluster split evenly between two epitopes takes whichever the incoming row order
put first -- and it sorts on ``p_value`` first, a constant 0.0 when ``compute_pvalue=False``, so
pandas' unstable quicksort permutes the frame before the tie is reached. ``precision``, ``recall``,
``f1`` and ``ami`` are all computed against that label, so all four depend on row order, while
``purity`` does not because it sums counts a tie leaves alone. Measured on the human TRB cohort with
the do-nothing partition, changing nothing but the row order: precision 0.9458, 0.9409, 0.9407,
0.9458. Two CI runs of the same build disagreed at 0.9459 against 0.9407, which is how it was found.

So the original intent below -- "bit-identical semantics" -- **was never achievable**: an order-
dependent number cannot be reproduced without reproducing the order, and nothing does. The tie-break
is now explicit (``fraction_matched``, then ``count_matched``, then the label name) and every affected
value **rose** by 0.0001 to 0.001, because preferring the larger matched count can never pick a
minority label where a majority exists. Filed upstream as antigenomics/tcremp#8. See the comment at
the divergence and ``ROADMAP_local.md`` §53.

Everything else is copied unchanged from `2026-vdjdb-update/review/rev2_benchmark/scripts/metrics_lib.py`, which is
itself a verbatim copy of the TCREMP source. Phase 11's acceptance criterion is a comparison against
numbers in that benchmark's `results/metrics_full.tsv`, and a comparison is only a comparison if both
sides are scored the same way.

⚠ **Not interchangeable with `mir.bench.metrics.cluster_metrics`**, which uses
`recall = tp / n_true_clustered` over clustered records only, where
:func:`precision_recall_fscore` here folds unclustered records into FN. The two metric families are
not comparable and mixing them silently rescales every number (ROADMAP section 8.9).

This module is the one place in the package that uses pandas; it is a validation dependency
(the ``motifs`` extra), never a build dependency.
"""
# Vendored from the TCREMP source so the benchmark is self-contained and env-independent, with the
# single tie-break divergence the module docstring names. Provenance:
#   binominal_test, count_clstr_purity  <- /Users/mikesh/vcs/code/tcremp/tcremp/ml_utils.py
#   precision_recall_fscore, get_clustermetrics <- /Users/mikesh/vcs/code/tcremp/tcremp/metrics.py
# These define purity / retention / consistency / AMI / precision / recall / F1 as used for the
# TCREMP paper (Kremlyakova et al. 2025, Table 1). Do not "improve" beyond the tie-break already
# fixed: every compared method has to be scored by the same code, and a second change would make
# our numbers incomparable with the benchmark for no stated reason.
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import multilabel_confusion_matrix, adjusted_mutual_info_score


def binominal_test(df, cluster, group, threshold=0.7, compute_pvalue=True):
    binom_df = df.copy()
    binom_df['total_cluster'] = binom_df.groupby(cluster)[cluster].transform('count')
    binom_df['total_group'] = binom_df.groupby(group)[group].transform('count')
    binom_df['count_matched'] = binom_df.groupby([group, cluster])[group].transform('count')
    binom_df['fraction_matched'] = binom_df['count_matched'] / binom_df['total_cluster']
    binom_df['fraction_matched_exp'] = binom_df['total_group'] / len(binom_df.index)
    if compute_pvalue:
        binom_df['p_value'] = binom_df.apply(
            lambda row: stats.binomtest(row['count_matched'], n=row['total_cluster'], p=row['fraction_matched_exp'],
                                        alternative='greater').pvalue, axis=1)
    else:
        binom_df['p_value'] = 0.0  # p_value is stored but not used by get_clustermetrics; skip for speed
    binom_df_cluster = binom_df[
        [group, cluster, 'total_cluster', 'total_group', 'count_matched', 'fraction_matched', 'fraction_matched_exp',
         'p_value']].drop_duplicates()
    binom_df_cluster['is_cluster'] = binom_df_cluster.apply(
        lambda x: 1 if (x.total_cluster > 1) and (x.cluster != -1) else 0, axis=1)
    binom_df_cluster['enriched_clstr'] = binom_df_cluster.apply(lambda x: 1
    if (x.fraction_matched >= threshold)
       and (x.is_cluster == 1) else 0, axis=1)
    # ⚠ DIVERGENCE FROM THE VENDORED ORIGINAL, and the only one: the tie-break is explicit.
    #
    # These two lines assign each cluster its majority `group` (the `label_cluster` that
    # `precision`, `recall`, `f1` and `ami` are all computed against). Upstream sorts on
    # `fraction_matched` alone and keeps the first row, so **a cluster split evenly between two
    # epitopes takes whichever one the incoming row order happened to put first** - and the original
    # also sorted on `p_value` first, which is a constant 0.0 whenever `compute_pvalue=False`, so
    # pandas' default unstable quicksort permuted the frame before the tie was even reached.
    #
    # Measured on the human TRB cohort with the do-nothing partition, shuffling the input rows and
    # changing nothing else: precision 0.9458, 0.9409, 0.9407, 0.9458 over four orders, with purity
    # identical at 0.9345 throughout, because purity sums `count_matched` over the kept rows and a
    # tie leaves those sums alone. Two CI runs of the same build disagreed at 0.9459 against 0.9407
    # for exactly this reason (ROADMAP_local section 53), which CLAUDE.md hard rule 7 forbids: the
    # same inputs must give the same bytes on another host at another core count. Filed upstream as
    # antigenomics/tcremp#8; this divergence goes away when that lands.
    #
    # So: highest `fraction_matched`, then the larger `count_matched`, then the `group` name
    # ascending. Every term is a property of the data, none is a property of the order, and the
    # third is only ever reached by a cluster tied on both counts, where any choice is arbitrary and
    # what matters is that it is the *same* arbitrary choice every time. `kind="stable"` so pandas
    # cannot reintroduce order dependence between equal keys.
    binom_df_cluster = binom_df_cluster.sort_values(
        ['fraction_matched', 'count_matched', group], ascending=[False, False, True], kind='stable')
    binom_df_cluster = binom_df_cluster.drop_duplicates('cluster', keep='first')
    return binom_df_cluster


def count_clstr_purity(binom_res):
    binom_res_clstr = binom_res[binom_res['is_cluster'] == 1]
    if len(binom_res_clstr) != 0:
        return sum(binom_res_clstr['count_matched']) / sum(binom_res_clstr['total_cluster'])


def precision_recall_fscore(df, ytrue, ypred, label):
    labels = np.unique(ytrue)
    precisions, recalls, accuracies, weights, supports = [], [], [], [], []
    fn_per_epi = df[df['cluster'].isnull()][label].value_counts()
    epmetrics = {'accuracy': {}, 'precision': {}, 'recall': {}, 'f1-score': {}, 'support': {}}
    for (i, cm) in enumerate(multilabel_confusion_matrix(ytrue, ypred, labels=labels)):
        tn = cm[0][0]; fn = cm[1][0]; tp = cm[1][1]; fp = cm[0][1]
        lbl = labels[i]
        missing_fn = fn_per_epi.get(lbl, 0)
        fn += missing_fn
        precision = 0.0 if tp + fp == 0 else tp / (tp + fp)
        recall = 0.0 if tp + fn == 0 else tp / (tp + fn)
        accuracy = 0 if tp + tn == 0 else (tp + tn) / (tp + tn + fp + fn)
        support = sum(ytrue == lbl)
        w = support / ytrue.shape[0]
        weights.append(w)
        precision *= w; recall *= w; accuracy *= w
        accuracies.append(accuracy); precisions.append(precision); recalls.append(recall)
        supports.append(support)
        epmetrics['accuracy'][labels[i]] = accuracy
        epmetrics['precision'][labels[i]] = precision
        epmetrics['recall'][labels[i]] = recall
        f = 0 if (precision * recall == 0) or (precision + recall == 0) else 2 * (precision * recall) / (precision + recall)
        epmetrics['f1-score'][labels[i]] = f
        epmetrics['support'][labels[i]] = support
    uncalled_epis = set(df[label]).difference(labels)
    for i in uncalled_epis:
        supports.append(sum(ytrue == i)); recalls.append(0); precisions.append(0)
        epmetrics['accuracy'][i] = 0; epmetrics['precision'][i] = 0
        epmetrics['recall'][i] = 0; epmetrics['f1-score'][i] = 0; epmetrics['support'][i] = sum(ytrue == i)
    recall = sum(recalls); precision = sum(precisions); support = sum(supports)
    f = 0 if (precision * recall == 0) or (precision + recall == 0) else 2 * (precision * recall) / (precision + recall)
    return 0, precision, recall, f, support, epmetrics


def get_clustermetrics(data_df, label):
    binom_res = data_df[['cluster', 'label_cluster', 'total_cluster', 'count_matched', 'fraction_matched', 'p_value',
                         'fraction_matched_exp', 'is_cluster']].drop_duplicates()
    df = data_df[data_df['is_cluster'] == 1]
    purity = count_clstr_purity(binom_res)
    ypred = df['label_cluster']
    ytrue = df[label]
    ami = adjusted_mutual_info_score(ytrue, ypred)
    accuracy, precision, recall, f1score, support, epmetrics = precision_recall_fscore(df, ytrue, ypred, label)
    counts = {k: v for k, v in df[label].value_counts().reset_index().values.tolist()}
    maincluster = {ep: df[df[label] == ep]['cluster'].value_counts().index[0] for ep in df[label].unique()}
    consistencymap = {ep: len(df[(df[label] == ep) & (df['cluster'] == maincluster[ep])]) / counts[ep] for ep in counts.keys()}
    return {
        'purity': round(purity, 4),
        'retention': round(len(data_df[data_df['is_cluster'] == 1]) / len(data_df), 4),
        'consistency': round(np.mean([(consistencymap[ep] * counts[ep]) / len(df) for ep in consistencymap.keys()]), 4),
        'ami': round(ami, 4),
        'precision': round(precision, 4),
        'recall': round(recall, 4),
        'f1': round(f1score, 4),
        'mean_clustsize': round(np.mean(list(binom_res[binom_res['is_cluster'] == 1]['total_cluster'])), 2),
    }


def metrics_from_assignments(assign_df, label='antigen.epitope', compute_pvalue=False):
    """assign_df must have columns [label, 'cluster'] (one row per record).
    Reproduces benchmark/models.py wiring (binominal_test -> merge -> get_clustermetrics).
    p_value is not used by get_clustermetrics, so compute_pvalue defaults False (fast)."""
    data = assign_df.copy()
    binom_res = binominal_test(data, 'cluster', label, compute_pvalue=compute_pvalue).rename({label: 'label_cluster'}, axis=1)
    data = data.merge(binom_res, on='cluster', how='left')
    data['is_cluster'] = data['is_cluster'].fillna(0)
    if (data['is_cluster'] == 1).sum() == 0:            # degenerate: nothing clustered
        m = {'purity': 0.0, 'retention': 0.0, 'consistency': 0.0, 'ami': 0.0,
             'precision': 0.0, 'recall': 0.0, 'f1': 0.0, 'mean_clustsize': 0.0}
    else:
        m = get_clustermetrics(data, label)
    m['n_records'] = len(data)
    m['n_epitopes'] = data[label].nunique()
    m['n_clusters_kept'] = int((binom_res['is_cluster'] == 1).sum())
    return m, binom_res, data


def per_epitope_metrics(data, label='antigen.epitope'):
    """Clean one-vs-rest per-epitope precision/recall/F1/retention.
    prediction = cluster majority epitope (label_cluster); FN includes unclustered records."""
    rows = []
    clustered = data[data['is_cluster'] == 1]
    for ep in sorted(data[label].unique()):
        n_ep = int((data[label] == ep).sum())
        retention = float((data[label] == ep).pipe(lambda s: (data.loc[s.index, 'is_cluster'] == 1).mean()))
        tp = int(((clustered[label] == ep) & (clustered['label_cluster'] == ep)).sum())
        fp = int(((clustered[label] != ep) & (clustered['label_cluster'] == ep)).sum())
        fn_clustered = int(((clustered[label] == ep) & (clustered['label_cluster'] != ep)).sum())
        fn_unclustered = int(((data[label] == ep) & (data['is_cluster'] != 1)).sum())
        fn = fn_clustered + fn_unclustered
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        rows.append({label: ep, 'n_records': n_ep, 'retention': round(retention, 4),
                     'precision': round(prec, 4), 'recall': round(rec, 4), 'f1': round(f1, 4)})
    return pd.DataFrame(rows)
