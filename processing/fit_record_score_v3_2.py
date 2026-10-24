"""Fit v3.2a/b offline from declared inputs, keeping assay data private.

Set TCRVDB_PATH or pass --tcrvdb-path. The selected-records CSV (or test.csv
member of a ZIP) uses anchor-inclusive junctions. No source paths or identities
are written into model parameters. Inputs are read anew on each invocation.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import time
from pathlib import Path
from zipfile import ZipFile

TCRVDB_PATH = os.environ.get('TCRVDB_PATH')


def labels(truth, selected=None):
    import polars as pl
    required = {'tcr_id', 'pos', 'ambiguous'}
    if not required <= set(truth.columns):
        raise ValueError('Validation CSV requires tcr_id, pos, ambiguous')
    if (truth['tcr_id'].null_count() or truth['tcr_id'].n_unique() != truth.height or
            truth.filter(pl.col('tcr_id').cast(pl.String).str.strip_chars() == '').height):
        raise ValueError('Validation IDs must be nonmissing and unique')
    if truth.filter(~pl.col('pos').is_in([0, 1]) | ~pl.col('ambiguous').is_in([0, 1]) |
                    pl.col('pos').is_null() | pl.col('ambiguous').is_null()).height:
        raise ValueError('pos and ambiguous must be binary')
    if truth.filter(pl.col('pos') + pl.col('ambiguous') > 1).height:
        raise ValueError('Positive and ambiguous flags must be mutually exclusive')
    if 'neg' in truth.columns:
        if truth.filter(~pl.col('neg').is_in([0, 1]) | pl.col('neg').is_null() |
                        (pl.col('pos') + pl.col('neg') + pl.col('ambiguous') != 1)).height:
            raise ValueError('Outcome flags must be exactly one of positive, negative, ambiguous')
    if selected is not None:
        fields = [n for n in ['epitope_aa', 'hla_long', 'cdr3_alpha_aa', 'cdr3_beta_aa',
                  'full_seq_alpha_aa', 'full_seq_beta_aa', 'TRAV_IMGT', 'TRAJ_IMGT', 'TRBV_IMGT', 'TRBJ_IMGT']
                  if n in truth.columns and n in selected.columns]
        matched = truth.join(selected.select('tcr_id', *fields), on='tcr_id',
                             how='left', suffix='_selected', validate='1:1')
        for field in fields:
            if matched.filter(pl.col(field).cast(pl.String).fill_null('') !=
                              pl.col(field + '_selected').cast(pl.String).fill_null('')).height:
                raise ValueError(f'Validation construct conflicts with selected input: {field}')
    return truth.select('tcr_id', pl.when(pl.col('ambiguous') == 0).then(pl.col('pos'))
                        .otherwise(-1).alias('label'))


def queries(test):
    import polars as pl
    from vdjdb.score.evidence import canonical
    if (test['tcr_id'].null_count() or test['tcr_id'].n_unique() != test.height or
            test.filter(pl.col('tcr_id').str.strip_chars() == '').height):
        raise ValueError('Selected IDs must be nonmissing and unique')
    return canonical(test.select(
        'tcr_id', 'epitope_aa', pl.lit('HomoSapiens').alias('species'),
        pl.col('epitope_aa').alias('antigen.epitope'),
        pl.concat_str(pl.lit('HLA-'), pl.col('hla_long')).alias('mhc.a'),
        pl.lit('B2M').alias('mhc.b'), pl.lit('MHCI').alias('mhc.class'),
        pl.col('cdr3_alpha_aa').alias('a_junction'), pl.col('cdr3_beta_aa').alias('b_junction'),
        *(pl.col(col).str.replace_all('_DV', '/DV', literal=True).alias(name) for col, name in [
            ('TRAV_IMGT', 'a_v'), ('TRAJ_IMGT', 'a_j'), ('TRBV_IMGT', 'b_v'), ('TRBJ_IMGT', 'b_j')])))


def retained(records, manifest):
    """Use the author's exact manifest, never redraw or alter motif geometry."""
    import polars as pl
    population = records.filter((pl.col('antigen.epitope') == 'SLLMWITQV') &
        (pl.col('reference.id') == 'PMID:40498839') &
        (pl.col('method.identification') == 'phage display, magnetic beads'))
    if manifest['record_id'].n_unique() != manifest.height:
        raise ValueError('Selection manifest IDs must be unique')
    selected = population.join(manifest.select('record_id'), on='record_id', how='semi')
    if selected.height != manifest.height or manifest.height != 90:
        raise ValueError('Selection manifest must select exactly 90 source records')
    if selected.group_by('meta.study.id').len().sort('meta.study.id')['len'].to_list() != [30, 30, 30]:
        raise ValueError('Selection manifest must retain 30 records per group')
    remove = population.join(manifest.select('record_id'), on='record_id', how='anti').select('record_id')
    return records.join(remove, on='record_id', how='anti')


def targets(records, chains):
    import polars as pl
    counted = chains.filter(pl.col('cdr3').fill_null('') != '').select('record_id', 'clonotype_id').join(
        records.select('record_id', 'pmhc_id', 'reference.id'), on='record_id', validate='m:1')
    counted = counted.filter(pl.col('reference.id').fill_null('') != '').with_columns(
        pl.col('reference.id').n_unique().over('clonotype_id', 'pmhc_id').alias('studies'))
    y = counted.group_by('record_id').agg((pl.col('studies') > 1).any().cast(pl.Int64).alias('label'))
    return records.select('record_id').join(y, on='record_id', how='left', validate='1:1').with_columns(
        pl.col('label').fill_null(0))


def prepare(tables, selected_records, tcrvdb_path, selection_manifest, motif_scores):
    import polars as pl
    from vdjdb.emit.vdjdb3 import read_table
    from vdjdb.score.evidence import FEATURES, KEY, interactions
    from vdjdb.score.v32 import collapse_confidence, fixed_features
    records, chains, evidence = [read_table(tables, n) for n in ['records', 'chains', 'evidence']]
    if selected_records.suffix == '.zip':
        with ZipFile(selected_records) as archive:
            test = pl.read_csv(io.BytesIO(archive.read('test.csv')), infer_schema_length=0)
    else:
        test = pl.read_csv(selected_records, infer_schema_length=0)
    q = queries(test).with_row_index('_row')
    truth = labels(pl.read_csv(tcrvdb_path), test)
    if truth.join(q.select('tcr_id'), on='tcr_id', how='anti').height:
        raise ValueError('Validation IDs absent from the selected cohort')
    q = q.join(truth, on='tcr_id', how='left', validate='1:1').with_columns(pl.col('label').fill_null(-1))
    full = interactions(records, chains)
    units = pl.concat([full.select(KEY), q.select(KEY)]).unique()
    features = fixed_features(records, chains, evidence, pl.read_parquet(motif_scores), units)
    selected = q.join(features.select('interaction_key', *FEATURES),
                      on='interaction_key', how='left', validate='m:1')
    source = collapse_confidence(test.select('tcr_id',
        pl.col('vdjdb_score').cast(pl.Float64).alias('legacy_score')))
    selected = selected.join(source.rename({'legacy_score': '_source_score'}),
        on='tcr_id', validate='1:1').with_columns(
            pl.coalesce('legacy_score', '_source_score').alias('legacy_score')).sort('_row')
    kept = retained(records, pl.read_csv(selection_manifest, separator='\t'))
    database = full.join(kept.select('record_id'), on='record_id', how='semi').join(
        features.select('interaction_key', *FEATURES), on='interaction_key', validate='m:1')
    prior = database.join(q.select('interaction_key').unique(), on='interaction_key', how='anti')
    return selected, prior, database


def metrics(scored):
    import numpy as np
    import polars as pl
    from sklearn.metrics import average_precision_score, auc, precision_recall_curve, roc_auc_score
    rows = []
    for (epitope,), group in scored.group_by('epitope', maintain_order=True):
        valid = group.filter(pl.col('label') >= 0)
        y, p = valid['label'].to_numpy(), valid['prediction'].to_numpy()
        row = dict(epitope=epitope, n=len(y), positive=int(y.sum()), negative=int(len(y)-y.sum()),
                   unknown=group.height-valid.height, ap=None, pr_auc=None, roc_auc=None)
        if len(np.unique(y)) == 2:
            precision, recall, _ = precision_recall_curve(y, p)
            row.update(ap=float(average_precision_score(y, p)), pr_auc=float(auc(recall, precision)),
                       roc_auc=float(roc_auc_score(y, p)))
        rows.append(row)
    return pl.DataFrame(rows)


def fit_assay(args):
    import numpy as np
    import pandas as pd
    import polars as pl
    from record_score_learning import fit_network, learn_structure
    from vdjdb.score.evidence import FEATURES
    from vdjdb.score.learning import predict_network
    from vdjdb.score.pipeline import load_model
    from vdjdb.score.v32 import encode
    selected, prior, _ = prepare(args.tables, args.selected_records, args.tcrvdb_path,
                                 args.selection_manifest, args.motif_scores)
    template, _ = load_model(version='3.2')
    x, u = [pd.DataFrame(encode(template, frame).to_dict(as_series=False)) for frame in (selected, prior)]
    y = selected['label'].to_numpy()
    states = {n: sorted(set(x[n]) | set(u[n]) | {'?'}) for n in FEATURES}
    # Reproduce the accepted bins/graphs by default. Fresh structure searches are explicit.
    for arm in args.arms:
        for method in args.methods:
            def fit_context(xx, yy, graph_template, extra_x=None):
                graph = graph_template['edges'] if method == 'frozen' else learn_structure(xx, yy, method=method)
                train_x, train_y = xx, yy
                if arm == 'labeled':
                    train_x, train_y = xx.loc[yy >= 0].reset_index(drop=True), yy[yy >= 0]
                elif arm == 'all-vdjdb':
                    train_x = pd.concat([xx, extra_x], ignore_index=True)
                    train_y = np.r_[yy, np.full(len(extra_x), -1)]
                model = fit_network(train_x, train_y, graph, states, max_iter=1000)
                if not model['converged']:
                    raise ValueError('Model fit did not converge')
                return model

            model = fit_context(x, y, template, u)
            local_models = {}
            for epitope, local_template in template['epitope_models'].items():
                if epitope == 'GILGFVFTL' and arm == 'all-vdjdb':
                    continue  # This requested alternative uses the full-corpus global model.
                mask = selected['epitope_aa'].to_numpy() == epitope
                if not mask.any():
                    continue
                if epitope == 'GILGFVFTL':
                    # Accepted GIL context uses all labeled global observations.
                    local_x, local_y = x.loc[y >= 0].reset_index(drop=True), y[y >= 0]
                    graph = local_template['edges'] if method == 'frozen' else learn_structure(local_x, local_y, method=method)
                    local = fit_network(local_x, local_y, graph, states, max_iter=1000)
                else:
                    local = fit_context(x.loc[mask].reset_index(drop=True), y[mask], local_template,
                                        u.loc[prior['antigen.epitope'].to_numpy() == epitope])
                if not local['converged']:
                    raise ValueError('Local model fit did not converge')
                local_models[epitope] = local
            model.update(epitope_models=local_models, motif_bounds=template['motif_bounds'],
                         iptm_bounds=template['iptm_bounds'], membership_mode='producer',
                         membership_width=template['membership_width'], score_version='3.2',
                         validation_score_version='3.2a', decision_threshold=.5)
            name = arm + '--' + method
            prediction = predict_network(model, x, epitopes=selected['epitope_aa'].to_list())
            write_fit(args.out, name, model, selected.select('tcr_id', pl.col('epitope_aa').alias('epitope'), 'label'),
                      prediction)


def write_fit(out, name, model, rows, prediction):
    import numpy as np
    import polars as pl
    if not np.isfinite(prediction).all():
        raise ValueError('Nonfinite predictions')
    (out / (name + '.json')).write_text(json.dumps(model, indent=2) + '\n')
    scored = rows.with_columns(pl.Series('prediction', prediction))
    scored.write_parquet(out / (name + '.parquet'), compression='zstd')
    metrics(scored).write_csv(out / (name + '-metrics.tsv'), separator='\t')


def fit_independent(args):
    import pandas as pd
    import polars as pl
    from record_score_learning import fit_network, learn_structure
    from vdjdb.emit.vdjdb3 import read_table
    from vdjdb.score.learning import predict_network
    from vdjdb.score.pipeline import load_model
    from vdjdb.score.v32 import encode, independent_features
    r, c, e = [read_table(args.tables, n) for n in ['records', 'chains', 'evidence']]
    r = retained(r, pl.read_csv(args.selection_manifest, separator='\t'))
    c, e = [f.join(r.select('record_id'), on='record_id', how='semi') for f in (c, e)]
    f = independent_features(r, c, e).join(targets(r, c), on='record_id', validate='1:1').sort('record_id')
    template, _ = load_model(version='3.2b')
    x = pd.DataFrame(encode(template, f).to_dict(as_series=False))
    states = {n: sorted(set(x[n]) | {'?'}) for n in x.columns}
    y = f['label'].to_numpy()
    for method in args.methods:
        graph = template['edges'] if method == 'frozen' else learn_structure(x, y, method=method)
        model = fit_network(x, y, graph, states)
        if not model['converged']:
            raise ValueError('Independent-study model did not converge')
        model.update(score_version='3.2b', iptm_bounds=template['iptm_bounds'])
        write_fit(args.out, 'independent--' + method, model,
                  f.select('record_id', pl.col('antigen.epitope').alias('epitope'), 'label'),
                  predict_network(model, x))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables', type=Path, required=True)
    parser.add_argument('--tcrvdb-path', type=Path, default=TCRVDB_PATH,
                        help='Private validation CSV; defaults to TCRVDB_PATH')
    parser.add_argument('--selected-records', type=Path)
    parser.add_argument('--selection-manifest', type=Path, required=True)
    parser.add_argument('--motif-scores', type=Path, help='Fixed continuous producer scores; not rebuilt')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', choices=['assay', 'independent', 'both'], default='assay')
    parser.add_argument('--methods', nargs='+', choices=['frozen', 'hc-bic', 'hc-aic', 'aracne'], default=['frozen'])
    parser.add_argument('--arms', nargs='+', choices=['labeled', 'selected', 'all-vdjdb'], default=['selected'])
    args = parser.parse_args()
    if args.task != 'independent' and not args.tcrvdb_path:
        parser.error('Set TCRVDB_PATH or pass --tcrvdb-path')
    inputs = [args.selection_manifest]
    if args.task != 'independent':
        if not args.selected_records or not args.motif_scores:
            parser.error('Assay fitting requires --selected-records and --motif-scores')
        inputs += [args.tcrvdb_path, args.selected_records, args.motif_scores]
    inputs += [args.tables / (n + '.parquet') for n in ['records', 'chains', 'evidence']]
    for path in inputs:
        if not path.is_file():
            parser.error(f'Required input file absent: {path.name}')
    if args.out.exists():
        parser.error('Use a new output destination; completed results are never overwritten')
    args.out.mkdir(parents=True)
    from threadpoolctl import threadpool_limits
    started = time.perf_counter()
    with threadpool_limits(limits=1):
        if args.task != 'independent':
            fit_assay(args)
        if args.task != 'assay':
            fit_independent(args)
    # Private provenance goes beside predictions, never inside distributed parameters.
    provenance = {'input_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                  'seconds': time.perf_counter()-started, 'no_folds': True,
                  'methods': args.methods, 'arms': args.arms}
    (args.out / 'SOURCES.json').write_text(json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    main()
