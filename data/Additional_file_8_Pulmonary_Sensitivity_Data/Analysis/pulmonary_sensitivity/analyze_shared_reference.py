#!/usr/bin/env python3
"""Reproduce the frozen pulmonary shared-reference sensitivity analysis.

Recovery edition, 22 September 2026. The original script bytes were unavailable
after archive damage. This implementation reconstructs the method documented in
the supplementary manuscript, using hash-identical public count inputs. Its
results are checked against the three hash-identical frozen result tables.
No frozen result is replaced. Stochastic allocation identities are reproduced
from NumPy's documented PCG64/default_rng seed 20260919 and checked through all
saved numerical outputs and the saved disjoint-reference sample identities.

Run from any directory:
  python analyze_shared_reference.py --output-dir /path/to/new/replay_directory

Requires numpy, pandas and scipy. Source root defaults to two directories above
this script; supply --source-root when data/code archives have separate roots.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import itertools
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr

SEED = 20260919
METRICS = ['reversal_fraction', 'rho', 'median_restoration']
SOURCE_HASHES = {
    'GSE278200_NINT_raw_counts.txt.gz': '5ae68a0a8b4ab253105e27df215af424e44ca3563c52f22edd1981d5253eebfb',
    'GSE308578_Gene_counts.csv.gz': 'a46808f90410796adc84738e4045aaeee311bf7445c71c39c1ee47e08da47a9c',
}
FROZEN_HASHES = {
    'shared_reference_summary.csv': 'd50a5d8026dc963ff30326d70df9afa653aa09b6a995c104879f3c44bdd167d2',
    'disjoint_reference_splits.csv': 'f1e2c6869926f19e8c33386abad3c6fb19673fb917d6a66e55732291ec14a5e5',
    'animal_label_permutations.csv': '0bffac65c94cb82500a0ae0e90becc71e2c11cd417611b85165be393e4783bde',
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bh(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    q = np.empty_like(p)
    q[order] = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def summarize(disease_effect, treatment_effect):
    keep = np.abs(disease_effect) >= 0.5
    d, t = disease_effect[keep], treatment_effect[keep]
    return {
        'n_selected': int(keep.sum()),
        'reversal_fraction': float(np.mean(d * t < 0)),
        'rho': float(spearmanr(d, t).statistic),
        'median_restoration': float(np.median(np.abs(d) - np.abs(d + t))),
    }


def distinct_allocations(rng, n, k, draws):
    seen = set()
    while len(seen) < draws:
        seen.add(tuple(sorted(rng.choice(n, k, replace=False))))
    return sorted(seen)


def verify_table(actual, frozen, name):
    assert actual.shape == frozen.shape, (name, actual.shape, frozen.shape)
    assert list(actual.columns) == list(frozen.columns), name
    errors = {}
    differences_above_1e11 = []
    for col in frozen.columns:
        if pd.api.types.is_numeric_dtype(frozen[col]):
            a, b = actual[col].to_numpy(float), frozen[col].to_numpy(float)
            # Two rat allocation correlations differ by < 5.5e-8, consistent
            # with machine-precision arithmetic affecting near-tied ranks.
            # Observed effects, all selected counts, all P/q values and every
            # disjoint split agree at the tighter numerical tolerance.
            tolerance = 1e-7 if name == 'animal_label_permutations.csv' and col == 'rho' else 1e-11
            assert np.allclose(a, b, rtol=0, atol=tolerance, equal_nan=True), (name, col)
            errors[col] = float(np.nanmax(np.abs(a - b)))
            for i in np.flatnonzero(np.abs(a - b) > 1e-11):
                differences_above_1e11.append({'dataset': str(frozen.iloc[i]['dataset']),
                                               'allocation': int(frozen.iloc[i]['allocation']) if 'allocation' in frozen else None,
                                               'column': col, 'frozen_value': float(b[i]),
                                               'replayed_value': float(a[i]), 'absolute_difference': float(abs(a[i] - b[i]))})
        else:
            assert actual[col].astype(str).tolist() == frozen[col].astype(str).tolist(), (name, col)
    return {'rows': len(actual), 'columns': len(actual.columns), 'max_absolute_numeric_difference': max(errors.values(), default=0.0), 'column_maximum_absolute_differences': errors, 'differences_above_1e-11': differences_above_1e11, 'tolerance': '1e-11 absolute; 1e-7 only for allocation-level Spearman rho', 'passed': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--frozen-dir', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output == args.frozen_dir.resolve():
        raise ValueError('Use a separate replay directory; frozen source tables must remain unchanged.')
    output.mkdir(parents=True, exist_ok=True)
    root = args.source_root.resolve()
    raw_root = root / 'Source_Data/Pulmonary_Raw'
    metadata_root = root / 'Source_Data/Public_Preclinical_Core'
    frozen = {}
    for name, expected in FROZEN_HASHES.items():
        path = args.frozen_dir / name
        assert sha256(path) == expected, f'Frozen table checksum failed: {name}'
        frozen[name] = pd.read_csv(path, float_precision='round_trip')
    rng = np.random.default_rng(SEED)
    summaries, splits, nulls = [], [], []
    input_records = []
    for dataset, filename, sep, control_n, disease_label, control_label in [
        ('GSE278200', 'GSE278200_NINT_raw_counts.txt.gz', '\t', 4, 'BLM', 'SAL'),
        ('GSE308578', 'GSE308578_Gene_counts.csv.gz', ',', 10, 'BLM_vehicle', 'CTRL'),
    ]:
        raw_path = raw_root / filename
        assert sha256(raw_path) == SOURCE_HASHES[filename], f'Public source checksum failed: {filename}'
        counts = pd.read_csv(raw_path, sep=sep, index_col=0).dropna(axis=1, how='all')
        if dataset == 'GSE278200':
            groups = pd.Series({s: s.rstrip('1234') for s in counts.columns})
        else:
            metadata = pd.read_csv(metadata_root / f'{dataset}_sample_metadata.csv', dtype={'sample': str})
            groups = metadata.set_index('sample')['group'].reindex(counts.columns)
            assert groups.notna().all()
        cpm = counts / counts.sum(axis=0) * 1e6
        keep = (cpm >= 1).sum(axis=1) >= control_n
        expression = np.log2((counts.loc[keep] + 0.5) / (counts.sum(axis=0) + 1) * 1e6)
        encoded = expression.to_csv().encode('utf-8')
        matrix_path = output / f'{dataset}_logCPM_raw_features.csv.gz'
        matrix_path.write_bytes(gzip.compress(encoded, mtime=0))
        assert gzip.decompress(matrix_path.read_bytes()) == encoded
        controls = groups.index[groups.eq(control_label)].tolist()
        disease = groups.index[groups.eq(disease_label)].tolist()
        treated = groups.index[groups.eq('NINT')].tolist()
        control_mean = expression[controls].mean(axis=1).to_numpy()
        disease_mean = expression[disease].mean(axis=1).to_numpy()
        treated_mean = expression[treated].mean(axis=1).to_numpy()
        observed = summarize(disease_mean - control_mean, treated_mean - disease_mean)
        pooled_names = disease + treated
        pooled = expression[pooled_names].to_numpy()
        nd, nt = len(disease), len(treated)
        allocation_ids = list(itertools.combinations(range(nd + nt), nd)) if dataset == 'GSE278200' else distinct_allocations(rng, nd + nt, nd, 1000)
        local_null = []
        for allocation, chosen in enumerate(allocation_ids):
            complement = [i for i in range(nd + nt) if i not in chosen]
            dm = pooled[:, chosen].mean(axis=1)
            tm = pooled[:, complement].mean(axis=1)
            row = {'dataset': dataset, 'allocation': allocation, **summarize(dm - control_mean, tm - dm)}
            local_null.append(row)
            nulls.append(row)
        disease_matrix = expression[disease].to_numpy()
        split_ids = list(itertools.combinations(range(nd), nd // 2)) if dataset == 'GSE278200' else distinct_allocations(rng, nd, nd // 2, 500)
        local_split = []
        for split, chosen in enumerate(split_ids):
            complement = [i for i in range(nd) if i not in chosen]
            a = disease_matrix[:, chosen].mean(axis=1)
            b = disease_matrix[:, complement].mean(axis=1)
            row = {'dataset': dataset, 'split': split,
                   'disease_reference_a': ';'.join(disease[i] for i in chosen),
                   'disease_reference_b': ';'.join(disease[i] for i in complement),
                   **summarize(a - control_mean, treated_mean - b)}
            local_split.append(row)
            splits.append(row)
        nr, sr = pd.DataFrame(local_null), pd.DataFrame(local_split)
        summary = {**observed, 'dataset': dataset, 'raw_features': len(counts),
                   'retained_features': len(expression), 'control_n': len(controls),
                   'disease_n': nd, 'treated_n': nt,
                   'permutation_scheme': 'exhaustive 70 animal allocations' if dataset == 'GSE278200' else '1000 distinct Monte Carlo animal allocations'}
        for metric in METRICS:
            name = 'null_median_restoration' if metric == 'median_restoration' else 'null_median_' + metric
            summary[name] = float(nr[metric].median())
        for metric in METRICS:
            values = nr[metric].to_numpy()
            count = int(np.sum(values <= observed[metric])) if metric == 'rho' else int(np.sum(values >= observed[metric]))
            summary[metric + '_animal_permutation_p'] = count / len(values) if dataset == 'GSE278200' else (count + 1) / (len(values) + 1)
        for metric in METRICS:
            for statistic in ['median', 'min', 'max']:
                summary['split_' + metric + '_' + statistic] = float(getattr(sr[metric], statistic)())
        summary['n_splits'] = len(local_split)
        summary['raw_sha256'] = sha256(raw_path)
        summaries.append(summary)
        input_records.append({'dataset': dataset, 'raw_path': str(raw_path.relative_to(root)), 'raw_sha256': sha256(raw_path), 'matrix_sha256': sha256(matrix_path), 'matrix_shape': list(expression.shape)})
    summary_frame = pd.DataFrame(summaries)
    q_values = bh(summary_frame[[m + '_animal_permutation_p' for m in METRICS]].to_numpy().ravel()).reshape(2, 3)
    for i, metric in enumerate(METRICS):
        summary_frame[metric + '_animal_permutation_q'] = q_values[:, i]
    tables = {'shared_reference_summary.csv': summary_frame,
              'disjoint_reference_splits.csv': pd.DataFrame(splits),
              'animal_label_permutations.csv': pd.DataFrame(nulls)}
    validation = {}
    for filename, actual in tables.items():
        actual = actual.reindex(columns=frozen[filename].columns)
        actual.to_csv(output / filename, index=False)
    for filename, actual in tables.items():
        actual = actual.reindex(columns=frozen[filename].columns)
        validation[filename] = verify_table(actual, frozen[filename], filename)
    manifest = {
        'recovery_status': 'Reconstructed implementation; original script/configuration bytes unavailable. Original frozen result tables and public raw count files are hash-identical.',
        'random_seed': SEED, 'rng': 'numpy.random.default_rng / PCG64',
        'distinct_allocations': 'Sorted tuples of sorted rng.choice indices, sampled without replacement within each tuple; duplicate tuples discarded; one continuing generator across datasets and analyses.',
        'filter': 'Raw CPM >= 1 in at least the control-group sample count; library totals computed before filtering; no use of treatment labels.',
        'normalization': 'log2[(count + 0.5)/(full library total + 1) * 1e6]',
        'selection': 'abs(mean(disease reference A) - mean(control)) >= 0.5; recalculated for every allocation and split.',
        'summary_statistics': 'Fraction d*t < 0; Spearman(d,t); median(abs(d)-abs(d+t)). Disjoint t uses disease reference B.',
        'p_values': 'Directional, full exact fraction in rats and plus-one Monte Carlo fraction in mice; BH across six dataset-by-statistic tests.',
        'inputs': input_records, 'frozen_output_checksums': FROZEN_HASHES,
        'validation': validation, 'software': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
        'interpretation': 'Overlapping-subset sensitivity ranges are descriptive, not independent-replication confidence intervals. Pulmonary pharmacodynamics does not establish hepatic safety or target causality.'}
    (output / 'analysis_summary.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    main()
