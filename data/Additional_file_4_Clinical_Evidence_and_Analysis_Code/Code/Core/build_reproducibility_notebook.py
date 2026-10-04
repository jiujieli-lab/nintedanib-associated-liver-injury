#!/usr/bin/env python3
"""Build the portable, executable audit notebook shipped with the submission.

The notebook is a release-capsule validator: it reads only files present in the
extracted submission directory. Heavy raw-data regeneration remains delegated
to the module scripts after third-party data retrieval.
"""

from pathlib import Path
import hashlib
import json


class _NotebookV4:
    """Small nbformat-v4 writer so the release builder has no extra dependency."""

    @staticmethod
    def new_notebook():
        return {"cells": [], "metadata": {}, "nbformat": 4, "nbformat_minor": 5}

    @staticmethod
    def new_markdown_cell(source):
        return {"cell_type": "markdown", "metadata": {}, "source": source}

    @staticmethod
    def new_code_cell(source):
        return {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": source,
        }


class _NotebookFormat:
    v4 = _NotebookV4()

    @staticmethod
    def validate(notebook):
        assert notebook.get("nbformat") == 4
        assert isinstance(notebook.get("cells"), list)
        ids = []
        for cell in notebook["cells"]:
            assert cell.get("cell_type") in {"markdown", "code"}
            assert isinstance(cell.get("source"), str)
            assert isinstance(cell.get("id"), str) and cell["id"]
            ids.append(cell["id"])
        assert len(ids) == len(set(ids))

    @staticmethod
    def write(notebook, path):
        Path(path).write_text(
            json.dumps(notebook, ensure_ascii=False, indent=1) + "\n",
            encoding="utf-8",
        )


nbf = _NotebookFormat()

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Public_Data_Reproducibility_Notebook.ipynb"

nb = nbf.v4.new_notebook()
nb["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.12"},
}

cells = []
cells.append(nbf.v4.new_markdown_cell("""# Role-separated public-data triangulation prioritizes target–phenotype hypotheses for nintedanib-associated liver injury

This notebook is the executable audit layer for a clean-slate study of fibrotic interstitial lung disease, nintedanib, disproportionate liver-injury reporting, human DILI proteomics, dose perturbation, liver-cell localization, virtual knockout, same-compartment evidence integration and public preclinical perturbation. It reads the frozen analysis products created by the dedicated scripts rather than an earlier manuscript or prior ranking.

**Interpretation guardrails.** FAERS estimates disproportionate reporting, not incidence or causality. The serum proteomics cohort contains multi-drug adjudicated DILI rather than nintedanib-specific cases, and its 12-protein panel was source-selected. LINCS is a HepG2 cell-line perturbation. scTenifold results are unsigned model-derived network displacement. GSE151374 has pooled libraries rather than independent-mouse single-cell inference. PXD052594 is an unconditioned pulmonary-response proteome, with a separate 5+5-animal phosphoproteome subset, rather than liver-toxicity evidence. The GSE299128 blood series contains no adjudicated DILI endpoint."""))

cells.append(nbf.v4.new_markdown_cell("""## tl;dr

- Nintedanib had disproportionate narrow and broad hepatic-event reporting relative to pirfenidone in the frozen openFDA release; this is a report-odds signal, not incidence.
- Among 12 source-selected proteins, FBP1 was the only FDR-significant separator of adjudicated DILI onset from non-drug acute liver injury in the reanalysed human serum dataset (Cliff's delta −0.505; q=4.69e-4), with higher values in non-drug injury; this is a phenotype anchor rather than an unbiased discovery claim.
- The six HepG2 dose signatures showed response structure, but no genome-wide concentration trend survived FDR, the prespecified DILI-anchor set was not enriched, and cross-cell-line concordance was weak.
- Human liver localization and three-seed virtual knockout define a testable target–phenotype map. FGFR1 was the only candidate meeting the prespecified exploratory T2 gate; the stricter all-seed sensitivity gate retained none, and no T1 target was assigned.
- GSE151374 pooled-library pseudobulk and PXD052594 unconditioned individual-animal proteomics provide model-boundary context for the prespecified target track while FBP1, GSTA1 and OTC remain strictly downstream phenotype anchors.
- No target is confirmed and no protective or harmful direction is inferred from the unsigned virtual-knockout model.
- Public preclinical perturbations add pulmonary-response and model-boundary constraints; pooled-slice and replicate-provenance-ambiguous datasets remain exploratory or descriptive rather than confirmatory."""))

cells.append(nbf.v4.new_code_cell("""from pathlib import Path
import hashlib, itertools, json, os, platform, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def locate_package_root(start):
    \"\"\"Find the extracted release without relying on its directory name.\"\"\"
    start = Path(start).resolve()
    candidates = []
    explicit = os.environ.get('NINTEDANIB_DILI_PACKAGE_ROOT')
    if explicit:
        candidates.append(Path(explicit).expanduser().resolve())
    candidates.extend([start, *start.parents])

    # Search at most three directory levels below the working directory. This
    # supports launching from a package parent or grandparent without an
    # unbounded filesystem crawl.
    frontier = [start]
    for _ in range(3):
        next_frontier = []
        for directory in frontier:
            try:
                children = [p for p in directory.iterdir() if p.is_dir()]
            except (OSError, PermissionError):
                children = []
            candidates.extend(children)
            next_frontier.extend(children)
        frontier = next_frontier
    for candidate in candidates:
        if ((candidate / 'Source_Data').is_dir()
                and (candidate / 'PACKAGE_VALIDATION.json').is_file()
                and (candidate / 'Figures').is_dir()):
            return candidate
    raise FileNotFoundError(
        'Could not locate the extracted package. Start Jupyter in the package '
        'directory or set NINTEDANIB_DILI_PACKAGE_ROOT to its absolute path.'
    )

ROOT = locate_package_root(Path.cwd())

class FrozenPackagePaths:
    \"\"\"Map frozen development paths used in audit cells to release paths.\"\"\"
    aliases = {
        'faers/results': 'Source_Data/FAERS',
        'proteomics': 'Source_Data/Human_DILI_Proteomics',
        'liver_singlecell': 'Source_Data/Liver_Single_Cell',
        'pharmacology': 'Source_Data/Pharmacology',
        'virtual_ko': 'Source_Data/Virtual_Knockout',
        'integration': 'Source_Data/Integration',
        'preclinical/results': 'Source_Data/Public_Preclinical_Core',
        'preclinical/gse151374': 'Source_Data/GSE151374',
        'preclinical/pxd052594': 'Source_Data/PXD052594',
        'blood_longitudinal': 'Source_Data/Whole_Blood',
        'figures': 'Figures',
    }

    def __init__(self, root):
        self.root = Path(root)

    def __truediv__(self, relative):
        relative = Path(relative).as_posix().strip('/')
        for prefix in sorted(self.aliases, key=len, reverse=True):
            if relative == prefix or relative.startswith(prefix + '/'):
                suffix = relative[len(prefix):].lstrip('/')
                return self.root / self.aliases[prefix] / suffix
        raise KeyError(f'No release-path mapping for: {relative}')

FA = FrozenPackagePaths(ROOT)

def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

release = json.loads((ROOT / 'PACKAGE_VALIDATION.json').read_text())
assert release['status'] == 'PASS'
print({
    'python': sys.version.split()[0],
    'platform': platform.platform(),
    'package_root': ROOT.name,
    'release_status': release['status'],
    'mode': 'portable frozen-output validation',
})"""))

cells.append(nbf.v4.new_markdown_cell("## Context & Methods\n\n### Source and inference-unit audit"))
cells.append(nbf.v4.new_code_cell("""sources = pd.DataFrame([
    ['openFDA/FAERS', 'report', 'nintedanib vs pirfenidone', 'spontaneous-report disproportionality'],
    ['DILI serum proteomics', 'source row', 'DO, NDO, HV, available follow-up blocks', 'multi-drug human injury phenotype'],
    ['LINCS GSE70138', 'Level-5 signature', 'six HepG2 doses at 24 h', 'observed dose perturbation'],
    ['GSE115469 liver atlas', 'donor', 'five healthy donors; cells descriptive', 'human liver localization'],
    ['ChEMBL', 'activity record/document', 'nintedanib single-protein assays', 'target pharmacology'],
    ['scTenifold', 'seed-specific inferred network', 'three cell compartments × three seeds', 'unsigned virtual-KO displacement'],
    ['Cross-domain integration', 'candidate within biological role', 'same-compartment localization and virtual-KO gate', 'hypothesis priority, not target confirmation'],
    ['GSE278200', 'animal', 'rat BLM lung; 4 per SAL, BLM and nintedanib group', 'independent-animal pulmonary-response context'],
    ['GSE308578', 'animal', 'mouse BLM lung; 10 control, 16 vehicle and 17 nintedanib', 'independent-animal model boundary'],
    ['GSE120804/GSE120679', 'pooled-tissue slice observation', 'rat liver/lung PCLS', 'exploratory context; not animal-population inference'],
    ['PXD024058', 'reported replicate column', 'mouse BLM lung TMT proteomics', 'descriptive; replicate provenance ambiguous'],
    ['GSE151374', 'pooled 10x library', 'three mice pooled/library; n=3 libraries per condition–time', 'exploratory pooled-library pseudobulk; not cell- or mouse-level inference'],
    ['PXD052594 proteome', 'individual mouse', 'all nintedanib n=10 vs all vehicle n=13', 'unconditioned pulmonary-response proteome; no radiomic response-cluster conditioning'],
    ['PXD052594 phosphoproteome', 'individual mouse in separate subset', 'nintedanib n=5 vs vehicle n=5', 'table-only pathway audit; not the full proteome cohort'],
    ['GSE299128', 'participant', 'seven longitudinally sampled participants', 'whole-blood sensitivity only'],
], columns=['source','inference_unit','design','permitted_interpretation'])
display(sources)"""))

cells.append(nbf.v4.new_markdown_cell("## Data and Results\n\n### Active-comparator pharmacovigilance"))
cells.append(nbf.v4.new_code_cell("""faers = pd.read_csv(FA/'faers/results/overall_active_comparator_signals.csv')
display(faers)
assert faers['endpoint'].isin(['narrow','broad','hy_proxy']).all()
assert (faers.filter(regex='count|reports', axis=1).select_dtypes('number') >= 0).all().all()
by_endpoint = faers.set_index('endpoint')
assert np.isclose(by_endpoint.loc['narrow','ROR'], 3.9966057509880484)
assert np.isclose(by_endpoint.loc['broad','ROR'], 3.3668703340517556)

heterogeneity = pd.read_csv(FA/'faers/results/year_effect_heterogeneity_tests.csv')
assert len(heterogeneity) == 8
heterogeneity_keyed = heterogeneity.set_index(['endpoint','window','test'])
narrow_bd = heterogeneity_keyed.loc[('narrow','all_years','Tarone-adjusted Breslow-Day')]
broad_q = heterogeneity_keyed.loc[('broad','all_years','Inverse-variance Cochran Q')]
assert np.isclose(narrow_bd.statistic, 88.25434621529466)
assert np.isclose(narrow_bd.p_value, 3.6637359812630166e-14)
assert int(narrow_bd.degrees_of_freedom) == 11
assert np.isclose(broad_q.statistic, 61.39384545688669)
assert np.isclose(broad_q.p_value, 1.256004261089914e-8)
assert np.isclose(broad_q.i2_percent, 80.45406683569463)

specifications = pd.read_csv(FA/'faers/results/report_level_specification_signals.csv')
suspect_analysis_ids = set(specifications.loc[
    specifications.analysis.str.contains('suspect_drug_report_filter'), 'analysis'
])
assert suspect_analysis_ids == {
    'suspect_drug_report_filter', 'ild_plus_suspect_drug_report_filter'
}
assert specifications.loc[
    specifications.analysis.str.contains('suspect_drug_report_filter'), 'filter_query'
].str.contains('drugcharacterization:\"1\"', regex=False).all()

roles = pd.read_csv(FA/'faers/results/drug_role_and_indication_profiles.csv')
role_lookup = roles.set_index(['endpoint','drug','variable','level'])['reports']
assert int(role_lookup.loc[('narrow','nintedanib','matching_drug_role','suspect')]) == 330
assert int(role_lookup.loc[('narrow','pirfenidone','matching_drug_role','suspect')]) == 105
assert int(role_lookup.loc[('broad','nintedanib','matching_drug_ild_indication','True')]) == 1373
assert int(role_lookup.loc[('broad','pirfenidone','matching_drug_ild_indication','True')]) == 604

display(heterogeneity)
print('FAERS audit passed: the primary estimand is all-indication drug-name-indexed report odds; suspect-role and ILD filters are sensitivity analyses.')
print('Official openFDA drugcharacterization codes distinguish suspect (1), concomitant (2), and interacting (3), but not primary versus secondary suspect.')"""))

cells.append(nbf.v4.new_markdown_cell("### Human DILI proteomic phenotype"))
cells.append(nbf.v4.new_code_cell("""eff = pd.read_csv(FA/'proteomics/effect_estimates.csv')
p = eff.query("cohort == 'confirmatory' and contrast == 'DO_vs_NDO'").copy()
cols = ['protein','n1_observed','n2_observed','cliffs_delta_group1_minus_group2',
        'cliffs_delta_ci_low','cliffs_delta_ci_high','mann_whitney_q_bh',
        'roc_auc_group1_positive','roc_auc_ci_low','roc_auc_ci_high']
display(p[cols].sort_values('mann_whitney_q_bh'))
fdr_hits = p.loc[p.mann_whitney_q_bh < 0.05, 'protein'].tolist()
assert fdr_hits == ['FBP1'], f'Unexpected FDR hit set: {fdr_hits}'
print('Confirmatory DO-vs-NDO FDR hit set:', fdr_hits)"""))

cells.append(nbf.v4.new_code_cell("""conc = pd.read_csv(FA/'proteomics/cross_cohort_concordance_summary.csv')
display(conc)
row = conc.query("contrast == 'DO_vs_NDO'").iloc[0]
assert int(row.n_direction_concordant) == 11 and int(row.n_proteins) == 12
print('Discovery-confirmatory directional concordance reproduced.')"""))

cells.append(nbf.v4.new_markdown_cell("### Nintedanib dose perturbation in HepG2"))
cells.append(nbf.v4.new_code_cell("""lincs_dir = ROOT/'Source_Data/LINCS'
lincs = pd.read_csv(lincs_dir/'gene_dose_response_metrics.csv.gz')
lincs_summary = json.loads((lincs_dir/'analysis_summary.json').read_text())
doses = lincs_summary['analysis_population']['doses_um']
assert doses == [0.04, 0.12, 0.37, 1.11, 3.33, 10.0]
assert len(lincs) == lincs_summary['analysis_population']['hepg2_level5_genes'] == 12328
assert int((lincs.bh_q_all_genes < 0.05).sum()) == lincs_summary['exact_dose_trend']['bh_q_lt_0_05'] == 0

anchors = ['FBP1','CES1','GSTA1','FAH','OTC','LECT2','ABCB11']
dose_columns = ['z_0p04uM','z_0p12uM','z_0p37uM','z_1p11uM','z_3p33uM','z_10p0uM']
anchor_view = lincs.loc[lincs.gene_symbol.isin(anchors),
    ['gene_symbol','measurement_stratum','spearman_rho','spearman_exact_p',
     'bh_q_all_genes','loo_rho_median','loo_sign_stability',*dose_columns]
].sort_values('spearman_rho', ascending=False)
display(anchor_view)

# Recompute each displayed exact six-dose Spearman statistic from packaged values.
for row in anchor_view.itertuples(index=False):
    values = np.asarray([getattr(row, c) for c in dose_columns], dtype=float)
    rho, _ = spearmanr(doses, values)
    assert np.isclose(rho, row.spearman_rho)
assert np.isclose(lincs_summary['dili_anchor_set']['permutation_p'], 0.5109395624175033)
print('Six-dose anchor trends reproduced from packaged values; no genome-wide trend passed FDR.')"""))

cells.append(nbf.v4.new_markdown_cell("### Human liver localization and target pharmacology"))
cells.append(nbf.v4.new_code_cell("""expr = pd.read_csv(FA/'liver_singlecell/target_expression_by_compartment.csv')
focus = ['FBP1','CES1','GSTA1','KDR','FLT1','FLT4','FGFR1','FGFR2','FGFR3','PDGFRA','PDGFRB']
display(expr[expr.gene_symbol.isin(focus)].sort_values(['gene_symbol','mean_expression'], ascending=[True,False]))
manifest = json.loads((FA/'liver_singlecell/analysis_manifest.json').read_text())
assert manifest['n_donors'] == 5 and manifest['n_cells'] == 8444
print('Liver atlas dimensions and donor count verified.')"""))

cells.append(nbf.v4.new_code_cell("""chem = pd.read_csv(FA/'pharmacology/chembl_target_summary.csv')
preferred = ['KDR','FLT1','FLT4','FGFR1','FGFR2','FGFR3','PDGFRA','PDGFRB']
display(chem[chem.gene_symbol.isin(preferred)].sort_values('max_pchembl', ascending=False))"""))

cells.append(nbf.v4.new_markdown_cell("### Virtual knockout stability"))
cells.append(nbf.v4.new_code_cell("""vko = pd.read_csv(FA/'virtual_ko/virtual_knockout_stability_summary.csv')
display(vko.sort_values('median_dili_anchor_percentile', ascending=False).head(25))
run = pd.read_csv(FA/'virtual_ko/virtual_knockout_run_manifest.csv')
seed_col = next(c for c in run.columns if 'seed' in c.lower())
assert run[seed_col].nunique() == 3
print('Three independent network seeds present. Values are unsigned displacement, not protection.')"""))

cells.append(nbf.v4.new_markdown_cell("### Same-compartment evidence integration"))
cells.append(nbf.v4.new_code_cell("""integration_dir = FA/'integration'
integration_manifest_path = integration_dir/'integration_analysis_manifest.json'
integration_manifest_observed = file_sha256(integration_manifest_path)
EXPECTED_INTEGRATION_MANIFEST_SHA256 = '00b4d15a4fd5ff65d1278bafd7e75be11297f2471b831ef141832c1ef4671688'
assert integration_manifest_observed == EXPECTED_INTEGRATION_MANIFEST_SHA256
integration_manifest = json.loads(integration_manifest_path.read_text())

figure5_files = {
    'svg': ROOT/'Figures/Figure_5_Multi_Evidence_Target_Prioritization.svg',
    'pdf': ROOT/'Figures/Figure_5_Multi_Evidence_Target_Prioritization.pdf',
    'png_600dpi': ROOT/'Figures/Figure_5_Multi_Evidence_Target_Prioritization.png',
    'tiff_600dpi': ROOT/'Figures/Figure_5_Multi_Evidence_Target_Prioritization.tiff',
}
for figure_format, figure_path in figure5_files.items():
    assert figure_path.is_file(), f'Missing Figure 5 artifact: {figure_path.name}'
    assert file_sha256(figure_path) == integration_manifest['figure_hashes']['figure5'][figure_format]

ranking = pd.read_csv(integration_dir/'candidate_consensus_ranking.csv')
tier2 = pd.read_csv(integration_dir/'tier2_candidate_shortlist.csv')
gate = pd.read_csv(integration_dir/'target_tier_gate_sensitivity.csv')
anchor_tiers = pd.read_csv(integration_dir/'phenotype_anchor_tiers.csv')

target_view = ranking.loc[
    ranking.candidate_class.eq('exposure_proximal_target'),
    ['gene_symbol','evidence_tier','consensus_rank','loso_rank_best','loso_rank_worst',
     'joint_context_compartment','joint_context_detection_fraction',
     'vko_within_compartment_rank','n_seeds_top_quartile','primary_T2_gate',
     'strict_all_seed_T2_sensitivity_gate','pharmacology_tier']
].sort_values('consensus_rank')
display(target_view.head(12))

assert set(ranking.candidate_class) == {'exposure_proximal_target','downstream_DILI_anchor'}
assert tier2.gene_symbol.tolist() == ['FGFR1']
assert tier2.candidate_class.tolist() == ['exposure_proximal_target']
fgfr1 = gate.set_index('gene_symbol').loc['FGFR1']
assert fgfr1.joint_context_compartment == 'Endothelial'
assert fgfr1.joint_context_detection_fraction >= 0.10
assert fgfr1.vko_within_compartment_rank >= 0.75
assert int(fgfr1.n_seeds) == 3 and bool(fgfr1.strong_pharmacology_gate)
assert bool(fgfr1.primary_T2_gate)
assert int(gate.strict_all_seed_T2_sensitivity_gate.sum()) == 0
assert set(gate.set_index('gene_symbol').loc[['FLT1','FGFR3'], 'evidence_tier']) == {
    'T3: screen/context target hypothesis'
}
assert anchor_tiers.loc[anchor_tiers.evidence_tier.str.startswith('P1:'), 'gene_symbol'].tolist() == ['FBP1']
anchor_roles = anchor_tiers.set_index('gene_symbol')
assert anchor_roles.loc[['FBP1','GSTA1','OTC'], 'candidate_class'].eq('downstream_DILI_anchor').all()
assert anchor_roles.loc['FBP1','evidence_tier'].startswith('P1:')
assert anchor_roles.loc[['GSTA1','OTC'], 'evidence_tier'].str.startswith('P2:').all()

integration_qa = json.loads((integration_dir/'integration_validation_report.json').read_text())
assert integration_qa['status'] == 'PASS'
assert integration_qa['n_checks'] == 44 == len(integration_qa['checks'])
assert all(check['passed'] for check in integration_qa['checks'])
print('Prespecified exploratory T2-gate set: FGFR1 only; strict all-seed sensitivity set: empty; no T1 target assigned.')
print('Downstream phenotype anchors: FBP1 (P1); GSTA1 and OTC (P2).')
print(f'Integration manifest SHA256: {integration_manifest_observed}; packaged Figure 5 hashes: 4/4 verified.')
print(f"Integration validation: {integration_qa['n_checks']}/{integration_qa['n_checks']} checks passed.")"""))

cells.append(nbf.v4.new_markdown_cell("""### Public preclinical pulmonary-response and tissue-context constraints

This section reads compact frozen products only. Independent-animal lung transcriptomes are kept separate from pooled-tissue slice observations and from the mouse-lung proteomic table whose replicate provenance is internally inconsistent. None of these models is treated as evidence of nintedanib-specific human DILI."""))
cells.append(nbf.v4.new_code_cell("""pre = FA/'preclinical/results'
pre_summary = json.loads((pre/'analysis_summary.json').read_text())
pre_scope = pd.DataFrame([
    ['GSE278200','Rattus norvegicus','lung RNA-seq','independent animal',
     pre_summary['GSE278200']['group_n']],
    ['GSE308578','Mus musculus','lung RNA-seq','independent animal',
     pre_summary['GSE308578']['group_n']],
    ['GSE120804','Rattus norvegicus','liver PCLS RNA-seq','pooled-tissue slice observation; exploratory',
     pre_summary['GSE120804']['group_n']],
    ['GSE120679','Rattus norvegicus','lung PCLS RNA-seq','pooled-tissue slice observation; exploratory',
     pre_summary['GSE120679']['group_n']],
    ['PXD024058','Mus musculus','lung TMT proteomics','reported columns; descriptive because provenance is ambiguous',
     {'proteins': pre_summary['PXD024058']['n_proteins']}],
], columns=['accession','species','assay','inference_unit','group_counts'])
display(pre_scope)

assert pre_summary['GSE278200']['group_n'] == {'BLM':4,'NINT':4,'SAL':4}
assert pre_summary['GSE308578']['group_n'] == {'BLM_vehicle':16,'CTRL':10,'NINT':17}
assert 'pooled-tissue' in pre_summary['GSE120804']['scope']
assert 'pooled-tissue' in pre_summary['GSE120679']['scope']
assert 'descriptive only' in pre_summary['PXD024058']['scope']
assert 'No sample was relabelled' in pre_summary['GSE308578']['metadata_resolution']

core_manifest_path = pre/'analysis_manifest.json'
core_manifest_observed = file_sha256(core_manifest_path)
EXPECTED_CORE_MANIFEST_SHA256 = '95e2f087db45a9ad2a09ea1bbddbf4333e4fb5c4db6eaf6e7761697002f555f9'
assert core_manifest_observed == EXPECTED_CORE_MANIFEST_SHA256
pre_manifest = json.loads(core_manifest_path.read_text())
pre_qa = json.loads((pre/'validation_report.json').read_text())
assert pre_qa['validation_passed'] and pre_qa['n_failed'] == 0
assert pre_qa['n_checks'] == 8 == len(pre_qa['checks'])
assert all(check['passed'] for check in pre_qa['checks'])
manifest_sha = {entry['path']: entry['sha256'] for entry in pre_manifest['files']}
selected_preclinical_files = [
    'GSE278200_candidate_results.csv','GSE308578_candidate_results.csv',
    'GSE278200_module_contrasts.csv','GSE308578_module_contrasts.csv',
    'GSE308578_leave_one_animal_out_summary.csv','GSE120804_candidate_results.csv',
    'GSE120679_candidate_results.csv','PXD024058_candidate_effects_descriptive.csv',
    'cross_model_candidate_module_effects.csv','analysis_summary.json','validation_report.json'
]
for filename in selected_preclinical_files:
    rel = f'fresh_analysis/preclinical/results/{filename}'
    assert rel in manifest_sha, f'Manifest entry missing: {rel}'
    assert file_sha256(pre/filename) == manifest_sha[rel], f'Hash mismatch: {filename}'
print(f'Frozen preclinical hashes verified: {len(selected_preclinical_files)}/{len(selected_preclinical_files)} files.')"""))

cells.append(nbf.v4.new_code_cell("""rat_lung = pd.read_csv(pre/'GSE278200_candidate_results.csv')
mouse_lung = pd.read_csv(pre/'GSE308578_candidate_results.csv')
rat_treatment = rat_lung.query("contrast == 'NINT_vs_BLM'")
mouse_treatment = mouse_lung.query("contrast == 'NINT_vs_BLM_vehicle'")
independent_animal_candidates = pd.concat([
    rat_treatment.assign(model='GSE278200 rat lung, n=4 vs 4'),
    mouse_treatment.assign(model='GSE308578 mouse lung, n=17 vs 16')
], ignore_index=True)
candidate_view = independent_animal_candidates.loc[
    independent_animal_candidates.gene.isin(['FGFR1','FGFR3','FLT1','FBP1','GSTA1','OTC','CES1D','COL1A1','POSTN','ABCB1B']),
    ['model','gene','detected_after_filtering','log2FC','CI95_low','CI95_high','q_value_BH_genome']
]
display(candidate_view)

for table in [rat_treatment, mouse_treatment]:
    observed_q = table.q_value_BH_genome.dropna()
    assert observed_q.between(0, 1).all()
assert rat_treatment.set_index('gene').loc['COL1A1','q_value_BH_genome'] < 0.05
assert rat_treatment.set_index('gene').loc['POSTN','q_value_BH_genome'] < 0.05
assert mouse_treatment.set_index('gene').loc['ABCB1B','q_value_BH_genome'] < 0.05

rat_modules = pd.read_csv(pre/'GSE278200_module_contrasts.csv').query("contrast == 'NINT_vs_BLM'")
mouse_modules = pd.read_csv(pre/'GSE308578_module_contrasts.csv').query("contrast == 'NINT_vs_BLM_vehicle'")
module_view = pd.concat([
    rat_modules.assign(model='GSE278200 rat lung'),
    mouse_modules.assign(model='GSE308578 mouse lung')
], ignore_index=True)[['model','module','effect_z_score','p_permutation','q_BH_within_contrast_modules','experimental_unit_scope']]
display(module_view)
assert module_view.experimental_unit_scope.eq('independent animal').all()
assert module_view.p_permutation.between(0, 1).all()
assert module_view.q_BH_within_contrast_modules.between(0, 1).all()

loo = pd.read_csv(pre/'GSE308578_leave_one_animal_out_summary.csv').set_index('feature')
assert loo.loc['ABCB1B','loo_min'] > 0 and loo.loc['ABCB1B','sign_concordance'] == 1.0
print('Independent-animal readout: pulmonary responses vary by model; no liver-toxicity inference is made.')

cross = pd.read_csv(pre/'cross_model_candidate_module_effects.csv')
cross_view = (cross[cross.feature.isin(['FGFR1','FLT1','FBP1','OTC','GSTA1','ECM/fibrosis'])]
              .pivot(index='model', columns='feature', values='nintedanib_vs_disease_log2_effect'))
display(cross_view)"""))

cells.append(nbf.v4.new_code_cell("""liver_slices = pd.read_csv(pre/'GSE120804_candidate_results.csv').query("contrast == 'BDL_Nint_vs_BDL'")
lung_slices = pd.read_csv(pre/'GSE120679_candidate_results.csv').query("contrast == 'TGFb_Nint_vs_TGFb'")
mouse_protein = pd.read_csv(pre/'PXD024058_candidate_effects_descriptive.csv')

for slice_table in [liver_slices, lung_slices]:
    detected_slices = slice_table.loc[slice_table.detected_after_filtering]
    assert len(detected_slices) > 0 and detected_slices.inferential_scope.notna().all()
    assert detected_slices.inferential_scope.str.contains('pooled before randomisation').all()
assert mouse_protein.inferential_scope.str.contains('descriptive only').all()
assert not any(column.startswith(('p_','q_')) for column in mouse_protein.columns)

context_view = pd.concat([
    liver_slices[liver_slices.gene.isin(['FGFR1','FLT1','FBP1','OTC'])]
        .assign(source='GSE120804 liver PCLS')[['source','gene','log2FC','inferential_scope']],
    lung_slices[lung_slices.gene.isin(['FGFR1','FLT1','FBP1','OTC'])]
        .assign(source='GSE120679 lung PCLS')[['source','gene','log2FC','inferential_scope']],
    mouse_protein[mouse_protein.gene.isin(['FBP1','OTC','GSTA1'])]
        .rename(columns={'NIB_vs_BLM_log2FC_descriptive':'log2FC'})
        .assign(source='PXD024058 mouse lung protein')[['source','gene','log2FC','inferential_scope']]
], ignore_index=True)
display(context_view)
print('Pooled-slice statistics remain exploratory; ambiguous proteomic replicate columns remain descriptive.')"""))

cells.append(nbf.v4.new_markdown_cell("""### GSE151374 pooled-library single-cell perturbation extension

The experimental unit is the pooled 10x library: three mice were pooled per library and there are three libraries per condition–time cell. The notebook therefore audits library-level pseudobulk estimates and exact 3-versus-3 permutation tests only; it does not treat cells or constituent mice as independent replicates. Raw H5 files are deliberately not opened here."""))
cells.append(nbf.v4.new_code_cell("""s6 = FA/'preclinical/gse151374'
s6_manifest_path = s6/'results/analysis_manifest.json'
s6_manifest_observed = file_sha256(s6_manifest_path)
EXPECTED_S6_MANIFEST_SHA256 = 'e558f28bdc62172eabf64f5c9e0e969c1a8cc9517d03afdd1947d3da08a749b8'
assert s6_manifest_observed == EXPECTED_S6_MANIFEST_SHA256

s6_manifest = json.loads(s6_manifest_path.read_text())
s6_audit = json.loads((s6/'validation/manifest_hash_validation.json').read_text())
s6_checks = json.loads((s6/'validation/validation_checks.json').read_text())
assert s6_manifest['statistical_unit'] == 'pooled 10x library (three mice pooled per library; n=3 libraries per condition/time point)'
assert s6_audit['manifest_sha256'] == s6_manifest_observed
assert s6_audit['entries_checked'] == 43 and s6_audit['all_recorded_hashes_match']
assert s6_audit['release_entries_checked'] == 23
assert s6_audit['omitted_public_raw_inputs_checked_upstream'] == 20
assert s6_audit['all_available_release_hashes_match']
for entry in s6_manifest['inputs'] + s6_manifest['outputs']:
    if entry['available_in_release']:
        released_file = ROOT/entry['release_path']
        assert released_file.is_file()
        assert released_file.stat().st_size == entry['bytes']
        assert file_sha256(released_file) == entry['sha256']
assert s6_checks == {
    'sample_count': 18,
    'all_six_condition_cells_have_three_libraries': True,
    'all_libraries_retain_at_least_1000_cells': True,
    'all_libraries_retain_at_least_100_macrophages': True,
    'no_cell_level_inferential_tests': True,
    'exact_permutation_minimum_two_sided_p_for_3v3': 0.1,
    'gene_level_BH_present_day7': True,
    'gene_level_BH_present_day14': True,
    'module_BH_present': True,
}

s6_libraries = pd.read_csv(s6/'results/library_manifest.csv')
library_cells = s6_libraries.groupby(['day','group']).size()
assert len(s6_libraries) == 18 and len(library_cells) == 6
assert library_cells.eq(3).all()
assert s6_libraries.n_mice_pooled.eq(3).all()
assert s6_libraries.statistical_unit.eq('pooled 10x library').all()
display(library_cells.rename('n_pooled_libraries').reset_index())
print(f'GSE151374 detached manifest audit: 43/43 PASS; SHA256 {s6_manifest_observed}.')"""))

cells.append(nbf.v4.new_code_cell("""s6_candidates = pd.read_csv(s6/'results/candidate_gene_contrasts.csv')
s6_expected_roles = {
    'Fgfr1': 'exposure_proximal_target',
    'Fbp1': 'DILI_phenotype_anchor_P1',
    'Gsta1': 'DILI_phenotype_anchor_P2',
    'Otc': 'DILI_phenotype_anchor_P2',
}
for gene, expected_role in s6_expected_roles.items():
    gene_rows = s6_candidates.loc[s6_candidates.gene.eq(gene)]
    assert len(gene_rows) == 4
    assert gene_rows.track_type.eq(expected_role).all()
    assert gene_rows.n_pool_each_group.eq(3).all()

assert s6_candidates.p_exact_permutation.ge(0.1).all()
assert s6_candidates.q_BH_exact_all_candidate_tests.between(0, 1).all()
s6_modules = pd.read_csv(s6/'results/macrophage_module_contrasts.csv')
assert len(s6_modules) == 10 and s6_modules.n_pool_each_group.eq(3).all()
assert s6_modules.p_exact_permutation.ge(0.1).all()
assert s6_modules.q_BH_exact_10_tests.ge(0.75).all()

s6_key = s6_candidates.loc[
    s6_candidates.gene.isin(s6_expected_roles),
    ['compartment','track_type','gene','day','log2CPM_difference_NINT_minus_BLM',
     'CI95_low','CI95_high','p_welch','p_exact_permutation',
     'q_BH_welch_all_candidate_tests','q_BH_exact_all_candidate_tests','n_pool_each_group']
].sort_values(['compartment','gene','day'])
display(s6_key)

s6_fgfr1_mac_d14 = s6_candidates.query(
    "compartment == 'Macrophage stratum' and gene == 'Fgfr1' and day == 14"
).iloc[0]
assert np.isclose(s6_fgfr1_mac_d14.log2CPM_difference_NINT_minus_BLM, 1.1859432853372582)
assert np.isclose(s6_fgfr1_mac_d14.p_exact_permutation, 0.1)
assert np.isclose(s6_fgfr1_mac_d14.q_BH_exact_all_candidate_tests, 1.0)
print('S6 role audit: FGFR1 is exposure-proximal; FBP1/GSTA1/OTC are downstream phenotype anchors.')
print('At n=3 pooled libraries per arm, the minimum attainable two-sided exact P is 0.1; no exact-test FDR claim is made.')"""))

cells.append(nbf.v4.new_markdown_cell("""### PXD052594 unconditioned individual-animal proteome and table-only phosphoproteome

The primary deposited lung proteome comparison includes every available animal in the processed matrix (nintedanib n=10 versus vehicle n=13) and is not conditioned on the post-treatment radiomic response cluster. The phosphoproteome is a distinct random 5-versus-5 subset and is audited from its compact summary table only. The original workbook and full phosphosite table are deliberately not opened here."""))
cells.append(nbf.v4.new_code_cell("""s7 = FA/'preclinical/pxd052594'
s7_manifest_path = s7/'results/analysis_manifest.json'
s7_manifest_observed = file_sha256(s7_manifest_path)
EXPECTED_S7_MANIFEST_SHA256 = '9ef9f38c22ad0da706cd9df455188257e215c82b6ebb76550d0b6a5d2c6d7647'
assert s7_manifest_observed == EXPECTED_S7_MANIFEST_SHA256

s7_manifest = json.loads(s7_manifest_path.read_text())
s7_audit = json.loads((s7/'validation/manifest_hash_validation.json').read_text())
s7_checks = json.loads((s7/'validation/validation_checks.json').read_text())
assert s7_manifest['statistical_unit'] == 'individual mouse'
assert s7_manifest['contrast'] == 'all nintedanib (n=10) versus all vehicle animals in processed proteome matrix (n=13)'
assert s7_manifest['selection_guardrail'] == 'No radiomic response cluster was used for conditioning, exclusion, stratification, or adjustment.'
assert s7_audit['manifest_sha256'] == s7_manifest_observed
assert s7_audit['entries_checked'] == 20 and s7_audit['all_recorded_hashes_match']
assert s7_audit['release_entries_checked'] == 17
assert s7_audit['omitted_public_raw_inputs_checked_upstream'] == 3
assert s7_audit['all_available_release_hashes_match']
for entry in s7_manifest['inputs'] + s7_manifest['outputs']:
    if entry['available_in_release']:
        released_file = ROOT/entry['release_path']
        assert released_file.is_file()
        assert released_file.stat().st_size == entry['bytes']
        assert file_sha256(released_file) == entry['sha256']

assert s7_checks['proteins'] == 7006 and s7_checks['animals'] == 23
assert s7_checks['nintedanib_animals'] == 10 and s7_checks['vehicle_animals'] == 13
assert s7_checks['matrix_has_no_missing_values'] and s7_checks['animal_ids_unique']
assert not s7_checks['post_treatment_cluster_used_in_analysis']
assert s7_checks['all_remaining_animals_included'] and s7_checks['BH_applied_across_all_proteins']
assert s7_checks['phosphosites'] == 20043
assert s7_checks['phosphoproteome_nintedanib_animals'] == 5
assert s7_checks['phosphoproteome_vehicle_animals'] == 5
assert s7_checks['phosphoproteome_has_no_missing_values'] and s7_checks['phosphosite_BH_applied']

s7_animals = pd.read_csv(s7/'results/animal_sample_manifest.csv')
assert len(s7_animals) == 23 and s7_animals.animal_id.nunique() == 23
assert s7_animals.statistical_unit.eq('individual mouse').all()
assert s7_animals.treatment.value_counts().to_dict() == {'vehicle':13, 'nintedanib':10}
s7_phospho_samples = pd.read_csv(s7/'results/phosphoproteome_PCA_and_sample_manifest.csv')
assert len(s7_phospho_samples) == 10 and s7_phospho_samples.animal_id.nunique() == 10
assert s7_phospho_samples.treatment.value_counts().to_dict() == {'nintedanib':5, 'vehicle':5}
display(pd.DataFrame([
    ['proteome','individual mouse',10,13,'all deposited processed-matrix animals; unconditioned'],
    ['phosphoproteome','individual mouse in separate subset',5,5,'compact table-only pathway audit'],
], columns=['assay','inference_unit','n_nintedanib','n_vehicle','analysis_scope']))
print(f'PXD052594 detached manifest audit: 20/20 PASS; SHA256 {s7_manifest_observed}.')"""))

cells.append(nbf.v4.new_code_cell("""s7_candidates = pd.read_csv(s7/'results/candidate_contrasts.csv')
s7_expected_roles = {
    'Fgfr1': 'exposure_proximal_target',
    'Fbp1': 'DILI_phenotype_anchor_P1',
    'Gsta1': 'DILI_phenotype_anchor_P2',
    'Otc': 'DILI_phenotype_anchor_P2',
}
for gene, expected_role in s7_expected_roles.items():
    row = s7_candidates.set_index('gene').loc[gene]
    assert row.track_type == expected_role
    assert int(row.n_NINT) == 10 and int(row.n_vehicle) == 13
assert s7_candidates.loc[s7_candidates.track_type.eq('exposure_proximal_target'), 'gene'].tolist() == ['Fgfr1']
assert bool(s7_candidates.set_index('gene').loc['Fgfr1','detected'])
assert not s7_candidates.set_index('gene').loc[['Fbp1','Gsta1','Otc'],'detected'].astype(bool).any()

s7_candidate_view = s7_candidates.loc[
    s7_candidates.gene.isin(s7_expected_roles),
    ['gene','track_type','detected','log2FC_NINT_vs_vehicle','CI95_low','CI95_high',
     'p_welch','p_mannwhitney','q_BH_welch_detected_tracks',
     'q_BH_mannwhitney_detected_tracks','n_NINT','n_vehicle']
]
display(s7_candidate_view)

s7_modules = pd.read_csv(s7/'results/module_contrasts.csv')
assert len(s7_modules) == 5
assert s7_modules.n_NINT.eq(10).all() and s7_modules.n_vehicle.eq(13).all()
for column in ['p_welch','p_mannwhitney','p_permutation_100k',
               'q_BH_p_welch','q_BH_p_mannwhitney','q_BH_p_permutation_100k']:
    assert s7_modules[column].between(0, 1).all()

s7_phospho_targets = pd.read_csv(s7/'results/phosphosite_target_pathway_summary.csv')
assert int(s7_phospho_targets.n_sites_q_lt_0_05.sum()) == 0
s7_fgfr1_phospho = s7_phospho_targets.set_index('gene').loc['Fgfr1']
assert int(s7_fgfr1_phospho.n_sites_quantified) == 0
display(s7_phospho_targets)
print('S7 role audit: FGFR1 is the sole prespecified exposure-proximal target row; FBP1/GSTA1/OTC remain phenotype anchors and were not quantified.')
print('No audited target/pathway phosphosite had BH q<0.05; the separate 5+5 subset is table-only context, not liver-toxicity proof.')"""))

cells.append(nbf.v4.new_markdown_cell("### Whole-blood negative-control/sensitivity layer"))
cells.append(nbf.v4.new_code_cell("""blood = pd.read_csv(FA/'blood_longitudinal/target_longitudinal_summary.csv')
blood_slopes = pd.read_csv(FA/'blood_longitudinal/participant_slopes.csv')
display(blood.sort_values('exact_signflip_p').head(15))

def exact_signflip_p(values):
    values = np.asarray(values, dtype=float)
    observed = abs(float(values.mean()))
    null = [abs(float(np.mean(values * np.asarray(signs))))
            for signs in itertools.product([-1.0, 1.0], repeat=len(values))]
    return sum(value >= observed - 1e-15 for value in null) / len(null)

recomputed = (blood_slopes.groupby('gene_symbol')['slope_per_visit']
              .apply(exact_signflip_p).rename('recomputed_exact_p'))
blood_check = blood.set_index('gene_symbol').join(recomputed)
assert np.allclose(blood_check.exact_signflip_p, blood_check.recomputed_exact_p)
assert (blood.n_participants == 7).all()
assert np.allclose(blood.exact_signflip_p * 128, np.round(blood.exact_signflip_p * 128))
assert (blood.q_bh >= 0.05).all()
assert np.isclose(blood.set_index('gene_symbol').loc['KDR','exact_signflip_p'], 0.0625)
assert np.isclose(blood.set_index('gene_symbol').loc['FGFR1','exact_signflip_p'], 0.421875)
assert np.isclose(blood.set_index('gene_symbol').loc['FBP1','exact_signflip_p'], 0.5625)
print('Exact 2^7 sign-flip P values reproduced: KDR=0.0625, FGFR1=0.421875, FBP1=0.5625.')
print('No prespecified blood target passed FDR; this is not a DILI validation cohort.')"""))

cells.append(nbf.v4.new_markdown_cell("""## Takeaways

1. Nintedanib has disproportionate liver-injury reporting relative to pirfenidone in the prespecified openFDA analysis; the estimand remains a report-odds signal.
2. Among the 12 source-selected proteins, FBP1 is the only FDR-significant discriminator of adjudicated DILI onset from alternative acute liver injury in the reanalysed human serum dataset, with higher values in non-drug injury. It is therefore a phenotype anchor, not an unbiased discovery or a proven susceptibility target.
3. Nintedanib produces dose-structured HepG2 transcriptional responses in several liver-injury anchors, but cross-cell-line concordance is weak and high doses are not treated as clinical exposures.
4. Role-separated integration did not splice expression and virtual-knockout evidence across compartments. FGFR1 alone met the prespecified exploratory cross-seed-median T2 gate in endothelium. JAK1 met the all-seed rank component but lacked strong-pharmacology support, so no target met the strict all-seed sensitivity gate and no T1 target was assigned; this is prioritization, not target confirmation.
5. Independent-animal rat and mouse lung transcriptomes define model-dependent pulmonary-response constraints. Rat lung COL1A1 and POSTN responses and mouse lung ABCB1B induction are directly auditable above; they do not establish liver toxicity or demonstrate preserved clinical efficacy.
6. GSE151374 contributes an exploratory pooled-library pseudobulk perturbation layer (n=3 libraries per condition–time, three mice pooled per library). The smallest possible two-sided 3-versus-3 exact P is 0.1, so nominal Welch results are not treated as confirmatory.
7. PXD052594 contributes an unconditioned individual-animal lung proteome (n=10 versus 13). FGFR1 is the sole prespecified exposure-proximal target row; FBP1, GSTA1 and OTC remain downstream phenotype anchors and were not quantified. The distinct 5-versus-5 phosphoproteome supports table-only pathway auditing and had no target/pathway site with BH q<0.05.
8. Pooled rat liver/lung slice experiments and the replicate-provenance-ambiguous PXD024058 mouse lung proteome provide only exploratory or descriptive context. Their anchor responses prevent us from treating circulating anchors as liver-exclusive intervention targets.
9. Exact participant-level sign-flip tests in GSE299128 reproduced on the full 2^7 assignment space; all 30 prespecified genes had q=1.0. The series remains a blood detectability sensitivity layer, not a DILI cohort.
10. The next confirmatory step is bidirectional perturbation with parent-drug/metabolite measurement and simultaneous retention of prespecified antifibrotic response readouts, as prespecified in the experimental validation protocol."""))

# nbformat 4.5 requires a unique cell identifier.  Derive deterministic IDs
# from cell order, type, and source so rebuilt notebooks remain warning-free
# and produce stable JSON when the analytical content is unchanged.
for index, cell in enumerate(cells, start=1):
    material = f"{index}\0{cell['cell_type']}\0{cell['source']}".encode("utf-8")
    cell["id"] = f"cell-{index:02d}-{hashlib.sha256(material).hexdigest()[:12]}"

nb["cells"] = cells
nbf.validate(nb)
nbf.write(nb, OUT)
print(OUT)
