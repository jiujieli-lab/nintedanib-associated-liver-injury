#!/usr/bin/env python3
import json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
p=Path(__file__).resolve().parent
s=pd.read_csv(p/'candidate_network_summary.csv');g=pd.read_csv(p/'network_specification_scores.csv');l=pd.read_csv(p/'network_donor_lodo_scores.csv');n=pd.read_csv(p/'network_null_matching_sensitivity.csv');a=pd.read_csv(p/'network_shrinkage_ablation.csv');c=pd.read_csv(p/'network_target_anchor_contributions.csv');m=pd.read_csv(p/'network_algorithm_comparison_metrics.csv')
checks={
 'fixed_17_candidate_universe_retained':s.gene_symbol.nunique()==17 and len(s)==51,
 'three_compartments_retained':s.compartment.nunique()==3,
 'fifty_primary_target_compartment_scores':s.network_available.sum()==50,
 'eighteen_specs_by_four_clinical_weights':len(g)==3600 and g.specification.nunique()==18 and g.weight_profile.nunique()==4,
 'five_real_donor_omissions':len(l)==250 and l.omitted_donor.nunique()==5,
 'all_primary_numerically_finite':s.loc[s.network_available,'matched_null_z'].notna().all(),
 'maximum_linear_solve_residual_below_1e_10':g.linear_solve_max_abs_residual.max()<1e-10,
 'tail_probabilities_respect_null_resolution':g.empirical_p_enrichment.between(.0005,1).all(),
 'q_values_between_zero_and_one':g.empirical_q_bh.between(0,1).all(),
 'all_main_scores_nonnegative':g.diffusion_score.ge(0).all(),
 'sum_anchor_contributions_equals_primary_score':np.allclose(c.groupby(['compartment','gene_symbol']).weighted_contribution.sum().sort_index(),s[s.network_available].set_index(['compartment','gene_symbol']).diffusion_score.sort_index(),rtol=1e-10,atol=1e-14),
 'no_detection_supported_significant_stable_target':int(s.expression_supported_stability_gate.sum())==0,
}
report={'decision':'PASS' if all(checks.values()) else 'FAIL','checks':{k:bool(v) for k,v in checks.items()},'n_checks':len(checks),'n_lodo_zero_null_variance':int(l.matched_null_z.isna().sum()),'runtime_numeric_precision':'64-bit floating point','warning':'A successful numerical validation is not validation of a liver-toxicity mechanism.'}
(p/'network_validation.json').write_text(json.dumps(report,indent=2))
comparisons=[]
for comp,x in a[a.weight_profile.eq('discovery_DO_vs_HV')].groupby('compartment'):
 r,pv=spearmanr(x.matched_null_z_unshrunk,x.matched_null_z_shrinkage)
 comparisons.append({'compartment':comp,'comparison':'unshrunk_Pearson_vs_LedoitWolf','spearman_rho':r,'n_candidates':len(x),'max_rank_difference':(x.candidate_rank_unshrunk-x.candidate_rank_shrinkage).abs().max()})
for comp,x in g[g.is_primary_spec].groupby('compartment'):
 w=x.pivot(index='gene_symbol',columns='weight_profile',values='matched_null_z')
 for label in w.columns:
  if label=='discovery_DO_vs_HV':continue
  r,pv=spearmanr(w['discovery_DO_vs_HV'],w[label]);comparisons.append({'compartment':comp,'comparison':label+'_vs_discovery_DO_vs_HV','spearman_rho':r,'n_candidates':len(w),'max_rank_difference':np.nan})
pd.DataFrame(comparisons).to_csv(p/'network_additional_benchmark_metrics.csv',index=False)
(p/'Methods_Results_Network_Diffusion.md').write_text('''# Clinically weighted, donor-blocked network propagation

## Methods

### Linking clinical injury proteins to drug-responsive liver networks

We implemented a clinically weighted network-propagation analysis to test whether the 17 pharmacologically supported candidate genes were unusually proximal to the 13-protein human DILI signature in liver reference networks. This analysis directly connects the circulating injury phenotype to candidate drug-binding proteins through observed cell-compartment expression structure. The target universe was fixed from the existing pharmacology table before fitting the new algorithm: FLT4, PDGFRA, FGFR1, KDR, JAK1, LYN, FLT1, FGFR2, FGFR3, ABL1, CDK4, FLT3, FGFR4, TYK2, AXL, LRRK2 and MET. All 17 candidates were retained in the output, including an explicit unavailable entry when a candidate was absent from a compartment network. No candidate-specific threshold was used and no source was removed for discordance with a preferred target.

We used the publicly deposited normalized GSE115469 expression values already archived in the reproducibility capsule. The matrix orientation was verified as genes by cells and cell identifiers were matched exactly to the deposited donor and cell annotations. Three reference compartments contained 705 genes in 3,501 hepatocytes, 717 genes in 844 endothelial cells and 716 genes in 1,192 macrophages. Each compartment contained cells from five biological donors. The supplied feature universe comprised the original 650 variance-ranked genes meeting at least 1% cell detection, augmented with the clinical, pharmacological and mechanism genes, after excluding genes with zero variance in that compartment. Thus, the reanalysis is conditional on this fixed, previously selected gene universe; leave-one-donor-out analysis does not repeat the upstream feature-universe selection. We did not reinterpret deposited normalized values as raw counts or repeat count normalization.

### Donor-aware covariance estimation and graph construction

For each compartment and donor separately, we centered and standardized each gene across that donor's cells using the population standard deviation. A gene constant within that donor was represented by a zero centered vector. We fitted Ledoit–Wolf shrinkage covariance to this matrix with an assumed zero mean, converted the covariance to a correlation matrix, and averaged the five donor matrices with equal coefficients. This procedure prevents a donor with thousands of cells from entering as thousands of independent biological replicates. Donor-specific shrinkage nevertheless regularizes estimates according to their information content, and an unshrunk Pearson-correlation ablation was fitted to measure its effect.

The primary graph retained positive correlation weights and the 20 strongest neighbors of each gene. An edge was retained if either endpoint selected the other; the resulting weighted graph was symmetric, self-edges were removed, and isolated nodes received a self-loop only to define a stochastic transition matrix. The fixed sensitivity grid combined positive versus absolute-correlation edges, 10, 20 or 40 neighbors, and restart probabilities of 0.3, 0.5 or 0.7, yielding 18 specifications. The primary setting was positive edges, 20 neighbors and a restart probability of 0.5. Absolute-correlation edges provide a sensitivity to connectivity independent of correlation sign; they do not confer an activating or inhibitory interpretation. Specifications were enumerated without optimizing target identity or clinical outcome performance.

### Clinical weighting and propagation statistic

The 13 seed proteins were ACO1, ASS1, FAH, CPS1, ALDOB, HPD, OTC, DMGDH, GSTA1, FBP1, PCK2, CES1 and LECT2. For the primary analysis, the seed weight for protein j was the absolute discovery-cohort Cliff's delta for DILI onset versus healthy volunteers, divided by the sum of the 13 absolute effect sizes. Cliff's delta was chosen instead of a fold change so that the weighting was comparable across the source assay scales. All 13 proteins contributed continuously; no significance threshold or target preference was used to select seeds. Three additional profiles used discovery DILI onset versus non-DILI liver injury onset, confirmatory DILI onset versus healthy volunteers, and confirmatory DILI onset versus non-DILI liver injury onset. These profiles evaluated dependence on clinical weighting and were not used to choose a winning parameter setting. Because weights were unsigned, the algorithm tests connectivity to an injury-associated protein pattern and does not determine whether target inhibition is protective or harmful.

Let W be the symmetric adjacency matrix, d its column sums, P = W diag(d)^(-1) the column-stochastic transition matrix, and v the normalized seed vector. The propagated score was f = r [I − (1−r)P]^(-1)v. Candidate-specific scores were evaluated by solving the equivalent adjoint linear systems in double precision. The maximum absolute residual was required to be below 10^(-10), column-stochasticity below 10^(-12), and calculated scores nonnegative. The contribution of an individual seed was retained as its seed weight multiplied by the corresponding propagation coefficient; these 13 contributions sum exactly to the candidate score. The method is a random walk with restart, building on established network-propagation methods, rather than a newly invented form of diffusion.

### Expression-matched nulls and graph-degree sensitivity

For each compartment, background genes excluded the 13 seeds and the 17 target candidates. Each seed was matched to its 50 nearest background genes using Euclidean distance in standardized log(1 + equal-donor mean expression) and logit detection, with a 0.005 pseudocount applied symmetrically to detection and non-detection. We drew 1,999 null seed sets with random seed 20260919. Within a null set, a matched gene could be used only once, and each sampled gene inherited the weight of its corresponding clinical seed. The target node remained fixed, so its own network degree was shared by the observed and null propagation statistics. Matching controls expression-dependent seed topology within the available selected gene universe; it is not a population-randomization test.

We reported the observed score, null mean and standard deviation, standardized score Z = (observed − null mean)/null SD, and the one-sided empirical enrichment P value [1 + number(null score ≥ observed score)]/2,000. Benjamini–Hochberg correction was applied across all 50 available target–compartment hypotheses jointly within each graph specification and clinical-weight profile. If the null variance was zero, Z was recorded as undefined, rather than converting numerical ties into ranks. The minimum empirical P value was 0.0005. Null Z values summarize an irregular conditional reference distribution and were not interpreted as standard-normal clinical test statistics.

To assess graph-degree sensitivity, we repeated primary-graph null generation with log(1 + weighted degree) as a third standardized matching feature, using matched pools of 25, 50 and 100 genes. The full sensitivity results were retained, including changes in significance. This analysis specifically checks whether unusually connected clinical seeds can explain apparent target proximity; it does not remove the need for an expression-support requirement.

### Donor robustness, parameter stability and method comparison

We omitted each of the five donors in turn, reconstructed the mean correlation graph from the remaining four donors, and recomputed the null matching features and matched sets using only those four donors. Clinical weights and the inherited gene universe remained fixed. This is genuine donor deletion for the network estimation and null matching, not a bootstrap of cells mislabeled as independent donor validation. Within each compartment, candidates were ranked by decreasing null-adjusted Z among finite values. The primary stability descriptor required primary BH q < 0.05, a top-quartile candidate rank in at least 75% of the 18 graph specifications, and a top-quartile rank in at least four of the five donor-deletion analyses. This descriptor represents rank stability; it does not require independent biological replication or confirm a pharmacological mechanism.

For a target to retain expression-supported network priority, we additionally applied the previously defined equal-donor mean detection threshold of 10% in the same compartment. This inherited threshold was not relaxed to rescue low-expression diffusion hits. We preserved lower-expression candidates and their complete scores for interpretation and follow-up. Algorithm comparison included the unshrunk covariance ablation, alternate clinical-weight profiles, full-versus-donor-deleted rank correlations, and comparison with numerically evaluable scTenifoldKnk outputs. A scTenifoldKnk result at floating-point noise scale cannot be treated as a valid displaced-gene ranking; the main benchmark must therefore apply its separate run-level numerical-evaluability audit before any old-method correlation is interpreted. Both methods use the same reference atlas, so their agreement is an algorithmic-consistency assessment, not independent biological validation. No diagnostic AUC, causal effect or experimental knockout effect was inferred from network proximity.

## Results

### Clinical-protein-seeded diffusion identifies narrow, expression-limited bridges

The new diffusion analysis connected the clinical 13-protein signature to 50 evaluable target–compartment combinations while retaining the full 17-candidate universe. All primary scores passed the numerical checks, and the maximum linear-system residual across the 18 graph specifications was 8.88 × 10^(-16). The complete grid produced 3,600 candidate estimates across the four clinical-weight profiles. Five donor deletions yielded 250 candidate estimates; one omitted-donor configuration produced an undefined null-adjusted value and was retained as unavailable rather than ranked. These properties distinguish a solved network-proximity calculation from a numerically null virtual perturbation.

In the primary analysis, endothelial LRRK2 (Z = 27.424), endothelial TYK2 (Z = 5.670), and macrophage PDGFRA (Z = 7.299) each reached an empirical P value of 0.0005 and a joint-family BH q value of 0.00833. Endothelial LRRK2 and macrophage PDGFRA ranked in the top quartile in all 18 specifications and in four and five of the donor-deletion analyses, respectively. However, their equal-donor mean detection fractions were only 4.184% and 0.139%; endothelial TYK2 was detected in 5.324%. None met the inherited 10% same-compartment expression threshold. Thus, no pharmacological candidate qualified as an expression-supported, significant and rank-stable clinical-signature bridge in this atlas.

Additive seed decomposition further narrowed the biological interpretation. OTC and ALDOB accounted for 47.43% and 45.74% of the endothelial LRRK2 score. PCK2 accounted for 92.68% of the endothelial TYK2 score, and HPD accounted for 91.41% of the macrophage PDGFRA score. These concentrated contributions show why an apparently strong omnibus score need not represent broad disruption of a 13-protein program. The macrophage PDGFRA signal arose in a candidate detected in only two cells across two donors and warrants targeted expression verification before mechanistic follow-up.

The endothelial LRRK2 and TYK2 enrichment results remained below q = 0.05 when matching also included graph degree across pool sizes of 25, 50 and 100. Macrophage PDGFRA failed this criterion with the 100-gene matching pool. Neither stability analysis resolved the low-expression limitation. Moreover, LRRK2 donor-deletion Z values ranged from −0.985 to 7.569 and TYK2 values from −0.934 to 42.265, showing that high full-atlas scores do not imply uniformly positive donor-deletion support. FGFR1 had negative null-adjusted primary scores in hepatocytes, endothelial cells and macrophages (−0.631, −1.238 and −1.007, respectively); the new algorithm therefore did not reproduce a specific FGFR1-to-injury-signature enrichment.

Changing the clinical contrast or cohort while keeping the graph fixed yielded candidate-rank correlations of 0.838–0.982 with the primary weighting. Unshrunk Pearson and shrinkage-conditioned diffusion ranks correlated at 0.517 in hepatocytes, 0.821 in endothelial cells and 0.924 in macrophages, with maximum candidate-rank changes of 13, 7 and 5 places. The covariance choice therefore affected ranking but did not establish predictive superiority. Full-versus-donor-deleted rank correlations ranged from 0.471 to 0.747 in hepatocytes, 0.182 to 0.947 in endothelial cells, and 0.449 to 0.951 in macrophages. These results motivate reporting an uncertainty-aware experimental-priority set through the combined pharmacology, exposure-response and clinical-protein analysis, rather than naming a molecular toxicity mechanism on the basis of a single network rank.

The methodological contribution is the reproducible connection of continuous clinical protein effects, cell-compartment expression, donor-level perturbation of the reference graph, topology-matched nulls and target-level parameter stability. Its value is measurable through numerical evaluability, transport of candidate ranks across real donor deletions and transparent rejection of low-support bridges. It does not depend on claiming that covariance shrinkage or random-walk propagation is itself newly invented.

After excluding scTenifoldKnk runs at numerical-zero scale, only one hepatocyte candidate, eight endothelial candidates and six macrophage candidates had at least one evaluable legacy run. Conditional rank comparisons with diffusion were 0.381 for the eight endothelial candidates and 0.086 for the six macrophage candidates. Requiring three evaluable legacy runs left no hepatocyte candidate and only two candidates in each other compartment, precluding an informative three-seed method-correlation estimate. This audit prevents numerical noise from being mistaken for corroborating perturbation evidence and explains why improved evaluability is the defensible algorithmic gain in this dataset.

## Method references

1. Ledoit O, Wolf M. A well-conditioned estimator for large-dimensional covariance matrices. J Multivar Anal. 2004;88:365–411. doi:10.1016/S0047-259X(03)00096-4.
2. Köhler S, Bauer S, Horn D, Robinson PN. Walking the interactome for prioritization of candidate disease genes. Am J Hum Genet. 2008;82:949–958. doi:10.1016/j.ajhg.2008.02.013.
3. Vanunu O, Magger O, Ruppin E, Shlomi T, Sharan R. Associating genes and protein complexes with disease via network propagation. PLoS Comput Biol. 2010;6:e1000641. doi:10.1371/journal.pcbi.1000641.
4. MacParland SA, et al. Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations. Nat Commun. 2018;9:4383. doi:10.1038/s41467-018-06318-7.
''')
(p/'README.md').write_text('''# Network diffusion reanalysis

Run `python run_network_diffusion.py --package /path/to/original/recovered/package --output /path/to/output`, then `python build_corrected_legacy_benchmark.py --package /path/to/original/recovered/package --audit /path/to/vko_numerical_run_audit.csv --output /path/to/output`, followed by `python build_network_report.py` with this report script placed alongside the outputs. The original package is read only. Dependencies: Python 3.12, numpy, pandas, scipy, scikit-learn, threadpoolctl. The actual reference input files are checksummed in `network_diffusion_summary.json`.

`candidate_network_summary.csv` is the 17-target × 3-compartment master table. `network_specification_scores.csv` retains all 18 graph specifications and four clinical seed-weight profiles. `network_donor_lodo_scores.csv` contains actual donor omissions, with null matching recomputed without the omitted donor. `network_null_matching_sensitivity.csv` tests additional graph-degree matching at three pool sizes. `network_shrinkage_ablation.csv` is the unshrunk-Pearson ablation. `network_target_anchor_contributions.csv` decomposes each main score by clinical protein. The complete candidate lists include low-expression and unsupported results.

`network_evaluable_legacy_comparison.csv` and `network_evaluable_legacy_metrics.csv` exclude numerical-zero legacy runs. The older raw-percentile comparison is retained for audit only and marked invalid for primary comparison and must be read with the new scTenifoldKnk numerical-evaluability audit, because a floating-point noise result cannot constitute a meaningful old-method rank. No new raw data, cell/animal interventions, clinical outcomes or causal effects were generated.
''')
files=[]
for f in sorted(p.iterdir()):
 if f.is_file() and f.name!='network_output_manifest.csv':files.append({'file':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
pd.DataFrame(files).to_csv(p/'network_output_manifest.csv',index=False)
print(json.dumps(report,indent=2))
