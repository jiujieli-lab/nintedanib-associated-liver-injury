PDGFRA / FLT4 model-conditional network perturbation extension

This folder contains executed quantitative analyses, complete positive and null outputs, figure sources and manuscript-ready text. It contains no State checkpoint predictions, measured genetic interventions, clinical causal effects, molecular dynamics or biological rescue experiments.

Scientific status
All nine primary matched-control magnitude tests were nonsignificant after BH correction (minimum q=0.5294). All three paired nonadditivity tests were nonsignificant after BH correction (minimum q=0.1011). The original candidate-prioritization conclusion remains unchanged. No significance threshold, target set or graph setting was changed to obtain a positive result. Numerical restoration of the original graph is an implementation check only.

Inputs
The supplied GSE115469 normalized compartment matrices contain 3,501 hepatocytes, 844 endothelial cells and 1,192 macrophages from five donors. They are inherited 705/717/716-gene analysis matrices, not a reprocessed full atlas. PDGFRA is sparsely detected in these compartments; the 37 original stellate cells were not available as a suitable network matrix and were not expanded or simulated. No new external cell measurements were generated. Input_manifest.csv records SHA-256 checksums.

Reproduction
Place this directory alongside an extracted data directory preserving the supplied package subfolders:
 data/Additional_file_5_Liver_and_Lung_Single_Cell_Data/Source_Data/Liver_Single_Cell
 data/Additional_file_4_Clinical_Evidence_and_Analysis_Code/Analysis/network_diffusion
Then run:
 python run_target_pair_perturbation.py
 python make_figure_s11.py
 python validate_extension.py
Python dependencies: numpy, pandas, scipy, scikit-learn, threadpoolctl, matplotlib, Pillow. The calculation uses CPU double precision and limits BLAS threads. The tested runtime is recorded in validation_summary.json.

Interpretation
The intervention deletes target-incident edges in an undirected donor-balanced transcript-covariance graph and recalculates clinical-seeded random-walk propagation. Its deltas are propagation-score changes, not expression fold changes. Reciprocal target-score shifts cannot identify receptor-to-receptor regulation. The exact dual-minus-single interaction contrast quantifies nonadditivity of the graph operator only. No-change and restored-graph controls establish numerical consistency, not predictive or biological validity. Donor omission is within-atlas sensitivity using an inherited gene universe, not held-out intervention validation or independent replication.

The attenuation series scales each edge incident to either target once, including an edge connecting both target nodes. This explicit definition avoids double attenuation of a shared edge. Primary 100% deletion values are unaffected by that partial-attenuation definition.

The accompanying Figure S11 contains nine panels; independent PNG panel exports are included in Individual_Panels. The complete figure is available as 1200-dpi TIFF, vector PDF and preview PNG. All nonsignificant findings and full null distributions are retained in the source tables. manuscript_insertion.txt contains the exact executed methods, results and full figure legend for integration into the parent manuscript.
