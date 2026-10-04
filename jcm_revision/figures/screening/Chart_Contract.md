# Figure 7 and Supplementary Figure S13 chart contract

Status: design contract; docking is still running and must not be represented as complete.

## Figure 7

Question: How do experimentally measured, same-study binding profiles and explicitly predictive computational measurements differentiate a fixed mechanistic comparator library?

Takeaway: Experimental endpoint-specific evidence and structure-dependent docking/property predictions provide complementary comparator profiles; none establishes hepatoprotection or receptor-specific clinical safety.

Surface: standalone static publication composite; Matplotlib and genuine molecular coordinates; approximately 7.2 × 8.8 in, 600 dpi RGB TIFF/PNG, PDF. No image-generated measurement charts or invented molecular poses.

| Panel | Comparison and chart | Grain and sufficiency | Encoding and safeguards |
|---|---|---|---|
| A | Same-study experimental Kd; paired horizontal log-scale point profiles | Davis 2011, six compounds × two targets; 12 records | PDGFRA teal circles and FLT4 coral squares; open symbols identify eight potential-duplicate records; >10,000 nM remains right-censored. No cross-study pooling and no conversion of censored values to exact estimates. |
| B | Predicted logD at pH 7.4; horizontal point and model-bound profiles | All 16 exact parent structures, same Inductive Bio model v1.7.0 | Index drug highlighted; fixed library order; bounds explicitly not CIs; no safety interpretation. |
| C | PDGFRA target-specific docking profile; individual-seed points and median/range | Exactly 16 compounds × 3 completed searches | Fixed candidate order, seed shapes retained, no significance testing; title identifies experimental 6JOL rigid receptor. |
| D | FLT4 target-specific docking profile; individual-seed points and median/range | Exactly 16 compounds × 3 completed searches | Separate axis; predicted AlphaFold receptor explicitly identified; no cross-target affinity or selectivity inference. |
| E | Crystallographic imatinib and redocked pose overlay | Actual 6JOL reference and top executed redock pose | Common receptor coordinate frame; actual symmetry-aware heavy-atom RMSD; no ligand alignment. |
| F | PDGFRA nintedanib predicted pose ambiguity | Predefined seed 104729, top-ranked and hinge-near retained alternative | Explicit pose ranks, actual score gap and hinge distances; do not depict only the canonical-looking alternative. |
| G | FLT4 nintedanib predicted pose | Predefined seed 104729 top-ranked pose in actual predicted receptor | Contact distances are heavy-atom proximities, not automatically hydrogen bonds; no experimental self-redock control available. |

## Supplementary Figure S13

Retain the full library and computational-search detail, grouped into reproducible panels: endpoint-separated biochemical evidence-count matrix; two target-specific full raw seed-score matrices; two top-pose RMSD consistency profiles; predicted acidic and basic pKa across all 16 with provider bounds, out-of-domain and low-confidence markings; all nintedanib retained-pose ranks versus score and hinge proximity. Full source tables accompany the figure for exact lookup. No missing experimental endpoint is encoded as inactivity.

## Palette and graphical rules

White background; near-black text; muted teal #177E89 for PDGFRA, coral #D66A52 for FLT4; navy #244F70 index highlight, neutral gray comparators. Target shapes differ in addition to color. Seeds differ by shape, not biological-group color. Native font size approximately 7.5–8.5 pt at final print size. No bar charts imply repeated computational seeds are biological replicates.

## Gates and QA

Before finalizing C/D and S13, require all 96 planned result.json outputs (16 compounds × 2 receptors × 3 seeds), no duplicate name/target/seed keys, and exactly three actual scores per compound–target pair. Capture source hashes. Inspect exported TIFF after the writer process exits, require nonzero IFD, verify PNG/TIFF pixel equality, and render/inspect PDF. Preserve the approved older workflow composites and their numbering.

## Revised main topology after paid Boltz authorization

Final assembly is on hold until both 16-compound Boltz jobs return actual outputs and field semantics are checked. The preferred main figure is six panels: A matched Davis Kd; B and C target-specific Boltz profiles with all 16 candidates; D imatinib redocking overlay; E PDGFRA nintedanib pose ambiguity; F FLT4 nintedanib predicted pose. Boltz scores/confidence will retain the exact returned definition and cannot be represented as experimental binding affinity. Quantitative Boltz panels require both complete job outputs and verified compound identity mapping.

Move full docking seeds, pose consistency, endpoint coverage and pose sensitivity to S13. Move full 16-compound logD, acidic/basic pKa and applicability/low-confidence flags to S14 if two supplement pages are needed for legibility. This supersedes the original dense main layout, while retaining all evidence and numerical gates.

Boltz update: returned engine is pipelineboltzmol v1.0, so panel labels will use BoltzMol, not Boltz-2. Main B/C use the provider field `binding_confidence` as a dimensionless score, not a calibrated binding probability, physical affinity, or percent. Final-disposition gate requires all 16 input identities, allowing explicit provider-filtered rows with no numerical marker. Four default-filtered rows per target have been reported while the jobs remain in progress; final counts must be verified from completed outputs. Filtered rows must never be plotted as zero.
