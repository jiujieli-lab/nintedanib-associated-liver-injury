# Workflow inspection and draft content specification

Status: inspection only. No figures have been changed. Revised story confirmation is pending.

## Verified edit scope

Only final Figure 1A, 5A, 6A, and 7A are workflow panels. **Figure 2A is a data plot**, not a workflow. The source pixel arrays were loaded successfully from the validated original TIFFs and extracted as temporary PNG previews.

| Final panel | Source panel | Original crop (L,T,R,B), px | Existing reassembled position (x,y,w,h), px |
|---|---|---|---|
| Figure_1A | Figure_1A | [0, 0, 2008, 211] | [12, 12, 2000, 198] |
| Figure_5A | Figure_4A | [0, 0, 2008, 228] | [12, 12, 2006, 216] |
| Figure_6A | Figure_5A | [0, 0, 2008, 219] | [12, 12, 2006, 216] |
| Figure_7A | Figure_6A | [0, 0, 2008, 221] | [12, 12, 2006, 214] |

## Recommended geometry

Use a 380-pixel-high icon strip instead of the present 198–216-pixel strip. On each reassembled figure, retain A's top-left coordinate and width. Increase the canvas height by the difference, and translate every lower panel downward by the same difference. Keep every lower panel's pixel data, dimensions, and x-coordinate unchanged. This preserves measurement figures exactly while allowing biomedical icons to be legible. A strict unchanged canvas is feasible but gives substantially smaller icons than the examples.

## Aesthetic direction

The references communicate through recognizable biomedical objects, ordered stages, and sparse captions. Use crisp organ/cell/protein silhouettes with restrained blue-gray structure and warm liver/serum accents; simple thin arrows; white background; short headings; consistent outline weight. Avoid decorative generated measurement plots, ornamental molecules, or dense repeated blue paragraph boxes. Use a wider biological illustration and shorter text in each stage.

## Scientific content and cautions

### Figure_1A

Clinical injury evidence and liver-cell context prioritize receptors for experimental investigation.

1. **Clinical injury** — FAERS signal / 13 serum proteins. Icons: human silhouette, serum tube, liver.
2. **Drug + liver cells** — Pharmacology + HepG2 / 5-donor liver atlas. Icons: nintedanib capsule, receptor on cell membrane, liver-cell cluster.
3. **Target priorities** — Calibrated integration / PDGFRA • FLT4. Icons: weighted gene network, two receptor nodes.
4. **Independent evaluation** — Liver exposure challenge / Pulmonary response context. Icons: mouse beside liver, lung.

No rescue, protection, causal DILI mechanism or successful validation claim. Final title/stages await updated analysis story.

### Figure_5A

Cell localization and raw-distance calibration constrain the interpretation of virtual knockouts.

1. **Liver-cell context** — 5 human donors / Compartment + detection. Icons: human liver, three cell phenotypes.
2. **Virtual knockout** — 210 archived runs / Raw network distances. Icons: small gene network with one gray removed node.
3. **Numerical audit** — Zero threshold ≤10⁻¹² / 132 numerical-zero runs. Icons: screen/filter with gray nodes rejected.
4. **Clinical projection** — 29 nonzero target runs / Matched-null comparison. Icons: protein anchors around retained network.

132/210 numerical-zero includes non-target runs; 29 valid pharmacologic-target runs are not 210−132. Network edges are computational associations; no successful biological enrichment claim.

### Figure_6A

Donor-balanced clinical-protein-guided integration yields stable receptor priorities.

1. **Liver + clinical anchors** — Donor-balanced networks / 13 protein weights. Icons: donor cells, serum tube, protein chain.
2. **Calibrated proximity** — 1,999 matched null sets / Expression + detection. Icons: weighted network, gray background nodes.
3. **Evidence integration** — Pharmacology · network / HepG2 · STRING. Icons: four small aligned evidence symbols, fusion ring.
4. **Frozen priorities** — 17 targets · 6 algorithms / PDGFRA • FLT4. Icons: two highlighted receptor nodes, lock.

The six algorithms include historical comparator and five revised rules. Stable rank is experimental priority, not causal proof. Do not show macrophage PDGFRA as robust expression.

### Figure_7A

An independent whole-animal liver dataset challenges frozen target ranks.

1. **Frozen candidates** — 17 targets / 6 algorithms. Icons: lock, small receptor panel.
2. **Independent mouse liver** — Healthy 3 · disease 9 / Nintedanib 10. Icons: three mouse group silhouettes, liver.
3. **Animal-level analysis** — Genome-wide + target family / 92,378 allocations. Icons: transcript reads on screen, mouse-linked samples.
4. **External rank challenge** — 12 measured targets / Priority vs. hepatic effect. Icons: ranked nodes beside liver response symbol.

22 animals in all groups; 19 in disease-vs-treatment contrast and leave-one-animal-out sensitivity. GSE125975 measures pharmacodynamics, not DILI rescue. Kdr significant only candidate-family, not genome-wide.

Source of truth: deliverables/Manuscript/Revised_Manuscript_Clean.docx; work/audit/panel_crop_manifest.csv; work/panels/reassembled/*_layout.json.
