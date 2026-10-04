All 24 returned Boltz complexes were inspected using the exact CIF and metrics paths in the two result indices. **None is ready for production MD as exported.** No inputs were changed and no PAE files were loaded.

| Target | Declared residues | Residues with coordinates | Missing native residue ranges |
|---|---:|---:|---|
| PDGFRA | 362 | 199 | 695–712, 728–779, 793–799, 855–878, 893–954 |
| FLT4 | 329 | 195 | 951–973, 985–995, 1074–1173 |

Each declared sequence matches its submitted construct, and every observed residue identity matches that sequence. The coordinate exports nevertheless contain four internal gaps for PDGFRA and two for FLT4. Their endpoints are separated by many angstroms. Normal peptide distances within the retained segments do not resolve these missing regions. The omission pattern is shared across all 12 ligands for each target.

**Confidence:** every atom B field is 100.000, and no local QA records are present. These fields cannot be used as pLDDT or colored as local confidence. Global confidence metrics are preserved in the JSON; local confidence in the pocket, motifs and insert cannot be determined from this audit.

| Priority complex | Ligand heavy atoms / bonds | Closest protein–ligand distance | Canonical hinge proximity |
|---|---:|---|---|
| PDGFRA–imatinib | 37 / 41 | Thr674 OG1–ligand N59: 2.807 Å | Cys677 backbone N–ligand N34: 2.930 Å |
| FLT4–tivozanib | 32 / 35 | Glu896 OE2–ligand N47: 2.823 Å | Cys930 backbone N–ligand N23: 2.975 Å |

These are distance-only contacts, not validated hydrogen bonds. PDGFRA numbering is CIF local residue +592; FLT4 is local residue +844. The PDGFRA canonical hinge is Glu675/Tyr676/Cys677, while Cys814/Val815 are separate catalytic-region residues. FLT4 hinge residues are Glu928/Phe929/Cys930. Both priority ligands make contacts in the kinase-pocket region; this does not validate their binding modes experimentally.

Across all 24 structures, ligand element counts, atom names, topology endpoints and input SMILES descriptors agree. No sub-2 Å protein–ligand heavy-atom overlaps or ligand bonds above 2 Å were found. This severe-overlap screen does not establish clash-free packing; atom-specific van der Waals contacts and force-field energies were not assessed. **Coordinate stereochemistry was not independently verified.** Matching metadata does not establish correct three-dimensional chirality, protonation or tautomer choice. Nerandomilast’s canonical descriptor can use charge-separated sulfoxide notation; do not treat that string difference alone as an identity failure.

Nintedanib, sunitinib, ponatinib and MAZ51 were filtered before returned structures were produced. Absence is not evidence of nonbinding. Nintedanib appears as an input reference ligand, but no nintedanib complex is available from these jobs for MD.

MD preparation would require reconstruction or an explicitly justified alternative construct, deliberate terminal chemistry, hydrogens/protonation, force fields and ligand parameterization. A complete sequence declaration, favorable global confidence or successful simulation setup cannot substitute for coordinate completeness or validate inhibition and liver-fibrosis efficacy. 3C7Q remains a KDR/VEGFR2 homolog reference.

Detailed per-model results and source hashes: Boltz_All24_Structure_QA.json. Native-numbered contacts: Boltz_All24_Contacts.csv. The priority-pair JSON also retains exact ligand atom IDs, coordinates and explicit CIF bond types for figure rendering.
