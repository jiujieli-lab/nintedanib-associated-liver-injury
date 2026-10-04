BOLTZMOL MAIN-FIGURE MOLECULAR PANELS

Panel_E_BoltzMol_PDGFRA_imatinib and Panel_F_BoltzMol_FLT4_tivozanib:
PNG, PDF, and SVG, each 2.25 × 2.4 inches; PNG 1350 × 1440 at 600 dpi.
No panel letters are embedded. Existing docking assets were retained unchanged.

SOURCE AND IDENTITY
Both views use the delivered BoltzMol predicted CIF coordinates directly, with
full decimal precision. Ligand topology is rebuilt from chem_comp_atom/bond,
sanitized in RDKit, and checked against the result metadata and library SMILES.
Protein residue identities match both the submitted sequence and the archived
canonical, UniProt-indexed full-length AlphaFold PDB sequence. Those AlphaFold
coordinates are used only for independent residue-identity verification, not
for the displayed complex geometry.

The PDGFRA–imatinib export contains 199 coordinate-bearing residues of the 362
submitted residues; the FLT4–tivozanib export contains 195 of 329. Native numbers
are local label_seq_id +592 and +844, respectively. No missing coordinates were
modeled or added. These are cropped predicted pocket views, not full structures.

RENDERING
A single shared rigid orthographic camera projects each complex. There is no
ligand alignment or coordinate editing. Uniform gray local C-alpha traces break
at missing residue numbers, abnormal C-alpha separation >4.5 Å, and view cutoffs.
Cubic smoothing is applied only within contiguous visible trace fragments.
All source B fields equal 100.000 and are NOT used as confidence or pLDDT.
Ligand carbon is teal, nitrogen navy, oxygen coral. Chlorine is explicitly labeled.
Actual bond orders are retained, with a Kekulé drawing of aromatic connectivity.
Selected protein residue bonds are gray. Molecular stereochemistry was not
independently validated from the 3-D coordinates and is not claimed as a QA result.

SELECTED CONTACTS (minimum heavy-atom distances; not hydrogen-bond assignments)
PDGFRA imatinib:
  Glu644 OE2 – ligand N58, 2.8324325266 Å
  Cys677 N   – ligand N34, 2.9300268100 Å
  Asp836 N   – ligand O32, 3.0323150302 Å
FLT4 tivozanib:
  Glu896 OE2 – ligand N47, 2.8225342189 Å
  Cys930 N   – ligand N23, 2.9745762274 Å
  Asp1055 N  – ligand O20, 2.9270787949 Å

Dashed segments join exact contact atom coordinates. Thin solid leaders identify
residues. All 37 priority-pair contact residues below 4 Å match the independently
provided Boltz_Priority_Pair_Contacts.csv by atom pair and numerical distance.
Only three selected contacts per view are labeled for readability.

BoltzMol_molecular_provenance.json records full source hashes, cameras, topology
verification, residue mapping, complete contact lists, and output dimensions.
BoltzMol_molecular_preview.png was visually inspected after final export.

Reproduce with:
work/jcm_revision/docking/env/bin/python \
  work/jcm_revision/figures/screening/molecular/render_boltzmol_panels.py

The script imports shared drawing helpers from render_molecular_panels.py.
Parent compositor owns final figure numbering, caption, and publication.
