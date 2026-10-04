MOLECULAR PANEL ASSETS — EXACT COORDINATE PROJECTIONS

Exports: 2.25 × 2.4 inches, white background, PNG 1350 × 1440 pixels at 600 dpi;
editable SVG and PDF companions. No panel letters are embedded. Original E/F/G
filenames are asset identifiers and may be remapped by the final compositor.

Panel_E_redocking: 6JOL crystal imatinib (teal) over seed-104729 top redocked
imatinib (coral), in the unchanged PDGFRA receptor coordinate frame. There was
NO ligand alignment. Recorded RMSD is 0.5742745619518652 Å. Recomputing from the
rounded exported SDF gives 0.5742743952859829 Å; both render as 0.574 Å.

Panel_F_PDGFRA_pose_ambiguity: seed-104729 nintedanib pose 1 (teal) and pose 3
(coral) in the same receptor frame, with scores −8.502 and −8.341 kcal/mol.
The rendered 8.43/3.38 Å distances are the nearest heavy atoms to Thr674,
the gatekeeper. They are NOT minimum distances to the actual 675–677 hinge.
Both poses are predictions; pose 3 is not a proven canonical binding pose.

Panel_G_FLT4_predicted_pocket: seed-104729 top nintedanib pose in the AlphaFold
FLT4 receptor. Four selected nearest heavy-atom contacts are drawn as dashed
segments at their actual coordinate endpoints: Glu896 3.26 Å, Val927 3.21 Å,
Asn934 2.87 Å, and Asp1055 3.23 Å. Thin solid leaders identify residues.
These are proximity contacts, with NO hydrogen-bond assignment. Carbon is teal,
nitrogen navy, oxygen coral; selected protein residue bonds are neutral gray.

All molecule coordinates remain unchanged. A single shared orthographic camera
rotation is applied to every object in each panel. The gray background is an
interpolated local C-alpha trace, not a claimed secondary-structure assignment.
Hydrogens are hidden. Actual SDF bond orders are retained; aromatic bonds use a
Kekulé representation. Callout positions are editorial offsets, not coordinates.

Reproduction:
work/jcm_revision/docking/env/bin/python \
  work/jcm_revision/figures/screening/molecular/render_molecular_panels.py

Scientific QA: RMSD verified against result.json; Thr674 distances reproduced
from PDB/SDF and checked against updated hinge QC; all four FLT4 distances match
Nintedanib_Docking_Contacts.csv to rounding and are <4 Å. Source hashes, camera
matrices, precise atom identities, and distances are in molecular_panel_provenance.json.
PNG outputs and the three-panel preview were visually inspected after export.

These are prespecified seed-104729 illustrations; they do not establish final
three-seed outcomes. Parent compositor owns final numbering and publication.
