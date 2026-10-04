Completed targeted analyses and exact input data

Unzip while preserving all relative folders. From this directory run:
python new_virtual/run_target_pair_perturbation.py
python new_virtual/make_figure_s11.py
python new_virtual/validate_extension.py
python binding_analysis/analyze_binding_evidence.py

The Python environment needs numpy, scipy, pandas, scikit-learn, threadpoolctl, matplotlib and Pillow. Actual package versions appear in ENVIRONMENT.json. Calculations use the supplied frozen inputs. Prediction responses from Inductive Bio are archived real service results; rerunning the local binding script does not silently call the service. No conventional docking, Boltz prediction or MD trajectory is included in these completed results.

Supplementary Tables N18 and N19 are new_virtual/gene_level_perturbations.csv.gz and new_virtual/matched_null_metrics.csv.gz.
