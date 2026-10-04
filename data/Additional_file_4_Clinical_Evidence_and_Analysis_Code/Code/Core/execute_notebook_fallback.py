#!/usr/bin/env python3
"""Execute a notebook top-to-bottom in-process when ZeroMQ is sandbox-blocked.

This runner preserves code-cell order, captures stdout and explicitly displayed
pandas objects, stops on the first exception, and writes execution metadata.
It is intentionally narrow and sufficient for the audit notebook, which does
not use magics, shell escapes, widgets, or asynchronous cells.
"""

from __future__ import annotations
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import traceback


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
NB_PATH = PACKAGE_ROOT / "Public_Data_Reproducibility_Notebook.ipynb"
REPORT_PATH = NB_PATH.with_name("Notebook_Validation_Report.md")
os.chdir(PACKAGE_ROOT)
nb = json.loads(NB_PATH.read_text(encoding="utf-8"))
namespace = {"__name__": "__main__"}
execution_count = 0


def write_notebook() -> None:
    NB_PATH.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")


def new_output(output_type: str, **kwargs):
    return {"output_type": output_type, **kwargs}


def display(value):
    outputs = namespace.setdefault("__cell_display_outputs__", [])
    data = {"text/plain": repr(value)}
    if hasattr(value, "to_html"):
        try:
            data["text/html"] = value.to_html(max_rows=40, max_cols=30)
        except TypeError:
            data["text/html"] = value.to_html()
    outputs.append(new_output(output_type="display_data", data=data, metadata={}))


namespace["display"] = display
for cell in nb["cells"]:
    if cell["cell_type"] != "code":
        continue
    execution_count += 1
    cell["execution_count"] = execution_count
    cell["outputs"] = []
    namespace["__cell_display_outputs__"] = []
    stdout = StringIO(); stderr = StringIO()
    try:
        with redirect_stdout(stdout), redirect_stderr(stderr):
            exec(compile(cell["source"], f"{NB_PATH.name}:cell{execution_count}", "exec"), namespace)
    except Exception as exc:
        tb = traceback.format_exc().splitlines()
        if stdout.getvalue():
            cell["outputs"].append(new_output(output_type="stream", name="stdout", text=stdout.getvalue()))
        if stderr.getvalue():
            cell["outputs"].append(new_output(output_type="stream", name="stderr", text=stderr.getvalue()))
        cell["outputs"].extend(namespace["__cell_display_outputs__"])
        cell["outputs"].append(new_output(
            output_type="error", ename=type(exc).__name__, evalue=str(exc), traceback=tb))
        nb["metadata"]["execution"] = {
            "status": "failed", "runner": "in-process ordered fallback",
            "failed_cell": execution_count,
            "utc": datetime.now(timezone.utc).isoformat(),
            "reason_standard_kernel_unavailable": "Jupyter and nbformat were not installed in the release-build runtime",
        }
        write_notebook()
        raise
    if stdout.getvalue():
        cell["outputs"].append(new_output(output_type="stream", name="stdout", text=stdout.getvalue()))
    if stderr.getvalue():
        cell["outputs"].append(new_output(output_type="stream", name="stderr", text=stderr.getvalue()))
    cell["outputs"].extend(namespace["__cell_display_outputs__"])

nb["metadata"]["execution"] = {
    "status": "passed",
    "runner": "in-process ordered fallback",
    "code_cells_executed": execution_count,
    "utc": datetime.now(timezone.utc).isoformat(),
    "reason_standard_kernel_unavailable": "Jupyter and nbformat were not installed in the release-build runtime",
}
assert nb.get("nbformat") == 4 and isinstance(nb.get("cells"), list)
cell_ids = [cell.get("id") for cell in nb["cells"]]
assert all(isinstance(cell_id, str) and cell_id for cell_id in cell_ids)
assert len(cell_ids) == len(set(cell_ids))
write_notebook()

error_outputs = [
    output
    for cell in nb["cells"]
    for output in cell.get("outputs", [])
    if output.get("output_type") == "error"
]
stderr_outputs = [
    output
    for cell in nb["cells"]
    for output in cell.get("outputs", [])
    if output.get("output_type") == "stream" and output.get("name") == "stderr"
]
output_count = sum(len(cell.get("outputs", [])) for cell in nb["cells"])
code_cells = [cell for cell in nb["cells"] if cell["cell_type"] == "code"]
expected_counts = list(range(1, len(code_cells) + 1))
observed_counts = [cell.get("execution_count") for cell in code_cells]
if error_outputs or stderr_outputs or observed_counts != expected_counts:
    raise RuntimeError("Post-execution notebook structure validation failed")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


notebook_sha = sha256(NB_PATH)
builder_path = PACKAGE_ROOT / "Code" / "Core" / "build_reproducibility_notebook.py"
runner_path = Path(__file__).resolve()
integration_manifest_sha = namespace.get("integration_manifest_observed", "not recorded")
core_manifest_sha = namespace.get("core_manifest_observed", "not recorded")
s6_manifest_sha = namespace.get("s6_manifest_observed", "not recorded")
s7_manifest_sha = namespace.get("s7_manifest_observed", "not recorded")
REPORT_PATH.write_text(
    "\n".join(
        [
            "# Notebook Validation Report",
            "",
            "## Overall assessment: PASS",
            "",
            f"- Notebook: `{NB_PATH.relative_to(PACKAGE_ROOT)}`",
            f"- Final execution (UTC): {nb['metadata']['execution']['utc']}",
            "- Runner: deterministic in-process ordered fallback; Jupyter and nbformat were not installed in the release-build runtime.",
            f"- Total cells: {len(nb['cells'])}",
            f"- Code cells executed: {len(code_cells)}",
            f"- Captured outputs: {output_count}",
            f"- Error outputs: {len(error_outputs)}",
            f"- Stderr outputs: {len(stderr_outputs)}",
            f"- Sequential execution counts: {observed_counts == expected_counts}",
            f"- Unique nbformat 4.5 cell IDs: {len(cell_ids) == len(set(cell_ids)) == len(nb['cells'])}",
            f"- Notebook SHA256: `{notebook_sha}`",
            f"- Builder SHA256: `{sha256(builder_path)}`",
            f"- Fallback runner SHA256: `{sha256(runner_path)}`",
            f"- Integration manifest SHA256: `{integration_manifest_sha}`",
            f"- Figure 6 core manifest SHA256: `{core_manifest_sha}`",
            f"- GSE151374/S6 manifest SHA256: `{s6_manifest_sha}`",
            f"- PXD052594/S7 manifest SHA256: `{s7_manifest_sha}`",
            "",
            "## Validation scope",
            "",
            "The notebook was rebuilt from its source script and executed top-to-bottom directly from the extracted submission directory. Its assertions verify the frozen FAERS endpoint table, the sole FDR-significant human DILI proteomic separator, six LINCS doses reconstructed from packaged per-gene values, liver-atlas dimensions, three virtual-knockout seeds, the same-compartment T2 gate and sensitivity analysis, the 44/44 integration validation report, all four packaged Figure 5 formats, and exact participant-level sign-flip P values recomputed from all 2^7 assignments. Public preclinical sections read frozen compact CSV/JSON products and detached manifest audits without rerunning heavy upstream analysis or opening raw H5 files/workbooks.",
            "",
            "GSE151374/S6 validates 18 pooled libraries (three mice pooled per library; n=3 libraries in each of six condition–time cells), library-level pseudobulk/exact-test constraints, and FGFR1 versus FBP1/GSTA1/OTC role separation. PXD052594/S7 validates the unconditioned individual-animal proteome (n=10 versus 13), the separate table-only phosphoproteome subset (n=5 versus 5), absence of radiomic response-cluster conditioning, and target/anchor detection boundaries.",
            "",
            "## Interpretation boundaries",
            "",
            "FAERS supports disproportionality, not incidence or causality. Human serum proteomics is multi-drug rather than nintedanib-specific. Virtual-knockout displacement is unsigned. GSE151374 supports pooled-library, not cell- or mouse-level, inference. PXD052594 is a pulmonary-response context and its phosphoproteome is a separate 5+5 subset, not liver-toxicity proof. Pooled-slice and replicate-provenance-ambiguous preclinical sources remain exploratory or descriptive. GSE299128 has no adjudicated DILI outcome.",
            "",
        ]
    ),
    encoding="utf-8",
)
print(f"PASS: executed {execution_count} code cells with {output_count} outputs: {NB_PATH}")
print(f"PASS: validation report written: {REPORT_PATH}")
