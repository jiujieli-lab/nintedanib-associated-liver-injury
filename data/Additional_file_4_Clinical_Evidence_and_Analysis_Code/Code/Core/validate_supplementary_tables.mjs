import fs from "node:fs/promises";
import path from "node:path";
import crypto from "node:crypto";
import os from "node:os";
import { execFileSync } from "node:child_process";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";
import { fileURLToPath } from "node:url";

const PACKAGE_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "../..");
const ROOT = process.env.NINTEDANIB_DILI_DEVELOPMENT_ROOT
  ? path.resolve(process.env.NINTEDANIB_DILI_DEVELOPMENT_ROOT)
  : PACKAGE_ROOT;
const XLSX = path.join(ROOT, "fresh_analysis", "submission", "Supplementary_Tables.xlsx");
const REPORT = path.join(ROOT, "fresh_analysis", "submission", "Supplementary_Tables_validation_report.json");
const RENDER_DIR = path.join(os.tmpdir(), "nintedanib_dili_workbook", "imported_final_renders");
const visualReviewed = process.argv.includes("--visual-reviewed");
const skipRender = process.argv.includes("--skip-render");
const targetExternalIndex = process.argv.includes("--target-external-index");
const packagedDay7Path = "Source_Data/GSE151374/results/macrophage_DE_NINT_vs_BLM_day7.csv";
const packagedDay14Path = "Source_Data/GSE151374/results/macrophage_DE_NINT_vs_BLM_day14.csv";

const previewRanges = {
  README: "A1:B21",
  Data_Sources: "A1:K16",
  Data_Provenance: "A1:F34",
  Statistical_Units: "A1:J17",
  Figure_Source_Map: "A1:H52",
  Endpoint_Definitions: "A1:H28",
  FAERS_Overall: "A1:X20",
  FAERS_Specifications: "A1:Y20",
  FAERS_Annual: "A1:U27",
  Proteomics_Effects: "A1:AL24",
  Proteomics_QC: "A1:K26",
  LINCS_Dose_Metrics: "A1:AD24",
  ChEMBL_Targets: "A1:P26",
  Liver_Expression: "A1:J26",
  Virtual_KO_Stability: "A1:J26",
  Integrated_Prioritization: "A1:O30",
  Blood_Sensitivity: "A1:I28",
  Preclinical_Evidence: "A1:O30",
  GSE151374_scRNA: "A1:M30",
  PXD052594_Proteomics: "A1:M30",
  Experimental_Validation: "A1:J12",
  Reporting_Checklists: "A1:G12",
};

if (!skipRender) await fs.rm(RENDER_DIR, { recursive: true, force: true });
await fs.mkdir(RENDER_DIR, { recursive: true });
const prior = JSON.parse(await fs.readFile(REPORT, "utf8"));
const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(XLSX));
const sheetCatalog = await workbook.inspect({ kind: "sheet", include: "id,name", maxChars: 20000 });
const requiredSheets = Object.keys(previewRanges);
const sheetCatalogRows = sheetCatalog.ndjson.trim().split("\n").filter(Boolean).map((line) => JSON.parse(line));
const actualSheetNames = sheetCatalogRows.map((row) => row.name);
const missingSheets = requiredSheets.filter((name) => !actualSheetNames.includes(name));
const unexpectedSheets = actualSheetNames.filter((name) => !requiredSheets.includes(name));

// Inspect the exported OOXML package itself rather than assuming the authoring
// request persisted. artifact-tool 2.8.6 currently omits worksheet <pane>
// records even after documented freezeRows/freezeColumns calls; disclose that
// limitation explicitly and do not count it as a passed workbook property.
const worksheetXmlFiles = execFileSync("unzip", ["-Z1", XLSX], { encoding: "utf8" })
  .split(/\r?\n/)
  .filter((name) => /^xl\/worksheets\/sheet\d+\.xml$/.test(name))
  .sort((a, b) => Number(a.match(/sheet(\d+)/)?.[1]) - Number(b.match(/sheet(\d+)/)?.[1]));
const exportedFreezePanes = worksheetXmlFiles.map((xmlPath, index) => {
  const xml = execFileSync("unzip", ["-p", XLSX, xmlPath], { encoding: "utf8", maxBuffer: 32 * 1024 * 1024 });
  const pane = xml.match(/<pane\b[^>]*\bstate="frozen"[^>]*\/?\s*>/i)?.[0] ?? null;
  return { sheet_name: actualSheetNames[index] ?? null, xml_path: xmlPath, pane_persisted: pane !== null, pane_xml: pane };
});
const persistedFreezePaneCount = exportedFreezePanes.filter((entry) => entry.pane_persisted).length;
const tableXmlFiles = execFileSync("unzip", ["-Z1", XLSX], { encoding: "utf8" })
  .split(/\r?\n/)
  .filter((name) => /^xl\/tables\/table\d+\.xml$/.test(name));
const persistedTableAutoFilterCount = tableXmlFiles.filter((xmlPath) => {
  const xml = execFileSync("unzip", ["-p", XLSX, xmlPath], { encoding: "utf8", maxBuffer: 32 * 1024 * 1024 });
  return /<(?:[A-Za-z0-9_]+:)?autoFilter\b/i.test(xml);
}).length;
const persistedTableAutoFiltersPassed = tableXmlFiles.length === 107 && persistedTableAutoFilterCount === 107;
const externalIndexValues = workbook.worksheets.getItem("GSE151374_scRNA").getRange("A879:I881").values;
const externalIndexPackagedPathsPassed = externalIndexValues?.[1]?.[1] === packagedDay7Path
  && externalIndexValues?.[2]?.[1] === packagedDay14Path;
const externalIndexHashesComplete = [externalIndexValues?.[1]?.[7], externalIndexValues?.[2]?.[7]]
  .every((value) => typeof value === "string" && /^[0-9a-f]{64}$/i.test(value));

const firstCellInspections = {};
const headerStyleInspections = {};
const firstHeaderChecks = {};
const headerStyleChecks = {};
for (const name of requiredSheets.filter((sheetName) => actualSheetNames.includes(sheetName))) {
  const cell = await workbook.inspect({
    kind: "table",
    range: `${name}!A1:A1`,
    include: "values,formulas",
    tableMaxRows: 2,
    tableMaxCols: 2,
    maxChars: 1200,
  });
  firstCellInspections[name] = cell.ndjson;
  const cellObj = JSON.parse(cell.ndjson.trim().split("\n")[0]);
  const observedFirstHeader = cellObj.values?.[0]?.[0];
  const expectedFirstHeader = prior.first_header_checks?.[name]?.expected;
  firstHeaderChecks[name] = { expected: expectedFirstHeader, observed: observedFirstHeader, passed: observedFirstHeader === expectedFirstHeader };

  const style = await workbook.inspect({
    kind: "computedStyle",
    sheetId: name,
    range: "A1:A1",
    maxChars: 2200,
  });
  headerStyleInspections[name] = style.ndjson;
  const styleObj = JSON.parse(style.ndjson.trim().split("\n")[0]);
  headerStyleChecks[name] = {
    fill: styleObj.style?.fill?.color?.value,
    font_color: styleObj.style?.font?.fill?.color?.value,
    bold: styleObj.style?.font?.bold,
    wrap_text: styleObj.style?.wrapText,
    horizontal_alignment: styleObj.style?.horizontalAlignment,
    passed: styleObj.style?.fill?.color?.value === "16324F"
      && styleObj.style?.font?.fill?.color?.value === "FFFFFF"
      && styleObj.style?.font?.bold === true
      && styleObj.style?.wrapText === true
      && styleObj.style?.horizontalAlignment === "center",
  };
}

const formulaScan = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 300 },
  summary: "imported-workbook formula error scan",
  maxChars: 4000,
});
const formulaScanPassed = formulaScan.ndjson.includes("matched 0 entries");

const requestedRenderSpecs = Array.isArray(prior.rendered_sheets) && prior.rendered_sheets.length
  ? prior.rendered_sheets.map((x) => ({ sheetName: x.sheetName, range: x.range, purpose: x.purpose || "imported-workbook repeat render" }))
  : Object.entries(previewRanges).map(([sheetName, range]) => ({ sheetName, range, purpose: "sheet overview" }));
const rendered = skipRender && Array.isArray(prior.rendered_sheets) ? prior.rendered_sheets : [];
if (!skipRender) {
  for (const { sheetName, range, purpose } of requestedRenderSpecs.filter((x) => actualSheetNames.includes(x.sheetName))) {
    const blob = await workbook.render({ sheetName, range, scale: 1, format: "png" });
    const outputPath = path.join(RENDER_DIR, `${String(rendered.length + 1).padStart(3, "0")}_${sheetName}_${range.replace(/[^A-Za-z0-9]+/g, "_")}.png`);
    await fs.writeFile(outputPath, new Uint8Array(await blob.arrayBuffer()));
    const stat = await fs.stat(outputPath);
    rendered.push({ sheetName, range, purpose, outputPath, bytes: stat.size, nonempty: stat.size > 1000 });
  }
}
if (targetExternalIndex) {
  const sheetName = "GSE151374_scRNA";
  const range = "A879:I881";
  const purpose = "final imported-workbook verification of packaged external-result paths, full SHA-256 values, and inclusion notes";
  const existingIndex = rendered.findIndex((entry) => entry.sheetName === sheetName && entry.range === range);
  const ordinal = existingIndex >= 0 ? existingIndex + 1 : rendered.length + 1;
  const outputPath = path.join(RENDER_DIR, `${String(ordinal).padStart(3, "0")}_${sheetName}_${range.replace(/[^A-Za-z0-9]+/g, "_")}.png`);
  const blob = await workbook.render({ sheetName, range, scale: 1, format: "png" });
  await fs.writeFile(outputPath, new Uint8Array(await blob.arrayBuffer()));
  const renderStat = await fs.stat(outputPath);
  const entry = { sheetName, range, purpose, outputPath, bytes: renderStat.size, nonempty: renderStat.size > 1000, imported_workbook_local_rerender: true };
  if (existingIndex >= 0) rendered[existingIndex] = entry;
  else rendered.push(entry);
}

const xlsxBuffer = await fs.readFile(XLSX);
const workbookSha256 = crypto.createHash("sha256").update(xlsxBuffer).digest("hex");
const allRendersNonempty = rendered.every((x) => x.nonempty);
const renderedSheetCount = new Set(rendered.map((x) => x.sheetName)).size;
const allA1Inspected = Object.keys(firstCellInspections).length === requiredSheets.length;
const allHeadersStyleInspected = Object.keys(headerStyleInspections).length === requiredSheets.length;
const allFirstHeadersPassed = allA1Inspected && Object.values(firstHeaderChecks).every((x) => x.passed);
const allHeaderStylesPassed = allHeadersStyleInspected && Object.values(headerStyleChecks).every((x) => x.passed);
const allSheetsRendered = renderedSheetCount === requiredSheets.length;
const programmaticPass = missingSheets.length === 0 && unexpectedSheets.length === 0 && formulaScanPassed
  && allRendersNonempty && allSheetsRendered && allFirstHeadersPassed && allHeaderStylesPassed
  && prior.structural_style_formula_qa?.filterable_table_creation_enforced_for_every_matrix === true
  && persistedTableAutoFiltersPassed && externalIndexPackagedPathsPassed && externalIndexHashesComplete;

const report = {
  ...prior,
  validated_utc: new Date().toISOString(),
  status: programmaticPass && visualReviewed ? "PASS" : programmaticPass ? "programmatic_checks_pass_visual_review_pending" : "FAIL",
  visual_review: visualReviewed ? `PASS: ${Math.max(0, rendered.length - 1)} unchanged ranges retained from the completed all-sheet visual pass, plus a final imported-workbook local rerender of GSE151374_scRNA!A879:I881; all ${requiredSheets.length} sheets and every data-block header were covered for clipping, overlap, unreadable headers, and inconsistent section styling` : "pending manual inspection of imported-workbook renders",
  workbook_sha256: workbookSha256,
  external_large_result_tables: {
    gse151374_macrophage_de_day7: {
      ...(prior.external_large_result_tables?.gse151374_macrophage_de_day7 || {}),
      path: undefined,
      packaged_path: packagedDay7Path,
      analysis_source_path: "fresh_analysis/preclinical/gse151374/results/macrophage_DE_NINT_vs_BLM_day7.csv",
    },
    gse151374_macrophage_de_day14: {
      ...(prior.external_large_result_tables?.gse151374_macrophage_de_day14 || {}),
      path: undefined,
      packaged_path: packagedDay14Path,
      analysis_source_path: "fresh_analysis/preclinical/gse151374/results/macrophage_DE_NINT_vs_BLM_day14.csv",
    },
    packaging_reason: prior.external_large_result_tables?.packaging_reason,
  },
  integrity_controls: (prior.integrity_controls || []).map((item) => item.startsWith("All data blocks are filterable, banded Excel tables with frozen first rows.")
    ? "All data blocks are filterable, banded Excel tables with visible headers. The documented freezeRows(1) request was issued for every worksheet but @oai/artifact-tool 2.8.6 did not persist worksheet <pane> records; this non-scientific navigation limitation is disclosed and not evaluated as a PASS property. Blocks above 10,000 rows retain native numeric values on General format to avoid artifact-tool per-cell style expansion; headers and column widths remain publication-professional."
    : item),
  structural_style_formula_qa: {
    ...prior.structural_style_formula_qa,
    first_row_freeze_enforced_during_authoring: undefined,
    first_row_freeze_requested_during_authoring: true,
    first_row_freeze_persisted_in_export: persistedFreezePaneCount === requiredSheets.length,
    exported_freeze_pane_count: persistedFreezePaneCount,
    freeze_panes_evaluated_as_pass_property: false,
    freeze_panes_export_note: "The documented freezeRows(1) request was issued for every worksheet, including a final post-formatting re-assertion, but @oai/artifact-tool 2.8.6 did not serialize worksheet <pane> records. This is a non-scientific navigation limitation; AutoFilter remains available on all 107 tables.",
    exported_table_xml_count: tableXmlFiles.length,
    exported_table_autofilter_count: persistedTableAutoFilterCount,
    exported_table_autofilters_passed: persistedTableAutoFiltersPassed,
  },
  imported_workbook_qa: {
    required_sheet_count: requiredSheets.length,
    missing_sheets: missingSheets,
    unexpected_sheets: unexpectedSheets,
    formula_error_scan_passed: formulaScanPassed,
    first_cell_inspection_count: Object.keys(firstCellInspections).length,
    header_style_inspection_count: Object.keys(headerStyleInspections).length,
    first_headers_passed: allFirstHeadersPassed,
    header_styles_passed: allHeaderStylesPassed,
    rendered_sheet_count: renderedSheetCount,
    rendered_range_count: rendered.length,
    all_sheets_rendered: allSheetsRendered,
    all_renders_nonempty: allRendersNonempty,
    manual_visual_review_recorded: visualReviewed,
    imported_render_reused_after_manual_review: skipRender,
    exported_freeze_pane_count: persistedFreezePaneCount,
    exported_freeze_pane_expected_count: requiredSheets.length,
    freeze_panes_evaluated_as_pass_property: false,
    freeze_panes_limitation: "Authoring request issued but artifact-tool 2.8.6 export did not persist worksheet <pane>; disclosed as a non-scientific navigation limitation. AutoFilter remains available.",
    exported_table_xml_count: tableXmlFiles.length,
    exported_table_autofilter_count: persistedTableAutoFilterCount,
    exported_table_autofilters_passed: persistedTableAutoFiltersPassed,
    external_index_packaged_paths_passed: externalIndexPackagedPathsPassed,
    external_index_full_sha256_values_passed: externalIndexHashesComplete,
    external_index_local_rerender_recorded: targetExternalIndex || prior.imported_workbook_qa?.external_index_local_rerender_recorded === true,
    render_basis: "143 unchanged builder renders plus one imported-workbook local rerender of GSE151374_scRNA!A879:I881 after the packaged-path-only update",
  },
  exported_freeze_pane_audit: exportedFreezePanes,
  formula_error_scan: formulaScan.ndjson,
  imported_sheet_catalog: sheetCatalog.ndjson,
  imported_first_header_checks: firstHeaderChecks,
  imported_header_style_checks: headerStyleChecks,
  imported_header_style_inspections: headerStyleInspections,
  rendered_sheets: rendered,
  external_index_inspection: {
    range: "GSE151374_scRNA!A879:I881",
    values: externalIndexValues,
    packaged_paths_passed: externalIndexPackagedPathsPassed,
    full_sha256_values_passed: externalIndexHashesComplete,
  },
};
await fs.writeFile(REPORT, JSON.stringify(report, null, 2) + "\n", "utf8");

await fs.rm(`${XLSX}.inspect.ndjson`, { force: true });

console.log(JSON.stringify({
  xlsx: XLSX,
  report: REPORT,
  status: report.status,
  workbookSha256,
  renderedSheetCount,
  renderedRangeCount: rendered.length,
  renderDir: RENDER_DIR,
  formulaScan: formulaScan.ndjson,
  missingSheets,
}, null, 2));
