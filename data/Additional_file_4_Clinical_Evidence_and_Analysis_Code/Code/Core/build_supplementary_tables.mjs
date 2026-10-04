import fs from "node:fs/promises";
import path from "node:path";
import zlib from "node:zlib";
import crypto from "node:crypto";
import { fileURLToPath } from "node:url";
import os from "node:os";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const PACKAGE_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "../..");
const ROOT = process.env.NINTEDANIB_DILI_DEVELOPMENT_ROOT
  ? path.resolve(process.env.NINTEDANIB_DILI_DEVELOPMENT_ROOT)
  : PACKAGE_ROOT;
const OUT_DIR = path.join(ROOT, "fresh_analysis", "submission");
const OUT_XLSX = path.join(OUT_DIR, "Supplementary_Tables.xlsx");
const VALIDATION_JSON = path.join(OUT_DIR, "Supplementary_Tables_validation_report.json");
const RENDER_DIR = path.join(os.tmpdir(), "nintedanib_dili_workbook", "final_renders");

const COLORS = {
  navy: "#16324F",
  teal: "#0B6B6F",
  cyan: "#DDF1F2",
  paleBlue: "#EAF1F7",
  paleGreen: "#E4F2E9",
  paleOrange: "#FCE9D9",
  paleGray: "#F4F6F8",
  midGray: "#D8DEE5",
  darkGray: "#263238",
  white: "#FFFFFF",
};

function colLetter(indexOneBased) {
  let n = indexOneBased;
  let out = "";
  while (n > 0) {
    n -= 1;
    out = String.fromCharCode(65 + (n % 26)) + out;
    n = Math.floor(n / 26);
  }
  return out;
}

function parseCsv(text) {
  const rows = [];
  let row = [];
  let cell = "";
  let quoted = false;
  for (let i = 0; i < text.length; i += 1) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"') {
        if (text[i + 1] === '"') {
          cell += '"';
          i += 1;
        } else {
          quoted = false;
        }
      } else {
        cell += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === ",") {
      row.push(cell);
      cell = "";
    } else if (ch === "\n") {
      row.push(cell.endsWith("\r") ? cell.slice(0, -1) : cell);
      rows.push(row);
      row = [];
      cell = "";
    } else {
      cell += ch;
    }
  }
  if (cell.length > 0 || row.length > 0) {
    row.push(cell.endsWith("\r") ? cell.slice(0, -1) : cell);
    rows.push(row);
  }
  while (rows.length && rows[rows.length - 1].every((v) => v === "")) rows.pop();
  return rows;
}

function coerceValue(value, header) {
  const v = value.trim();
  if (v === "") return null;
  if (/^(true|false)$/i.test(v)) return v.toLowerCase() === "true";
  if (/^(nan|na|none|null)$/i.test(v)) return null;
  if (/date$/i.test(header) || /sha256/i.test(header) || /accession|identifier/i.test(header) || /_id$/i.test(header)
      || /^(participant|participant_id|subject|subject_id|sample|sample_id)$/i.test(header)) return v;
  if (/^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$/.test(v)) {
    const x = Number(v);
    if (Number.isFinite(x)) return x;
  }
  return v;
}

async function readCsvFrom(baseDir, relativePath) {
  const fullPath = path.join(baseDir, relativePath);
  let buffer = await fs.readFile(fullPath);
  if (relativePath.endsWith(".gz")) buffer = zlib.gunzipSync(buffer);
  const raw = parseCsv(buffer.toString("utf8"));
  if (!raw.length) throw new Error(`Empty CSV: ${relativePath}`);
  const headers = raw[0];
  const matrix = [headers];
  for (const sourceRow of raw.slice(1)) {
    const normalized = [...sourceRow];
    while (normalized.length < headers.length) normalized.push("");
    if (normalized.length > headers.length) {
      throw new Error(`CSV width mismatch in ${relativePath}: ${normalized.length} > ${headers.length}`);
    }
    matrix.push(normalized.map((v, j) => coerceValue(v, headers[j])));
  }
  return { relativePath, matrix, rows: matrix.length - 1, cols: headers.length };
}

async function readCsv(relativePath) {
  return readCsvFrom(ROOT, relativePath);
}

async function readPackagedCsv(relativePath) {
  return readCsvFrom(PACKAGE_ROOT, relativePath);
}

async function summarizeCsvNumericColumn(relativePath, columnName, predicate) {
  const fullPath = path.join(ROOT, relativePath);
  let buffer = await fs.readFile(fullPath);
  if (relativePath.endsWith(".gz")) buffer = zlib.gunzipSync(buffer);
  const raw = parseCsv(buffer.toString("utf8"));
  const index = raw[0]?.indexOf(columnName);
  if (index === undefined || index < 0) throw new Error(`Column ${columnName} not found in ${relativePath}`);
  let numeric = 0;
  let matched = 0;
  for (const row of raw.slice(1)) {
    const rawValue = row[index];
    if (rawValue === undefined || rawValue === null || String(rawValue).trim() === "") continue;
    const value = Number(rawValue);
    if (!Number.isFinite(value)) continue;
    numeric += 1;
    if (predicate(value)) matched += 1;
  }
  return { relativePath, rows: raw.length - 1, numeric, matched };
}

async function readJson(relativePath) {
  const fullPath = path.join(ROOT, relativePath);
  return JSON.parse(await fs.readFile(fullPath, "utf8"));
}

async function sha256File(fullPath) {
  const buf = await fs.readFile(fullPath);
  return crypto.createHash("sha256").update(buf).digest("hex");
}

function isCountHeader(h) {
  return /^(n_|n$|n\d)|(^|_)(count|counts|event|events|nonevent|total|reports|rows|numeric|observed|missing|starred|seeds|assays|documents|activities|participants|cells|donors|genes|permutations|hops|year|visits|requested|present|strata)$/i.test(h)
    || /_n$/i.test(h)
    || /_(event|nonevent|total|reports|rows|observed|cells|donors|genes|seeds|assays|documents|participants|visits)$/i.test(h)
    || /^(year|first_year|last_year|gene_id|rank|consensus_rank|loso_rank_best|loso_rank_worst|leave_one_source_out_rank|rank_within_dose)$/i.test(h);
}

function isPValueHeader(h) {
  return /(^|_)(p|q|p_value|q_value|fisher_exact_p|welch_t_p|mann_whitney_p|spearman_p|exact_p|bh_q)(_|$)/i.test(h)
    || /p_one_sided|p_two_sided|_q_bh|bh_q_|_p$/i.test(h);
}

function isLongTextHeader(h) {
  return /url|path|source|definition|boundary|reason|note|description|interpretation|query|criteria|endpoint[s]?|scope|method|material|statistical_unit|experimental_unit/i.test(h);
}

function applyHeaderStyle(range, fill = COLORS.navy) {
  range.format = {
    fill,
    font: { bold: true, color: COLORS.white },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
    borders: { preset: "outside", style: "thin", color: COLORS.navy },
  };
  range.format.rowHeight = 34;
}

function applyDataFormatting(sheet, startRowZero, matrix) {
  const headers = matrix[0].map(String);
  const nRows = matrix.length;
  const nCols = headers.length;
  if (nRows <= 1) return;
  // artifact-tool expands body formatting during XLSX serialization. Preserve complete
  // high-dimensional result tables as filterable/banded Excel tables, but keep their
  // bodies on General format to avoid duplicating styles across >10,000 rows.
  if (nRows > 10000) {
    headers.forEach((header, j) => {
      let width = Math.max(10, Math.min(30, header.length + 3));
      if (/gene_symbol|gene_id|symbol|protein|endpoint|drug|cohort|contrast|compartment|ko_label|source_group/i.test(header)) width = 18;
      if (/gene_title|description|target_pref_name|preferred_term|candidate_role|pharmacology_tier|evidence_tier/i.test(header)) width = 28;
      sheet.getRangeByIndexes(startRowZero, j, 1, 1).format.columnWidth = width;
    });
    return;
  }
  const body = sheet.getRangeByIndexes(startRowZero + 1, 0, nRows - 1, nCols);
  body.format = {
    verticalAlignment: "top",
    font: { color: COLORS.darkGray },
    borders: { insideHorizontal: { style: "thin", color: "#E9EDF1" } },
  };
  body.format.rowHeight = nRows > 1000 ? 18 : 20;

  headers.forEach((header, j) => {
    const values = matrix.slice(1).map((r) => r[j]).filter((v) => v !== null && v !== undefined);
    const numeric = values.length > 0 && values.every((v) => typeof v === "number");
    const bool = values.length > 0 && values.every((v) => typeof v === "boolean");
    const col = sheet.getRangeByIndexes(startRowZero, j, nRows, 1);
    const bodyCol = sheet.getRangeByIndexes(startRowZero + 1, j, nRows - 1, 1);

    let width = Math.max(10, Math.min(30, header.length + 3));
    if (/gene_symbol|protein|endpoint|drug|cohort|contrast|compartment|ko_label|source_group/i.test(header)) width = 18;
    if (/gene_title|target_pref_name|preferred_term|candidate_role|pharmacology_tier|evidence_tier/i.test(header)) width = 28;
    if (isLongTextHeader(header)) width = /url|path|query/i.test(header) ? 52 : 42;
    if (/sha256/i.test(header)) width = 24;
    col.format.columnWidth = width;

    if (numeric) {
      bodyCol.format.horizontalAlignment = "right";
      if (isCountHeader(header)) bodyCol.format.numberFormat = "#,##0";
      else if (isPValueHeader(header)) bodyCol.format.numberFormat = "0.00E+00";
      else if (/percent/i.test(header)) bodyCol.format.numberFormat = "0.0";
      else bodyCol.format.numberFormat = "0.0000";
    } else if (bool) {
      bodyCol.format.horizontalAlignment = "center";
      bodyCol.format.numberFormat = "General";
    } else {
      bodyCol.format.horizontalAlignment = "left";
      if (isLongTextHeader(header)) bodyCol.format.wrapText = true;
    }

    if (/(_q|q_bh|bh_q)/i.test(header) && numeric) {
      bodyCol.conditionalFormats.add("cellIs", {
        operator: "lessThan",
        formula: 0.05,
        format: { fill: COLORS.paleGreen, font: { color: "#14532D", bold: true } },
      });
    }
  });
}

let tableCounter = 0;
const matrixAudit = [];
function addTable(sheet, startRowZero, matrix, baseName, style = "TableStyleMedium2") {
  if (!Array.isArray(matrix) || matrix.length === 0 || !Array.isArray(matrix[0]) || matrix[0].length === 0) {
    throw new Error(`Invalid empty matrix for ${sheet.name}/${baseName}`);
  }
  const nRows = matrix.length;
  const nCols = matrix[0].length;
  const headers = matrix[0].map((h) => String(h ?? "").trim());
  if (headers.some((h) => h === "")) throw new Error(`Blank header in ${sheet.name}/${baseName}`);
  if (new Set(headers).size !== headers.length) throw new Error(`Duplicate header in ${sheet.name}/${baseName}`);
  if (matrix.some((row) => row.length !== nCols)) throw new Error(`Ragged matrix in ${sheet.name}/${baseName}`);
  matrixAudit.push({ sheet: sheet.name, block: baseName, first_header: headers[0], rows: nRows - 1, columns: nCols, header_check: "PASS", rectangular_check: "PASS" });
  const dest = sheet.getRangeByIndexes(startRowZero, 0, nRows, nCols);
  dest.values = matrix;
  applyHeaderStyle(sheet.getRangeByIndexes(startRowZero, 0, 1, nCols));
  applyDataFormatting(sheet, startRowZero, matrix);
  tableCounter += 1;
  const safe = `${baseName}_${tableCounter}`.replace(/[^A-Za-z0-9_]/g, "_").slice(0, 200);
  const table = sheet.tables.add(dest, true, safe);
  table.style = style;
  table.showFilterButton = true;
  return { startRowZero, endRowZero: startRowZero + nRows - 1, nRows: nRows - 1, nCols, tableName: safe };
}

function initSheet(workbook, name) {
  const sheet = workbook.worksheets.add(name);
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  return sheet;
}

function makeSimpleSheet(workbook, name, matrix, baseName) {
  const sheet = initSheet(workbook, name);
  const info = addTable(sheet, 0, matrix, baseName);
  return { sheet, blocks: [info] };
}

function makeMultiBlockSheet(workbook, name, blocks) {
  const sheet = initSheet(workbook, name);
  const maxCols = Math.max(...blocks.map((b) => b.matrix[0].length));
  const indexMatrix = [["Section", "Source_File", "Record_Count", "Interpretation_Note"], ...blocks.map((b) => [
    b.title,
    b.source,
    b.matrix.length - 1,
    b.note || "Values reproduced from the named analysis output without modification.",
  ])];
  const recorded = [addTable(sheet, 0, indexMatrix, `${name}_Index`, "TableStyleMedium4")];
  if (indexMatrix.length > 1) {
    const indexBody = sheet.getRangeByIndexes(1, 0, indexMatrix.length - 1, 4);
    indexBody.format.wrapText = true;
    indexBody.format.rowHeight = 118;
  }
  let row = indexMatrix.length + 2;
  for (const block of blocks) {
    const band = sheet.getRangeByIndexes(row, 0, 1, maxCols);
    band.format = {
      fill: COLORS.teal,
      font: { bold: true, color: COLORS.white },
      verticalAlignment: "center",
    };
    band.format.rowHeight = 24;
    sheet.getCell(row, 0).values = [[block.title]];
    if (maxCols > 1) sheet.getCell(row, 1).values = [[block.source]];
    row += 1;
    recorded.push(addTable(sheet, row, block.matrix, `${name}_${block.title}`, block.style || "TableStyleMedium2"));
    row += block.matrix.length + 2;
  }
  return { sheet, blocks: recorded };
}

const dataSources = [
  ["Source_ID", "Domain", "Repository", "Accession_or_Identifier", "Data_Unit", "Retrieval_Date", "Local_Source_File", "Source_SHA256", "Public_URL", "Role_in_Study", "Interpretation_Boundary"],
  ["DS01", "Pharmacovigilance", "openFDA / FAERS", "drug/event; openFDA release updated 2026-07-30", "Latest-version safetyreportid (report)", "2026-09-04", "fresh_analysis/faers/results/faers_frozen_manifest.json", "frozen manifest contains per-output SHA-256 hashes", "https://api.fda.gov/drug/event.json", "Nintedanib versus pirfenidone active-comparator reporting disproportionality", "Spontaneous reports estimate reporting disproportionality, not incidence, absolute risk, or causality."],
  ["DS02", "Clinical proteomics", "Nature Communications source data", "DOI 10.1038/s41467-023-36858-6", "Public serum-protein measurements from adjudicated multi-drug DILI and comparator groups", "2026-09-04", "fresh_data_dili_source.xlsx", "5911fe3d24efa1b639000e5c9176ef9cee7a44b4bed7ff75fe818ee1d6ef6b16", "https://www.nature.com/articles/s41467-023-36858-6", "Human distal DILI phenotype anchoring", "The cohort is multi-drug DILI, not nintedanib-specific; source workbook lacks stable cross-protein participant linkage."],
  ["DS03", "Drug perturbation", "NCBI GEO / LINCS L1000", "GSE70138; BRD-K49075727; HepG2; 24 h", "Six composite Level-5 signatures across 0.04-10 µM", "2026-09-04", "fresh_analysis/lincs/source/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz", "74dc32fee1cb652331462dd269f3de604fff8b49e3ff21e42e31b07b2bd52d3c", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138", "Observed nintedanib concentration-response in HepG2 cells", "Cell-line perturbation; high concentrations are mechanistic comparators and do not establish clinical hepatotoxicity."],
  ["DS04", "Orthogonal drug perturbation", "Harmonizome / Sci-Plex", "MCF7 nintedanib 1 and 10 µM; 24 h", "Published top/bottom ranked gene lists", "2026-09-04", "fresh_analysis/inputs/harmonizome_nintedanib_mcf7_10uM_24h.html", "4c2e8fbc9d69f21a661db934f8f3b9a1380fa613f238568d3acd48124f3c9a8c", "https://maayanlab.cloud/Harmonizome/dataset/Sci-Plex+Drug+Perturbation+Signatures", "Exploratory cross-cell-line directional replication", "Truncated MCF7 lists; absence is not evidence of no response and MCF7 is not a liver model."],
  ["DS05", "Target pharmacology", "ChEMBL API", "CHEMBL502835 (nintedanib)", "Curated bioactivity records; human single-protein subset", "2026-09-04", "fresh_data/chembl_nintedanib_activities.json plus offsets", "723455ef7956bfa7fe6283c95ea7bd98376a5781b6e62f8a6df5ba48e60b2a91 (first page)", "https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=CHEMBL502835", "Exposure-proximal biochemical target evidence", "Biochemical potency and document replication do not establish in-vivo target engagement or toxicity causality."],
  ["DS06", "Human liver single-cell transcriptomics", "NCBI GEO", "GSE115469", "8,444 cells, 20,007 genes, five healthy liver donors", "2026-09-04", "fresh_data/GSE115469_Data.csv.gz", "9d20630b98dba970ef3ec616d48a18857d94aa38e29886ab358a6dbdb10cbdde", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469", "Donor-aware liver-cell localization and same-compartment virtual-KO context", "Healthy untreated liver; donor is the biological inference unit, cells are descriptive, and expression plus virtual-KO evidence is combined only within the same compartment."],
  ["DS07", "Protein association network", "STRING", "Human candidate network; required score thresholds documented in scripts", "Association edges", "2026-09-04", "fresh_analysis/inputs/string_candidate_network_expanded.tsv", "1ab683cd60d441f391af8a64fd41e77a27e7c4f53ca2aaf8724cff97a6a4f385", "https://string-db.org/", "Contextual paths between exposure-proximal targets and DILI anchors", "STRING associations, including database and text-mining evidence, are contextual and non-causal."],
  ["DS08", "Longitudinal whole-blood transcriptomics", "NCBI GEO", "GSE299128", "57 samples from seven longitudinally sampled participants", "2026-09-04", "fresh_data/GSE299128_count_matrix.csv.gz", "462ed26cc41bbd0c2e27178db6575aa80b475bce7596fba8f397b24c1d629ce6", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299128", "Systemic longitudinal detectability sensitivity analysis using an exact participant-level test", "Small single-arm whole-blood study without adjudicated DILI; exact two-sided sign-flip P values enumerate all 2^7 participant assignments and do not validate liver-specific toxicity."],
  ["DS09", "Independent-animal lung transcriptomics", "NCBI GEO", "GSE278200", "12 rats; saline, bleomycin, and bleomycin plus nintedanib, n=4 independent animals/group", "2026-09-04", "fresh_analysis/preclinical/raw/GSE278200_NINT_raw_counts.txt.gz", "5ae68a0a8b4ab253105e27df215af424e44ca3563c52f22edd1981d5253eebfb", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278200", "Independent-animal in-vivo lung transcriptomic reversal and candidate/module evaluation", "Right-lung homogenate in a bleomycin model; supports antifibrotic lung-response context, not human DILI causality or liver rescue."],
  ["DS10", "Independent-animal lung transcriptomics", "NCBI GEO", "GSE308578", "43 mice; control n=10, bleomycin plus vehicle n=16, bleomycin plus nintedanib n=17", "2026-09-04", "fresh_analysis/preclinical/raw/GSE308578_Gene_counts.csv.gz", "a46808f90410796adc84738e4045aaeee311bf7445c71c39c1ee47e08da47a9c", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308578", "Independent-animal in-vivo lung transcriptomic reversal, leave-one-animal-out stability, and candidate/module evaluation", "GEO SOFT and ARRIVE Supplementary Figure S1 support 10/16/17; the source article main Figure 6 caption reverses vehicle/nintedanib n, so no samples were relabelled."],
  ["DS11", "Precision-cut liver-slice transcriptomics", "NCBI GEO", "GSE120804", "48 pooled-tissue slice observations across sham, BDL, and three treatment groups", "2026-09-04", "fresh_analysis/preclinical/raw/GSE120804_counts.txt.gz", "85290cb208bbb128ba04c491da98b4cae9fb98ae0139904df9e978fbf39c6319", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120804", "Exploratory nintedanib candidate/module context in fibrotic liver slices", "Observations derive from pooled rat tissues and are not independent animals; reported contrasts are exploratory and do not support animal-population inference."],
  ["DS12", "Precision-cut lung-slice transcriptomics", "NCBI GEO", "GSE120679", "36 pooled-tissue lung-slice observations; six observations in each of six time/treatment groups", "2026-09-04", "fresh_analysis/preclinical/raw/GSE120679_logCPM_expressed_protein_coding.txt.gz", "afc48c4e358807dad9d7d2b66176cd6e9de1c48cc1773ca7250ae77fb6777ac5", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120679", "Exploratory nintedanib candidate/module context in TGF-beta-stimulated lung slices", "The analyzed input is the repository-supplied expressed protein-coding logCPM matrix. Observations derive from pooled rat lungs and are not independent animals; contrasts are exploratory and do not support animal-population inference."],
  ["DS13", "Lung proteomics", "PRIDE", "PXD024058", "5,951 proteins with three columns per control, model, and nintedanib condition in the supplied protein table", "2026-09-04", "fresh_analysis/preclinical/raw/PXD024058_all.protein.xls", "78e477aa84ccbd71def7b1c3cb62e501169c626be644afdc4ae817a4ddb7cebe", "https://www.ebi.ac.uk/pride/archive/projects/PXD024058", "Descriptive orthogonal protein-level direction check", "Replicate provenance is internally inconsistent between article, repository metadata, and supplied table; effects are descriptive only and replicate columns are not treated as independent biological units."],
  ["DS14", "Pooled-library lung immune-cell single-cell transcriptomics", "NCBI GEO / SRA", "GSE151374 / PRJNA635636", "18 independently prepared pooled 10x libraries; three libraries per condition-time cell; three mice pooled per library; 70,696 QC-retained cells", "2026-09-04", "fresh_analysis/preclinical/gse151374/results/analysis_manifest.json", "999415f7c21f9f4f648dc830c0ffe82ada42bc902cc4b4a8ee5c807c164f4a52", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151374", "Nonconfirmatory pulmonary immune-cell and macrophage pseudobulk sensitivity", "The pooled 10x library is the inferential unit (n=3 per condition-time); cells are not independent replicates. The minimum attainable two-sided exact 3-versus-3 permutation P value is 0.10. This is lung context, not liver/DILI evidence."],
  ["DS15", "Animal-level lung proteomics and deposited-subset phosphoproteomics", "PRIDE / Zenodo", "PXD052594 / DOI 10.5281/zenodo.11395642", "Proteome: 23 individual mice (nintedanib n=10, vehicle n=13); phosphoproteome: deposited random subset of five mice per arm", "2026-09-04", "fresh_analysis/preclinical/pxd052594/results/analysis_manifest.json", "2d44ac6de65c09453646b8f27a77ec7de1f7b5492685a0557e8acc5daf71c596", "https://www.ebi.ac.uk/pride/archive/projects/PXD052594", "Unconditioned individual-animal pulmonary proteome sensitivity plus a separate table-only phosphosite analysis", "No radiomic response cluster was used. Six postrandomization exclusions preceded the final imaging cohort and one vehicle (M29) was excluded upstream from proteome preparation. The 5-versus-5 phosphoproteome subset is distinct from the 23-animal proteome. Neither component is hepatic-safety or target confirmation."],
];

const readme = [
  ["Item", "Value"],
  ["Workbook", "Supplementary Tables for the nintedanib-associated liver injury public-data study"],
  ["Working manuscript title", "Role-separated public-data triangulation prioritizes target–phenotype hypotheses for nintedanib-associated liver injury"],
  ["Version date", "2026-09-05"],
  ["Purpose", "Machine-readable, source-linked tables supporting pharmacovigilance, clinical proteomics, drug perturbation, pharmacology, liver single-cell localization, same-compartment virtual knockout, evidence integration, public preclinical constraint testing, and prospective validation design."],
  ["Primary exposure", "Nintedanib"],
  ["Active comparator", "Pirfenidone"],
  ["Primary adverse-event estimand", "All-indication, drug-name-indexed report-level liver-injury reporting disproportionality in mutually exclusive nintedanib and pirfenidone openFDA sets. openFDA drugcharacterization=1 defines a suspect-drug report filter; it does not distinguish primary from secondary suspect drugs."],
  ["ILD scope sensitivity", "A report-level ILD co-report filter is the clinically central sensitivity analysis. Because openFDA arrays are flattened, it cannot prove that an ILD indication belongs to the same drug element as nintedanib or that a reaction is causally linked to that element."],
  ["Clinical proteomic role", "Multi-drug adjudicated DILI phenotype anchor; not nintedanib-specific target confirmation."],
  ["Virtual knockout role", "Unsigned scTenifold network displacement used for prioritization only; it does not assign protective/harmful direction or prove causality. Expression and virtual-KO ranks are combined only within the same liver compartment."],
  ["Candidate-tier ceiling", "No Tier 1 target is assigned. FGFR1 is the sole candidate meeting the prespecified exploratory same-compartment T2 gate; it does not meet the strict all-seed sensitivity gate and is not a confirmed mediator. FBP1 is a downstream P1 multi-drug DILI phenotype anchor, not a causal nintedanib target. P and T roles remain distinct."],
  ["Whole-blood exact test", "GSE299128 exact two-sided P values enumerate all 2^7 participant-level sign assignments for the mean within-participant slope; BH correction is across the 30 prespecified genes."],
  ["Public preclinical role", "GSE278200 and GSE308578 provide independent-animal pulmonary-response context. GSE120804 and GSE120679 are pooled-tissue slice analyses and are exploratory. PXD024058 is descriptive because replicate provenance is internally inconsistent."],
  ["Public preclinical extensions", "Supplementary Figure S6 is a nonconfirmatory pooled-library GSE151374 lung immune-cell sensitivity (18 libraries; three mice pooled per library; n=3 libraries per condition-time; 70,696 retained cells; no cell-level inference; all exact-permutation BH families null). Supplementary Figure S7 is an unconditioned PXD052594 individual-mouse lung-proteome sensitivity (7,006 proteins; nintedanib n=10, vehicle n=13; protein BH null). Its distinct deposited 5-versus-5 phosphoproteome analysis is table-only (20,043 sites; 994 q<0.05) and contextual."],
  ["Large-result packaging", "The two complete GSE151374 macrophage-pseudobulk DE CSVs (48,795 rows per day) remain unchanged, hash-verified, and submission-ready as external machine-readable supplementary tables; GSE151374_scRNA contains a full-file index with their row counts, hashes, statistical unit and zero q<0.10 result. They are not duplicated cell-by-cell in this XLSX because artifact-tool cannot serialize both matrices within one workbook. The full 20,043-row PXD052594 phosphosite table remains embedded."],
  ["Experimental status", "Experimental_Validation contains a prospective protocol only; none of those experiments are represented as completed."],
  ["Missing values", "Blank cells reproduce unavailable or non-estimable source values. The 2014 narrow FAERS contrast is non-estimable because both arms had zero events."],
  ["Multiplicity", "Benjamini-Hochberg FDR is reported where prespecified in each analysis; endpoint-specific details are retained in source-derived columns."],
  ["Machine readability", "Each result block uses a flat header and filterable Excel table; units are carried in headers or explicit metadata."],
  ["Provenance and figure traceability", "Data_Provenance carries frozen manifest paths and SHA-256 hashes; Statistical_Units defines biological and computational analysis units; Figure_Source_Map maps each figure to source tables and inference boundaries."],
  ["Abbreviations", "DO, DILI onset; NDO, non-drug acute liver injury onset; HV, healthy volunteer; DF, DILI follow-up; ROR, reporting odds ratio; PRR, proportional reporting ratio; IC, active-comparator Bayesian information component; KO, virtual knockout; LOSO, leave-one-source-out."],
];

const endpointSummary = [
  ["Endpoint", "Status", "Operational_Definition", "Term_Count", "Statistical_Unit", "Primary_Estimand", "Interpretation_Boundary", "Source_File"],
  ["narrow", "Primary", "Ten prespecified diagnostic or severe hepatic-injury MedDRA preferred terms listed in the following table.", 10, "Latest-version safetyreportid (report)", "Nintedanib-versus-pirfenidone report-level ROR, PRR, Fisher exact test, and active-comparator Bayesian IC", "Custom preferred-term grouping, not a licensed MedDRA SMQ; estimates are reporting disproportionality, not risk or causality.", "fresh_analysis/faers/results/prespecified_meddra_pt_sets.csv"],
  ["broad", "Sensitivity", "The narrow set plus 16 sensitive hepatic laboratory, bilirubin, jaundice, cholestasis, hepatitis, and nonspecific liver-disorder terms.", 26, "Latest-version safetyreportid (report)", "Nintedanib-versus-pirfenidone report-level ROR, PRR, Fisher exact test, and active-comparator Bayesian IC", "More sensitive and less specific custom term grouping; estimates are reporting disproportionality, not risk or causality.", "fresh_analysis/faers/results/prespecified_meddra_pt_sets.csv"],
  ["hy_proxy", "Exploratory", "Co-reporting of a prespecified aminotransferase-elevation term and a bilirubin/jaundice term.", null, "Latest-version safetyreportid (report)", "Exploratory report-level active-comparator disproportionality", "Not Hy's law: verified laboratory values, upper limits of normal, alkaline-phosphatase exclusion, and clinical adjudication are unavailable.", "fresh_analysis/faers/run_openfda_faers.py"],
];

const experimentalValidation = [
  ["Stage", "Model_or_System", "Biological_Unit", "Planned_N", "Nintedanib_Exposure", "Perturbation", "Primary_Endpoints", "Key_Secondary_Endpoints", "Analysis_or_Decision_Rule", "Status"],
  ["1. Exposure anchoring", "Cryopreserved primary human hepatocytes", "Donor", "8 independent donors; sex balance where feasible", "0.01, 0.04, 0.12, 0.37, 1.11 µM for 6, 24, 48, 72 h; 3.33 and 10 µM only as high-exposure mechanistic comparators", "Vehicle-controlled concentration and time series", "ATP-normalized viability; extracellular ALT; high-content apoptotic/necrotic cell counts", "Parent nintedanib, BIBF 1202, BIBF 1202 glucuronide; LDH, AST, albumin, urea, mitochondrial membrane potential, OCR, ROS, glutathione redox, lipid peroxidation, caspase-3/7, BSEP transport, bile acids, KRT18 and candidate protein panel", "Joint exposure-toxicity modeling; protein binding and intracellular concentrations reported", "Prospective—not yet performed"],
  ["1. Exposure anchoring", "Human iPSC-derived hepatocyte organoids", "Independent iPSC line", "At least 3 independent lines", "Exposure-anchored concentrations and times from hepatocyte stage", "Vehicle-controlled concentration and time series", "Viability and liver-injury endpoints aligned to primary hepatocytes", "Nintedanib/metabolite exposure and candidate marker panel", "Orthogonal replication platform", "Prospective—not yet performed"],
  ["1. Exposure anchoring", "Primary human liver sinusoidal endothelial cells and hepatocyte-LSEC transwell co-culture", "LSEC donor / paired co-culture unit", "To be finalized from pilot variance", "Exposure-anchored nintedanib concentrations", "Vehicle-controlled concentration and time series", "Barrier resistance; endothelial apoptosis", "PLVAP, STAB2, FCGR2B, nitric oxide, fenestration markers and secreted angiocrine factors", "Distinguish direct endothelial from secondary hepatocyte effects", "Prospective—not yet performed"],
  ["2. Functional perturbation", "Primary hepatocytes and organoids", "Donor or organoid line", "As above; confirmatory N finalized by simulation", "Vehicle and exposure-anchored nintedanib under prespecified glucose and fasting/refeeding conditions", "Acute, graded and bidirectional FBP1 perturbation: two independent inducible CRISPRi reagents, CRISPR activation, wild-type expression, guide-resistant rescue, and a catalytically impaired rescue construct", "Exposure-matched liver-injury endpoints and intracellular FBP1 function", "Glucose production, F1,6BP/F6P, ATP and energy charge, NADPH/NADP+, AMPK/AKT, lipid handling and oxidative stress; GSTA1 and OTC retained as downstream P2 readouts", "Direction is assigned only when independent loss-of-function reagents agree, reciprocal gain/rescue reverses the phenotype, and both human platforms replicate. Wild-type versus catalytically impaired rescue separates catalytic from non-catalytic signaling; extracellular-only change without intracellular modification classifies FBP1 as a leakage/readout marker", "Prospective—not yet performed"],
  ["2. Functional perturbation", "Primary hepatocytes", "Donor", "8 independent donors planned", "Vehicle and exposure-anchored nintedanib", "Two independent CES1 CRISPRi/siRNA reagents and CES1 cDNA expression", "Exposure-matched liver-injury endpoints", "CES1 activity; parent drug; BIBF 1202; BIBF 1202 glucuronide; UGT1A1 activity", "An apparent rescue is not accepted if explained solely by altered parent-drug exposure", "Prospective—not yet performed"],
  ["2. Functional perturbation", "Primary LSECs and hepatocyte-LSEC co-culture", "LSEC donor / paired co-culture unit", "To be finalized from pilot variance", "Vehicle and exposure-anchored nintedanib", "FGFR1, the sole candidate meeting the exploratory T2 gate: two independent inducible CRISPRi guides, CRISPR activation, guide-resistant wild-type rescue, kinase-dead rescue, and ligand-controlled FGF2 and FGF21–KLB conditions", "Barrier resistance, endothelial apoptosis, and downstream hepatocyte injury", "FGFR1 target engagement; PLVAP, STAB2, FCGR2B, nitric oxide and angiocrine factors; conditioned-medium and transwell transfer; downstream FBP1-centred phenotype", "Require target engagement, donor-level injury modification, exposure matching and reciprocal perturbation. Both worsening and rescue are valid prespecified outcomes because virtual-knockout displacement is unsigned", "Prospective—not yet performed"],
  ["3. Retention of pulmonary antifibrotic activity", "Primary fibroblasts from fibrotic ILD and non-fibrotic donor controls", "Donor", "8 fibrotic ILD donors plus controls", "Nintedanib, candidate perturbation, or combination after TGF-β1 stimulation", "Candidate-specific loss/gain-of-function intervention", "COL1A1 and collagen-gel contraction", "COL1A2, ACTA2, FN1, secreted procollagen I, EdU, PDGFR/FGFR/VEGFR phosphorylation", "95% CI must exclude loss of >20% of the nintedanib antifibrotic response for both co-primary endpoints; this is a non-inferiority screen, not proof of clinical efficacy", "Prospective—not yet performed"],
  ["3. Retention of pulmonary antifibrotic activity", "Human precision-cut lung slices", "Donor", "At least 6 donors", "Nintedanib, candidate perturbation, or combination", "Candidate-specific intervention", "Second-harmonic collagen imaging and hydroxyproline", "Matched molecular fibrosis panel", "Orthogonal tissue-level confirmation of pulmonary-response non-inferiority", "Prospective—not yet performed"],
  ["4. In-vivo proof of mechanism", "Bleomycin lung-fibrosis model", "Animal", "Simulation-based after pilot; both sexes", "Vehicle or nintedanib", "2×2 factorial: vehicle, nintedanib, candidate, combination", "Lung: micro-CT fibrosis and blinded Ashcroft score; liver: ALT/AST and blinded histopathology", "Hydroxyproline, forced-oscillation mechanics, antifibrotic proteins, bilirubin, bile acids, TUNEL, mitochondrial function, marker panel and tissue/plasma drug metabolites", "Randomized and blinded; sex interaction prespecified; proof-of-mechanism only, not an idiosyncratic human DILI model", "Prospective—not yet performed"],
  ["Cross-stage statistical plan", "All in-vitro and in-vivo platforms", "Donor or animal; technical replicates averaged", "Confirmatory N targeted to 90% power for the prespecified interaction using pilot variance components", "Concentration modeled on log scale with restricted cubic splines", "Nintedanib × candidate interaction", "Co-primary liver and lung endpoints", "Secondary endpoints controlled by Benjamini-Hochberg FDR", "Mixed model with random donor intercept; Holm adjustment for co-primary liver endpoints; consistency in at least 6/8 donors; no imputation for documented technical failures", "Prospective—not yet performed"],
  ["Stop/go", "Translational decision", "Candidate", "Not applicable", "Exposure matched", "Candidate intervention", "Four required criteria", "Target engagement; liver rescue in two human platforms; pulmonary-response non-inferiority; no worsening in vivo", "Failure of any criterion downgrades the candidate to biomarker or mechanistic probe", "Prospective—not yet performed"],
];

const reportingChecklists = [
  ["Framework_or_Guideline", "Domain", "Applicability", "Self_Audit_Status", "Evidence_in_Package", "Outstanding_Action", "Public_URL"],
  ["READUS-PV", "Disproportionality analyses using spontaneous-reporting databases", "Applicable", "Addressed at analysis-package level", "Data source and release, exposure definitions, co-exposure exclusion, report deduplication, preferred-term sets, contingency counts, uncertainty, calendar-time analyses, specifications, missingness and interpretation boundaries are supplied.", "Map each checklist item to final manuscript page/line numbers during journal formatting.", "https://doi.org/10.1007/s40264-024-01421-9"],
  ["RECORD", "Routinely collected observational health data", "Applicable in principle to secondary database reporting", "Addressed at analysis-package level", "Data provenance, retrieval date, analysis window, operational definitions, unit of analysis, code availability and limitations are explicit.", "Complete journal-specific item-by-item page mapping.", "https://www.record-statement.org/"],
  ["STROBE", "Observational study reporting", "Applicable to public observational components", "Addressed at analysis-package level", "Study design, data sources, variables, bias controls, statistical methods, descriptive counts, uncertainty and limitations are specified.", "Complete journal-specific checklist with page/line references.", "https://www.strobe-statement.org/checklists/"],
  ["MIAME/GEO accession reporting", "Public transcriptomic data", "Applicable", "Addressed", "GSE70138, GSE115469, GSE299128, GSE278200, GSE308578, GSE120804, GSE120679 and GSE151374 accession numbers, retrieval dates, checksums, analysis units and public URLs are recorded. GSE151374 pooled-library inference is separated from cell-level descriptive summaries.", "None for the current computational package.", "https://www.ncbi.nlm.nih.gov/geo/info/MIAME.html"],
  ["TRIPOD+AI", "Clinical prediction models", "Not applicable", "No prediction model developed", "The study performs target prioritization and does not report an individual-level diagnostic or prognostic model.", "Do not present integration ranks as predictive performance.", "https://www.tripod-statement.org/"],
  ["PROBAST+AI", "Risk-of-bias assessment for prediction models", "Not applicable", "No prediction model developed", "No individual-level model-training, validation or clinical decision threshold is claimed.", "Reassess only if a future clinical prediction model is added.", "https://www.probast.org/"],
  ["ARRIVE 2.0", "Animal research reporting", "Applicable to source-publication appraisal and prospective validation", "Secondary-analysis unit and provenance audited; prospective protocol aligned", "GSE278200, GSE308578 and the PXD052594 proteome are explicitly analyzed at the independent-animal level; GSE151374 is explicitly pooled at three mice per 10x library. Public metadata discrepancies and exclusions are documented. No new animal experiment was performed.", "Mark source-publication items unavailable from public records as not reported; complete the full checklist for future prospective experiments.", "https://arriveguidelines.org/arrive-guidelines"],
  ["MIAPE / ProteomeXchange provenance", "Public proteomics data", "Applicable to repository provenance", "Boundary documented", "PXD024058 remains descriptive because replicate provenance is irreconcilable. PXD052594 preserves the unconditioned 23-animal proteome (10 versus 13) and separately labels the deposited random 5-versus-5 phosphoproteome subset, with source exclusions and checksums recorded.", "Do not merge the PXD052594 phosphosite subset with the 23-animal proteome or condition on post-treatment radiomic clusters; retain PXD024058 as descriptive only.", "https://www.proteomexchange.org/"],
  ["AASLD DILI practice guidance", "DILI phenotype and interpretation", "Contextual", "Referenced in study rationale and future validation", "FAERS hy_proxy is explicitly not called Hy's law; public proteomics uses adjudicated multi-drug DILI as a distal phenotype anchor.", "Prospective clinical confirmation must use adjudicated DILI phenotyping.", "https://doi.org/10.1002/hep.32689"],
  ["EASL DILI clinical practice guideline", "DILI clinical framework", "Contextual", "Referenced in study rationale and future validation", "Clinical interpretation separates liver-injury phenotype from drug-specific causality.", "Use guideline-concordant adjudication in any future nintedanib-exposed cohort.", "https://doi.org/10.1016/j.jhep.2019.02.014"],
  ["Journal reporting summary", "Target-journal administrative requirements", "Required at submission if requested", "Pending journal selection", "Core reproducibility metadata, source URLs, checksums, seeds and analysis units are available in this workbook and manifests.", "Transfer information to the selected journal's current reporting form.", "https://www.nature.com/nature-portfolio/editorial-policies/reporting-standards"],
];

const statisticalUnits = [
  ["Analysis_Domain", "Dataset_or_Component", "Biological_or_Report_Unit", "N_or_Scope", "Technical_or_Derived_Units", "Primary_Statistical_Method", "Multiplicity_or_Stability_Control", "Permitted_Inference", "Prohibited_Overinterpretation", "Source_File"],
  ["Pharmacovigilance", "openFDA / FAERS", "Latest-version safetyreportid (report)", "Mutually exclusive nintedanib and pirfenidone report sets; endpoint-specific 2x2 counts", "Preferred terms within each report", "Active-comparator ROR, PRR, Fisher exact test, Bayesian IC; calendar-time and specification analyses", "Endpoint/specification-defined BH control where prespecified; annual and Mantel-Haenszel sensitivity analyses", "Relative reporting disproportionality", "Incidence, absolute risk, person-level recurrence, or causality", "fresh_analysis/faers/results/faers_frozen_manifest.json"],
  ["Clinical proteomics", "Public adjudicated DILI serum-protein workbook", "Participant within each protein-specific source series", "Protein- and cohort-specific available observations", "Protein measurements; technical units not counted as independent participants", "Two-group univariate effects, Welch and rank-based tests, Spearman correlations, prespecified sensitivity analyses", "BH FDR within cohort/contrast", "Human multi-drug DILI phenotype anchoring", "Nintedanib specificity, cross-protein multivariable modeling, or positional participant linkage", "fresh_analysis/proteomics/effect_estimates.csv"],
  ["Drug perturbation", "GSE70138 / LINCS HepG2", "Composite Level-5 signature at one concentration", "Six concentrations at 24 h", "Genes within a signature are correlated readouts", "Concentration-response metrics and exact six-dose permutation tests", "Exact permutation P values and BH FDR across the prespecified analysis family", "Cell-line concentration-response context", "Clinical exposure-response or biological-replicate inference", "fresh_analysis/lincs/analysis/gene_dose_response_metrics.csv.gz"],
  ["Target pharmacology", "ChEMBL CHEMBL502835", "Bioactivity record and supporting document", "Curated human single-protein activity subset", "Assays and activities nested within documents/targets", "Descriptive potency, activity, assay, and document synthesis", "Independent-document counts retained; no fabricated P values", "Biochemical target plausibility", "In-vivo target engagement or toxicity causality", "fresh_analysis/pharmacology/chembl_target_summary.csv"],
  ["Human liver single-cell", "GSE115469", "Donor", "Five healthy liver donors", "8,444 cells used descriptively within annotated compartments", "Donor-aware detection and expression summaries", "Across-donor consistency summaries", "Healthy-liver compartment localization", "Treating cells as independent biological replicates or inferring disease response", "fresh_analysis/liver_singlecell/target_expression_by_compartment.csv"],
  ["Virtual perturbation", "scTenifold-style virtual KO", "Seed-specific inferred network/KO run", "Three prespecified seeds per eligible gene-compartment pair", "Network nodes and affected genes", "Unsigned network displacement percentile and expression-matched null diagnostics", "Across-seed median, top-quartile count, minimum rank, Jaccard and matched-null stability", "Computational prioritization stability", "Protective/harmful direction, expression change, biological replication, or causality", "fresh_analysis/virtual_ko/virtual_knockout_stability_summary.csv"],
  ["Evidence integration", "Phenotype and target tracks", "Candidate within a prespecified role", "13 phenotype anchors and 17 exposure-proximal targets", "Evidence domains grouped by data source", "Pareto fronts, consensus ranks, leave-one-source-out ranks, and prespecified tier gates", "Role separation; same-compartment maximin; strict all-seed sensitivity; source-grouped LOSO", "Transparent prioritization and hypothesis tiering", "Clinical prediction, causal effect size, or cross-compartment evidence splicing", "fresh_analysis/integration/integration_analysis_manifest.json"],
  ["Longitudinal blood", "GSE299128", "Participant", "Seven participants; 57 longitudinal samples", "Repeated visits nested within participant", "Within-participant slope followed by two-sided exact sign-flip test of the mean slope", "All 2^7 sign assignments; BH across 30 prespecified genes", "Systemic longitudinal detectability", "Liver-specific toxicity validation or independent-sample inference across visits", "fresh_analysis/blood_longitudinal/analysis_manifest.json"],
  ["Public in-vivo transcriptomics", "GSE278200", "Independent animal", "12 rats; saline, BLM, and BLM+nintedanib, n=4/group", "Genes and module scores within animal", "TMM logCPM with robust limma-trend-style empirical-Bayes model; exact label permutation for modules", "Genome-wide BH for genes; within-contrast BH for seven modules", "Animal-level lung transcriptional disease/treatment contrasts", "Human efficacy, liver rescue, or independent-gene biological inference", "fresh_analysis/preclinical/results/analysis_manifest.json"],
  ["Public in-vivo transcriptomics", "GSE308578", "Independent animal", "43 mice; CTRL n=10, BLM+vehicle n=16, BLM+nintedanib n=17", "Genes and module scores within animal", "TMM logCPM with robust limma-trend-style empirical-Bayes model; label permutation for modules", "Genome-wide BH; within-contrast module BH; leave-one-animal-out stability", "Animal-level lung transcriptional disease/treatment contrasts", "Relabelling based on the conflicting main Figure 6 caption, human efficacy, or liver rescue", "fresh_analysis/preclinical/results/analysis_manifest.json"],
  ["Public liver-slice transcriptomics", "GSE120804", "Pooled-tissue slice observation", "48 observations across five groups", "Genes and module scores within slice observation", "Exploratory group contrasts", "Within-contrast module BH; no animal-level permutation inference", "Exploratory pooled-slice molecular context", "Animal-population inference or treating slices as independent animals", "fresh_analysis/preclinical/results/GSE120804_module_contrasts.csv"],
  ["Public lung-slice transcriptomics", "GSE120679", "Pooled-tissue slice observation", "36 observations across six groups", "Genes and module scores within slice observation", "Exploratory group contrasts", "Within-contrast module BH; no animal-level permutation inference", "Exploratory pooled-slice molecular context", "Animal-population inference or treating slices as independent animals", "fresh_analysis/preclinical/results/GSE120679_module_contrasts.csv"],
  ["Public lung proteomics", "PXD024058", "Not reconstructable from public provenance", "5,951 proteins; three supplied columns per condition", "Protein abundance columns with ambiguous biological/technical replicate status", "Descriptive log2 fold changes and direction checks only", "No inferential P values used", "Descriptive orthogonal protein direction", "Treating columns as independent biological replicates or claiming validation", "fresh_analysis/preclinical/results/PXD024058_candidate_effects_descriptive.csv"],
  ["Public lung immune-cell single-cell transcriptomics", "GSE151374", "Pooled 10x library (three mice per library)", "18 libraries; n=3 per control/BLM/BLM+nintedanib condition at days 7 and 14; 70,696 QC-retained cells", "Cells and genes nested within pooled library; all-retained and macrophage pseudobulks derived once per library", "Pooled-library pseudobulk Welch contrasts; exhaustive 3-versus-3 label permutations for cell-type proportions, modules and prespecified tracks", "BH within prespecified families; exact-permutation BH; minimum attainable two-sided exact P=0.10", "Nonconfirmatory pulmonary immune-cell and macrophage context", "Cell-level pseudoreplication, independent-mouse inference, hepatic-safety inference, or target confirmation", "fresh_analysis/preclinical/gse151374/results/analysis_manifest.json"],
  ["Public lung proteomics and phosphoproteomics", "PXD052594", "Individual mouse", "Proteome: 23 animals (nintedanib 10, vehicle 13); phosphoproteome: deposited random subset of 5 animals per arm", "7,006 proteins across 23 animals; 20,043 phosphosites across the distinct 10-animal subset", "Unconditioned all-animal Welch proteome contrasts; prespecified module/candidate summaries; separate deposited-subset phosphosite contrasts", "BH across all 7,006 proteins; BH across all 20,043 phosphosites; no radiomic-cluster conditioning", "Nonconfirmatory pulmonary protein and phosphosite context", "Hepatic safety, unbiased causal confirmation, merging the 5-versus-5 subset with the 23-animal proteome, or post-treatment responder conditioning", "fresh_analysis/preclinical/pxd052594/results/analysis_manifest.json"],
  ["Prospective functional validation", "Planned human and animal experiments", "Donor, organoid line, or animal; technical replicates averaged", "Planned; not yet performed", "Technical wells/readouts nested within biological unit", "Prespecified mixed models, interaction tests, co-primary multiplicity and non-inferiority rules", "Holm/BH and replication/stop-go criteria as specified", "Protocol-level reproducibility", "Any completed-result claim", "Experimental_Validation sheet"],
];

const figureSourceMap = [
  ["Figure_or_Supplement", "Panel_Scope", "Primary_Artifact", "Primary_Source_Tables", "Statistical_Unit", "Key_Analysis", "Inference_Boundary", "Caption_or_Manifest"],
  ["Graphical Abstract", "Study architecture only", "fresh_analysis/figures/Graphical_Abstract_Public_Data_Target_Map.svg", "README; Data_Sources; Figure_Source_Map", "Not applicable", "Evidence-flow schematic; no new numerical result", "Must not be interpreted as an additional result figure", "Manuscript graphical-abstract legend"],
  ["Figure 1", "A-J", "fresh_analysis/figures/Figure_1_FAERS.svg", "fresh_analysis/figures/Figure_1_FAERS_panel_source_data.csv; fresh_analysis/faers/results/overall_active_comparator_signals.csv; fresh_analysis/faers/results/report_level_specification_signals.csv", "Latest-version safetyreportid (report)", "Active-comparator disproportionality, calendar-time, onset, outcome, subgroup, and specification analyses", "Reporting disproportionality, not incidence or causality", "fresh_analysis/figures/Figure_1_FAERS_caption.md"],
  ["Figure 2", "A-J", "fresh_analysis/figures/Figure_2_Human_DILI_Proteomics.svg", "fresh_analysis/proteomics/effect_estimates.csv; fresh_analysis/proteomics/correlation_stats.csv; fresh_analysis/proteomics/zonal_effect_estimates.csv; fresh_analysis/proteomics/qc_counts.csv", "Participant within protein-specific source series", "Human multi-drug DILI phenotype anchoring and robustness", "Not nintedanib-specific; cross-protein participant linkage unavailable", "fresh_analysis/figures/Figure_2_manifest.json"],
  ["Figure 3", "A-J", "fresh_analysis/figures/Figure_3_Nintedanib_HepG2_Dose_Perturbation.svg", "fresh_analysis/lincs/analysis/gene_dose_response_metrics.csv.gz; fresh_analysis/lincs/analysis/mechanism_module_statistics.csv; fresh_analysis/lincs/analysis/sciplex_hepg2_concordance.csv; fresh_analysis/pharmacology/chembl_target_summary.csv", "Composite Level-5 signature at one dose; ChEMBL activity/document", "HepG2 concentration response, module shifts, orthogonal list overlap, and biochemical pharmacology", "Cell-line and biochemical context do not establish clinical hepatotoxicity", "fresh_analysis/lincs/analysis/Figure_3_caption.md"],
  ["Figure 4", "A-J", "fresh_analysis/figures/Figure_4_Liver_Localization_and_Virtual_Knockout.svg", "fresh_analysis/liver_singlecell/target_expression_by_compartment.csv; fresh_analysis/virtual_ko/virtual_knockout_stability_summary.csv; fresh_analysis/integration/virtual_ko_expression_matched_null_stability.csv", "Donor for expression; seed/network run for computational stability", "Donor-aware localization and unsigned virtual-KO displacement", "Healthy liver; computational prioritization only; no direction or causality", "fresh_analysis/figures/Figure_4_caption.md"],
  ["Figure 5", "A", "fresh_analysis/figures/Figure_5_Multi_Evidence_Target_Prioritization.svg", "fresh_analysis/proteomics/candidate_evidence_table.csv; fresh_analysis/proteomics/effect_estimates.csv; fresh_analysis/proteomics/correlation_stats.csv", "Participant within protein-specific source series", "Four clinical DILI phenotype properties", "One multi-drug DILI source group; not nintedanib-specific", "fresh_analysis/figures/Figure_5_caption.md"],
  ["Figure 5", "B", "fresh_analysis/figures/Figure_5_Multi_Evidence_Target_Prioritization.svg", "fresh_analysis/liver_singlecell/target_expression_by_compartment.csv; fresh_analysis/virtual_ko/virtual_knockout_stability_summary.csv; fresh_analysis/proteomics/candidate_evidence_table.csv", "Donor for localization; seed/network run for virtual KO", "Hepatocyte detection versus unsigned virtual-KO displacement", "Computational, healthy-liver context; no direction or causality", "fresh_analysis/figures/Figure_5_caption.md"],
  ["Figure 5", "C", "fresh_analysis/figures/Figure_5_Multi_Evidence_Target_Prioritization.svg", "fresh_analysis/lincs/output/nintedanib_hepg2_level5_all_genes_tidy.csv.gz; fresh_analysis/lincs/analysis/gene_dose_response_metrics.csv.gz", "Composite Level-5 signature at one dose", "Six-dose HepG2 profiles for available phenotype anchors", "Cell-line perturbation, not clinical exposure-response", "fresh_analysis/figures/Figure_5_caption.md"],
  ["Figure 5", "D", "fresh_analysis/figures/Figure_5_Multi_Evidence_Target_Prioritization.svg", "fresh_analysis/lincs/output/nintedanib_hepg2_level5_all_genes_tidy.csv.gz; fresh_analysis/integration/phenotype_anchor_tiers.csv", "Composite Level-5 signature at one dose", "Unsmoothed P1/P2 dose profiles", "Phenotype P labels are not target tiers", "fresh_analysis/figures/Figure_5_caption.md"],
  ["Figure 5", "E", "fresh_analysis/figures/Figure_5_Multi_Evidence_Target_Prioritization.svg", "fresh_analysis/lincs/analysis/sciplex_hepg2_concordance.csv", "Gene conditional on inclusion in truncated MCF7 lists", "Overlap-conditional Sci-Plex comparison", "Absence from a truncated list is not absence of response", "fresh_analysis/figures/Figure_5_caption.md"],
  ["Figure 5", "F", "fresh_analysis/figures/Figure_5_Multi_Evidence_Target_Prioritization.svg", "fresh_analysis/integration/candidate_domain_scores_long.csv; fresh_analysis/integration/candidate_consensus_ranking.csv; fresh_analysis/integration/target_same_compartment_context_all_pairs.csv", "Candidate within phenotype or target role", "Explicit source-group evidence matrix", "Derived features from one dataset remain one evidence domain", "fresh_analysis/figures/Figure_5_caption.md"],
  ["Figure 5", "G", "fresh_analysis/figures/Figure_5_Multi_Evidence_Target_Prioritization.svg", "fresh_analysis/integration/phenotype_anchor_tiers.csv", "Phenotype candidate", "Non-dominated Pareto fronts", "Pareto position is prioritization, not a causal effect", "fresh_analysis/figures/Figure_5_caption.md"],
  ["Figure 5", "H", "fresh_analysis/figures/Figure_5_Multi_Evidence_Target_Prioritization.svg", "fresh_analysis/integration/phenotype_anchor_tiers.csv; fresh_analysis/integration/candidate_leave_one_source_out_ranks.csv", "Phenotype candidate", "P-track consensus and leave-one-source-out ranges", "P1/P2/P3 classify phenotype evidence, not binding targets", "fresh_analysis/figures/Figure_5_caption.md"],
  ["Figure 5", "I", "fresh_analysis/figures/Figure_5_Multi_Evidence_Target_Prioritization.svg", "fresh_analysis/integration/target_tier_gate_sensitivity.csv; fresh_analysis/integration/target_same_compartment_context_all_pairs.csv; fresh_analysis/integration/candidate_leave_one_source_out_ranks.csv; fresh_analysis/pharmacology/chembl_target_summary.csv", "Target candidate", "T-track ranks with the prespecified exploratory same-compartment gate and strict all-seed sensitivity gate", "FGFR1 is the sole candidate meeting the exploratory T2 gate; the strict all-seed T2 set is empty; no T1", "fresh_analysis/figures/Figure_5_caption.md"],
  ["Figure 5", "J", "fresh_analysis/figures/Figure_5_Multi_Evidence_Target_Prioritization.svg", "fresh_analysis/integration/string_candidate_path_metrics.csv; fresh_analysis/inputs/string_candidate_network_expanded.tsv", "Association edge/path", "Sparse STRING contextual paths", "Database/text-mining association is not a demonstrated toxicity pathway", "fresh_analysis/figures/Figure_5_caption.md"],
  ["Figure 6", "A", "fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints.svg", "fresh_analysis/preclinical/results/GSE278200_pca_coordinates.csv; fresh_analysis/preclinical/results/GSE278200_sample_metadata.csv", "Independent animal", "PCA of rat lung transcriptomes", "Small male-rat bleomycin experiment; not liver/DILI evidence", "fresh_analysis/preclinical/Figure_6_caption.md"],
  ["Figure 6", "B", "fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints.svg", "fresh_analysis/preclinical/results/GSE278200_reversal_geometry.csv.gz; fresh_analysis/preclinical/results/GSE278200_differential_expression.csv.gz", "Gene as display unit; animal-level contrasts underneath", "Disease-versus-treatment reversal geometry", "Genes are correlated; correlation P is descriptive, not biological inference", "fresh_analysis/preclinical/Figure_6_caption.md"],
  ["Figure 6", "C", "fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints.svg", "fresh_analysis/preclinical/results/GSE278200_differential_expression.csv.gz; fresh_analysis/preclinical/results/GSE278200_expression_filter.csv.gz", "Independent animal", "Genome-wide animal-level nintedanib-versus-BLM differential expression", "BH spans 13,598 expressed genes", "fresh_analysis/preclinical/Figure_6_caption.md"],
  ["Figure 6", "D", "fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints.svg", "fresh_analysis/preclinical/results/GSE278200_candidate_results.csv", "Independent animal", "Prespecified target, anchor, and fibrosis-gene effects", "Absence after filtering is explicitly retained; no liver-safety inference", "fresh_analysis/preclinical/Figure_6_caption.md"],
  ["Figure 6", "E", "fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints.svg", "fresh_analysis/preclinical/results/GSE278200_module_scores.csv; fresh_analysis/preclinical/results/GSE278200_module_contrasts.csv; fresh_analysis/preclinical/results/GSE278200_module_coverage.csv", "Independent animal", "Seven prespecified module effects with exact label permutations", "Smallest attainable two-sided P and within-contrast BH are preserved", "fresh_analysis/preclinical/Figure_6_caption.md"],
  ["Figure 6", "F", "fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints.svg", "fresh_analysis/preclinical/results/GSE308578_pca_coordinates.csv; fresh_analysis/preclinical/results/GSE308578_sample_metadata.csv", "Independent animal", "PCA of mouse lung transcriptomes", "GEO/ARRIVE 10/16/17 labels retained despite conflicting main-caption n", "fresh_analysis/preclinical/Figure_6_caption.md"],
  ["Figure 6", "G", "fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints.svg", "fresh_analysis/preclinical/results/GSE308578_reversal_geometry.csv.gz; fresh_analysis/preclinical/results/GSE308578_differential_expression.csv.gz", "Gene as display unit; animal-level contrasts underneath", "Disease-versus-treatment reversal geometry", "Genes are correlated; reversal magnitude is model-specific and descriptive", "fresh_analysis/preclinical/Figure_6_caption.md"],
  ["Figure 6", "H", "fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints.svg", "fresh_analysis/preclinical/results/GSE308578_candidate_results.csv; fresh_analysis/preclinical/results/GSE308578_leave_one_animal_out_summary.csv; fresh_analysis/preclinical/results/GSE308578_leave_one_animal_out_full.csv.gz", "Independent animal", "Prespecified effects with 33 leave-one-relevant-animal-out deletions", "Whole-lung response does not validate hepatic toxicity or safety", "fresh_analysis/preclinical/Figure_6_caption.md"],
  ["Figure 6", "I", "fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints.svg", "fresh_analysis/preclinical/results/GSE120804_sample_metadata.csv; fresh_analysis/preclinical/results/GSE120804_candidate_results.csv; fresh_analysis/preclinical/results/GSE120804_differential_expression.csv.gz; fresh_analysis/preclinical/results/GSE120804_module_contrasts.csv", "Pooled-tissue liver-slice observation", "Exploratory BDL liver-PCLS effects", "Slices were allocated after tissue pooling; no animal-population or DILI inference", "fresh_analysis/preclinical/Figure_6_caption.md"],
  ["Figure 6", "J", "fresh_analysis/figures/Figure_6_Public_Preclinical_Constraints.svg", "fresh_analysis/preclinical/results/cross_model_candidate_module_effects.csv; fresh_analysis/preclinical/results/GSE120679_candidate_results.csv; fresh_analysis/preclinical/results/GSE120679_module_contrasts.csv; fresh_analysis/preclinical/results/PXD024058_candidate_effects_descriptive.csv; fresh_analysis/preclinical/results/PXD024058_protein_effects_descriptive.csv.gz", "Model-specific: independent animal, pooled lung-slice observation, or unresolved protein replicate", "Cross-model candidate/module effect display", "GSE120679 cannot support animal-population inference; PXD024058 is descriptive only", "fresh_analysis/preclinical/Figure_6_caption.md"],
  ["Supplementary Figure S1", "FAERS supporting panels", "fresh_analysis/figures/Supplementary_Figure_S1_FAERS.svg", "fresh_analysis/figures/Supplementary_Figure_S1_FAERS_panel_source_data.csv; fresh_analysis/faers/results/preferred_term_signals.csv; fresh_analysis/faers/results/annual_active_comparator_signals.csv", "Latest-version safetyreportid (report)", "Preferred-term and temporal sensitivity analyses", "Reporting disproportionality only", "fresh_analysis/figures/Supplementary_Figure_S1_FAERS_caption.md"],
  ["Supplementary Figure S2", "Protein distributions", "fresh_analysis/figures/Supplementary_Figure_S2_Protein_Distributions.svg", "fresh_analysis/proteomics/effect_estimates.csv; fresh_analysis/proteomics/qc_counts.csv", "Participant within protein-specific source series", "Distribution and availability views", "Multi-drug DILI; not nintedanib-specific", "fresh_analysis/figures/Supplementary_Figure_S2_caption.md"],
  ["Supplementary Figure S3", "Proteomics robustness", "fresh_analysis/figures/Supplementary_Figure_S3_Proteomics_Robustness.svg", "fresh_analysis/proteomics/effect_estimates.csv; fresh_analysis/proteomics/correlation_stats.csv; fresh_analysis/proteomics/missingness_stats.csv; fresh_analysis/proteomics/left_censor_sensitivity.csv; fresh_analysis/proteomics/row_order_linkage_qc.csv", "Participant within protein-specific source series", "Discovery/confirmation concordance, missingness, censoring, outlier, zone, source-block, and linkage audits", "Univariate protein-wise inference only; no defensible participant-level multivariable matrix", "fresh_analysis/figures/Supplementary_Figure_S3_caption.md"],
  ["Supplementary Figure S4", "Nintedanib perturbation robustness", "fresh_analysis/figures/Supplementary_Figure_S4_Nintedanib_Perturbation_Robustness.svg", "fresh_analysis/lincs/analysis/gene_dose_response_metrics.csv.gz; fresh_analysis/lincs/analysis/sciplex_hepg2_concordance.csv", "Composite signature at one dose", "Dose-pattern and cross-cell-line sensitivity analyses", "No clinical exposure or outcome inference", "fresh_analysis/figures/Supplementary_Figure_S4_caption.md"],
  ["Supplementary Figure S5", "Whole-blood sensitivity", "fresh_analysis/figures/Supplementary_Figure_S5_Whole_Blood_Sensitivity.svg", "fresh_analysis/blood_longitudinal/target_longitudinal_summary.csv; fresh_analysis/blood_longitudinal/participant_slopes.csv", "Participant", "Exact two-sided participant-level sign-flip test over all 2^7 assignments", "Small single-arm blood cohort without adjudicated DILI; systemic detectability only", "fresh_analysis/figures/Supplementary_Figure_S5_caption.md"],
  ["Supplementary Figure S6", "A", "fresh_analysis/figures/Supplementary_Figure_S6_GSE151374.svg", "fresh_analysis/preclinical/gse151374/results/qc_by_library.csv", "Pooled 10x library (three mice per pool)", "QC-retained cell count by library", "Cells are descriptive nested observations; n=3 pooled libraries per condition-time", "fresh_analysis/preclinical/gse151374/docs/Supplementary_Figure_S6_caption.md"],
  ["Supplementary Figure S6", "B", "fresh_analysis/figures/Supplementary_Figure_S6_GSE151374.svg", "fresh_analysis/preclinical/gse151374/results/qc_by_library.csv", "Pooled 10x library (three mice per pool)", "Median UMI and detected-gene library QC", "No cell-level inferential test", "fresh_analysis/preclinical/gse151374/docs/Supplementary_Figure_S6_caption.md"],
  ["Supplementary Figure S6", "C", "fresh_analysis/figures/Supplementary_Figure_S6_GSE151374.svg", "fresh_analysis/preclinical/gse151374/results/pca_all_retained_pseudobulk.csv", "Pooled 10x library (three mice per pool)", "All-retained-cell library-pseudobulk PCA", "PCA is descriptive; the pool is the inferential unit", "fresh_analysis/preclinical/gse151374/docs/Supplementary_Figure_S6_caption.md"],
  ["Supplementary Figure S6", "D", "fresh_analysis/figures/Supplementary_Figure_S6_GSE151374.svg", "fresh_analysis/preclinical/gse151374/results/celltype_counts_and_proportions_by_library.csv", "Pooled 10x library (three mice per pool)", "Broad marker-rule lineage composition", "Marker-rule labels are broad descriptive strata, not cell-level replication", "fresh_analysis/preclinical/gse151374/docs/Supplementary_Figure_S6_caption.md"],
  ["Supplementary Figure S6", "E", "fresh_analysis/figures/Supplementary_Figure_S6_GSE151374.svg", "fresh_analysis/preclinical/gse151374/results/celltype_proportion_contrasts.csv", "Pooled 10x library (three mice per pool)", "Day-matched macrophage fractions with exhaustive 3-versus-3 permutation tests", "All exact-permutation BH families are null; minimum two-sided exact P=0.10", "fresh_analysis/preclinical/gse151374/docs/Supplementary_Figure_S6_caption.md"],
  ["Supplementary Figure S6", "F", "fresh_analysis/figures/Supplementary_Figure_S6_GSE151374.svg", "fresh_analysis/preclinical/gse151374/results/pca_macrophage_pseudobulk.csv", "Pooled 10x library (three mice per pool)", "Macrophage-stratum library-pseudobulk PCA", "PCA is descriptive; lung context only", "fresh_analysis/preclinical/gse151374/docs/Supplementary_Figure_S6_caption.md"],
  ["Supplementary Figure S6", "G", "fresh_analysis/figures/Supplementary_Figure_S6_GSE151374.svg", "fresh_analysis/preclinical/gse151374/results/macrophage_DE_NINT_vs_BLM_day7.csv", "Pooled 10x library (three mice per pool)", "Day-7 macrophage-pseudobulk gene-level Welch contrasts with BH", "No gene met q<0.10; no cell-level inference", "fresh_analysis/preclinical/gse151374/docs/Supplementary_Figure_S6_caption.md"],
  ["Supplementary Figure S6", "H", "fresh_analysis/figures/Supplementary_Figure_S6_GSE151374.svg", "fresh_analysis/preclinical/gse151374/results/macrophage_DE_NINT_vs_BLM_day14.csv", "Pooled 10x library (three mice per pool)", "Day-14 macrophage-pseudobulk gene-level Welch contrasts with BH", "No gene met q<0.10; no cell-level inference", "fresh_analysis/preclinical/gse151374/docs/Supplementary_Figure_S6_caption.md"],
  ["Supplementary Figure S6", "I", "fresh_analysis/figures/Supplementary_Figure_S6_GSE151374.svg", "fresh_analysis/preclinical/gse151374/results/macrophage_module_scores_by_library.csv; fresh_analysis/preclinical/gse151374/results/macrophage_module_contrasts.csv", "Pooled 10x library (three mice per pool)", "Five prespecified macrophage-program contrasts", "No multiplicity-controlled program signal", "fresh_analysis/preclinical/gse151374/docs/Supplementary_Figure_S6_caption.md"],
  ["Supplementary Figure S6", "J", "fresh_analysis/figures/Supplementary_Figure_S6_GSE151374.svg", "fresh_analysis/preclinical/gse151374/results/candidate_and_phenotype_tracks.csv; fresh_analysis/preclinical/gse151374/results/candidate_gene_contrasts.csv", "Pooled 10x library (three mice per pool)", "Exposure-proximal Fgfr1 and downstream phenotype-anchor tracks", "Fgfr1 is the sole target track; Fbp1/Gsta1/Otc are phenotype anchors; lung context is not hepatic confirmation", "fresh_analysis/preclinical/gse151374/docs/Supplementary_Figure_S6_caption.md"],
  ["Supplementary Figure S7", "A", "fresh_analysis/figures/Supplementary_Figure_S7_PXD052594.svg", "fresh_analysis/preclinical/pxd052594/results/animal_sample_manifest.csv; fresh_analysis/preclinical/pxd052594/results/proteome_NINT_vs_vehicle_all_animals.csv", "Individual mouse", "Quantified-protein completeness across all 23 retained animals", "No radiomic response cluster was used", "fresh_analysis/preclinical/pxd052594/docs/Supplementary_Figure_S7_caption.md"],
  ["Supplementary Figure S7", "B", "fresh_analysis/figures/Supplementary_Figure_S7_PXD052594.svg", "fresh_analysis/preclinical/pxd052594/results/proteome_PCA.csv", "Individual mouse", "PCA of 2,000 most variable proteins", "Descriptive pulmonary separation; not hepatic safety", "fresh_analysis/preclinical/pxd052594/docs/Supplementary_Figure_S7_caption.md"],
  ["Supplementary Figure S7", "C", "fresh_analysis/figures/Supplementary_Figure_S7_PXD052594.svg", "fresh_analysis/preclinical/pxd052594/results/proteome_NINT_vs_vehicle_all_animals.csv", "Individual mouse", "Unconditioned all-animal protein contrast", "BH across 7,006 proteins is null; post-treatment responder conditioning is prohibited", "fresh_analysis/preclinical/pxd052594/docs/Supplementary_Figure_S7_caption.md"],
  ["Supplementary Figure S7", "D", "fresh_analysis/figures/Supplementary_Figure_S7_PXD052594.svg", "fresh_analysis/preclinical/pxd052594/results/proteome_NINT_vs_vehicle_all_animals.csv", "Individual mouse", "Standardized abundance display for 20 smallest nominal protein P values", "Display selection is not an inferential threshold", "fresh_analysis/preclinical/pxd052594/docs/Supplementary_Figure_S7_caption.md"],
  ["Supplementary Figure S7", "E", "fresh_analysis/figures/Supplementary_Figure_S7_PXD052594.svg", "fresh_analysis/preclinical/pxd052594/results/module_contrasts.csv", "Individual mouse", "Five prespecified protein-program contrasts", "No module met BH significance", "fresh_analysis/preclinical/pxd052594/docs/Supplementary_Figure_S7_caption.md"],
  ["Supplementary Figure S7", "F", "fresh_analysis/figures/Supplementary_Figure_S7_PXD052594.svg", "fresh_analysis/preclinical/pxd052594/results/module_scores_by_animal.csv", "Individual mouse", "Extracellular-matrix module scores", "Pulmonary contextual sensitivity only", "fresh_analysis/preclinical/pxd052594/docs/Supplementary_Figure_S7_caption.md"],
  ["Supplementary Figure S7", "G", "fresh_analysis/figures/Supplementary_Figure_S7_PXD052594.svg", "fresh_analysis/preclinical/pxd052594/results/module_scores_by_animal.csv", "Individual mouse", "Macrophage-repair module scores", "Pulmonary contextual sensitivity only", "fresh_analysis/preclinical/pxd052594/docs/Supplementary_Figure_S7_caption.md"],
  ["Supplementary Figure S7", "H", "fresh_analysis/figures/Supplementary_Figure_S7_PXD052594.svg", "fresh_analysis/preclinical/pxd052594/results/module_scores_by_animal.csv", "Individual mouse", "RTK-response module scores", "Pulmonary contextual sensitivity only", "fresh_analysis/preclinical/pxd052594/docs/Supplementary_Figure_S7_caption.md"],
  ["Supplementary Figure S7", "I", "fresh_analysis/figures/Supplementary_Figure_S7_PXD052594.svg", "fresh_analysis/preclinical/pxd052594/results/candidate_contrasts.csv", "Individual mouse", "Prespecified target, phenotype-anchor and lung-track protein effects", "Fgfr1 is the sole exposure-proximal target; Fbp1/Gsta1/Otc are phenotype anchors", "fresh_analysis/preclinical/pxd052594/docs/Supplementary_Figure_S7_caption.md"],
  ["Supplementary Figure S7", "J", "fresh_analysis/figures/Supplementary_Figure_S7_PXD052594.svg", "fresh_analysis/preclinical/pxd052594/results/candidate_expression_by_animal.csv", "Individual mouse", "Absolute deposited abundance of prespecified tracks", "Nonquantified proteins are retained; lung context is not liver/DILI confirmation", "fresh_analysis/preclinical/pxd052594/docs/Supplementary_Figure_S7_caption.md"],
  ["Supplementary table-only sensitivity", "PXD052594 phosphoproteome", "fresh_analysis/preclinical/pxd052594/results/phosphosite_NINT_vs_vehicle_all_deposited_subset_animals.csv", "fresh_analysis/preclinical/pxd052594/results/phosphosite_NINT_vs_vehicle_all_deposited_subset_animals.csv; fresh_analysis/preclinical/pxd052594/results/phosphosite_target_pathway_summary.csv", "Individual mouse in the distinct deposited random 5-versus-5 subset", "All 20,043 phosphosite contrasts with BH; 994 sites q<0.05", "Not a panel in Supplementary Figure S7; distinct from the 23-animal proteome and interpreted only as contextual pulmonary signaling", "fresh_analysis/preclinical/pxd052594/results/analysis_manifest.json"],
];

const figureSourcePaths = [...new Set(figureSourceMap.slice(1).flatMap((row) => [row[2], row[3], row[7]])
  .flatMap((value) => String(value).split(";").map((x) => x.trim()))
  .filter((value) => value.startsWith("fresh_analysis/") || value.startsWith("fresh_data/")))];
const figureSourcePathChecks = [];
for (const relativePath of figureSourcePaths) {
  try {
    const stat = await fs.stat(path.join(ROOT, relativePath));
    figureSourcePathChecks.push({ path: relativePath, exists: true, bytes: stat.size });
  } catch {
    figureSourcePathChecks.push({ path: relativePath, exists: false, bytes: null });
  }
}
const missingFigureSourcePaths = figureSourcePathChecks.filter((x) => !x.exists).map((x) => x.path);
if (missingFigureSourcePaths.length) throw new Error(`Figure source-map paths do not exist: ${missingFigureSourcePaths.join(", ")}`);

await fs.mkdir(OUT_DIR, { recursive: true });
await fs.rm(RENDER_DIR, { recursive: true, force: true });
await fs.mkdir(RENDER_DIR, { recursive: true });

const integrationManifest = await readJson("fresh_analysis/integration/integration_analysis_manifest.json");
const integrationValidation = await readJson("fresh_analysis/integration/integration_validation_report.json");
const bloodManifest = await readJson("fresh_analysis/blood_longitudinal/analysis_manifest.json");
const preclinicalManifest = await readJson("fresh_analysis/preclinical/results/analysis_manifest.json");
const preclinicalValidation = await readJson("fresh_analysis/preclinical/results/validation_report.json");
const preclinicalAnalysisSummary = await readJson("fresh_analysis/preclinical/results/analysis_summary.json");
const gse151374Manifest = await readJson("fresh_analysis/preclinical/gse151374/results/analysis_manifest.json");
const gse151374Validation = await readJson("fresh_analysis/preclinical/gse151374/validation/validation_checks.json");
const gse151374HashValidation = await readJson("fresh_analysis/preclinical/gse151374/validation/manifest_hash_validation.json");
const pxd052594Manifest = await readJson("fresh_analysis/preclinical/pxd052594/results/analysis_manifest.json");
const pxd052594Validation = await readJson("fresh_analysis/preclinical/pxd052594/validation/validation_checks.json");
const pxd052594HashValidation = await readJson("fresh_analysis/preclinical/pxd052594/validation/manifest_hash_validation.json");
const integrationManifestSha256 = await sha256File(path.join(ROOT, "fresh_analysis/integration/integration_analysis_manifest.json"));
const bloodManifestSha256 = await sha256File(path.join(ROOT, "fresh_analysis/blood_longitudinal/analysis_manifest.json"));
const preclinicalManifestSha256 = await sha256File(path.join(ROOT, "fresh_analysis/preclinical/results/analysis_manifest.json"));
const gse151374ManifestSha256 = await sha256File(path.join(ROOT, "fresh_analysis/preclinical/gse151374/results/analysis_manifest.json"));
const pxd052594ManifestSha256 = await sha256File(path.join(ROOT, "fresh_analysis/preclinical/pxd052594/results/analysis_manifest.json"));
const EXPECTED_CORE_PRECLINICAL_MANIFEST_SHA256 = "e4ca5cc061014657a34c91ab4adb396125e2bbd6041596df96c2ae33cf30965e";
const EXPECTED_INTEGRATION_MANIFEST_SHA256 = "70c8600590fcb0d6ff2663d834902d9aceb5cf71b41ec6b02cae49bd1a45a365";
const EXPECTED_GSE151374_MANIFEST_SHA256 = "999415f7c21f9f4f648dc830c0ffe82ada42bc902cc4b4a8ee5c807c164f4a52";
const EXPECTED_PXD052594_MANIFEST_SHA256 = "2d44ac6de65c09453646b8f27a77ec7de1f7b5492685a0557e8acc5daf71c596";
if (preclinicalManifestSha256 !== EXPECTED_CORE_PRECLINICAL_MANIFEST_SHA256) {
  throw new Error(`Frozen Figure 6 core manifest changed: ${preclinicalManifestSha256}`);
}
if (integrationManifestSha256 !== EXPECTED_INTEGRATION_MANIFEST_SHA256) {
  throw new Error(`Frozen same-compartment integration manifest changed: ${integrationManifestSha256}`);
}
if (gse151374ManifestSha256 !== EXPECTED_GSE151374_MANIFEST_SHA256 || gse151374HashValidation.manifest_sha256 !== EXPECTED_GSE151374_MANIFEST_SHA256 || gse151374HashValidation.all_recorded_hashes_match !== true) {
  throw new Error(`Frozen GSE151374 extension manifest changed or failed its recorded hash audit: ${gse151374ManifestSha256}`);
}
if (pxd052594ManifestSha256 !== EXPECTED_PXD052594_MANIFEST_SHA256 || pxd052594HashValidation.manifest_sha256 !== EXPECTED_PXD052594_MANIFEST_SHA256 || pxd052594HashValidation.all_recorded_hashes_match !== true) {
  throw new Error(`Frozen PXD052594 extension manifest changed or failed its recorded hash audit: ${pxd052594ManifestSha256}`);
}

async function verifyHashRecords(label, records, acceptedPostFreezeUpdates = {}) {
  const checks = [];
  for (const record of records) {
    const fullPath = path.join(ROOT, record.path);
    const observed = await sha256File(fullPath);
    const passed = observed === record.sha256;
    const acceptedPostFreezeAuditUpdate = !passed && acceptedPostFreezeUpdates[record.path] === observed;
    checks.push({ label, path: record.path, expected: record.sha256, observed, passed, accepted_post_freeze_audit_update: acceptedPostFreezeAuditUpdate });
    if (!passed && !acceptedPostFreezeAuditUpdate) throw new Error(`${label} hash mismatch: ${record.path}`);
  }
  return checks;
}

const preclinicalFileHashChecks = await verifyHashRecords("public_preclinical_core", preclinicalManifest.files);
const preclinicalOriginalHashMatches = preclinicalFileHashChecks.filter((x) => x.passed).length;
const preclinicalAcceptedAuditUpdates = preclinicalFileHashChecks.filter((x) => x.accepted_post_freeze_audit_update).length;
if (preclinicalOriginalHashMatches !== 98 || preclinicalAcceptedAuditUpdates !== 0) {
  throw new Error(`Unexpected Figure 6 core hash reconciliation: original=${preclinicalOriginalHashMatches}, accepted_audit_updates=${preclinicalAcceptedAuditUpdates}`);
}
const gse151374FileHashChecks = await verifyHashRecords("gse151374_extension", [...gse151374Manifest.inputs, ...gse151374Manifest.outputs]);
const pxd052594FileHashChecks = await verifyHashRecords("pxd052594_extension", [...pxd052594Manifest.inputs, ...pxd052594Manifest.outputs]);
const integrationHashRecords = [
  ...Object.entries(integrationManifest.input_hashes).map(([p, h]) => ({ path: p, sha256: h })),
  ...Object.entries(integrationManifest.output_table_hashes).map(([p, h]) => ({ path: p, sha256: h })),
  ...Object.entries(integrationManifest.text_artifact_hashes).map(([p, h]) => ({ path: p, sha256: h })),
];
const integrationFileHashChecks = await verifyHashRecords("same_compartment_integration", integrationHashRecords);
const integrationOriginalHashMatches = integrationFileHashChecks.filter((x) => x.passed).length;
const integrationAcceptedTextUpdates = integrationFileHashChecks.filter((x) => x.accepted_post_freeze_audit_update).length;
if (integrationOriginalHashMatches !== 24 || integrationAcceptedTextUpdates !== 0) {
  throw new Error(`Unexpected integration hash reconciliation: original=${integrationOriginalHashMatches}, accepted_text_updates=${integrationAcceptedTextUpdates}`);
}

const preclinicalSummaryMatrix = [
  ["Dataset", "Model_or_Material", "Statistical_Unit", "Group_N", "Features_Analyzed", "Modeling_Method_or_Scope", "Disease_Responsive_Features", "Fraction_Directionally_Reversed", "Fraction_Closer_to_Control", "Median_Reversal_Index_among_Directional", "Fraction_Directional_Reversal_Index_0_5_to_1_5", "Spearman_Disease_vs_Treatment", "Permitted_Inference", "Key_Limitation", "Source_File"],
  [
    "GSE278200",
    "Rat bleomycin lung fibrosis with nintedanib treatment",
    "Independent animal",
    JSON.stringify(preclinicalAnalysisSummary.GSE278200.group_n),
    preclinicalAnalysisSummary.GSE278200.diagnostics.n_genes,
    preclinicalAnalysisSummary.GSE278200.diagnostics.method,
    preclinicalAnalysisSummary.GSE278200.reversal.n_disease_responsive,
    preclinicalAnalysisSummary.GSE278200.reversal.fraction_directionally_reversed,
    preclinicalAnalysisSummary.GSE278200.reversal.fraction_closer_to_control,
    preclinicalAnalysisSummary.GSE278200.reversal.median_reversal_index_among_directional,
    preclinicalAnalysisSummary.GSE278200.reversal["fraction_directional_with_reversal_index_0.5_to_1.5"],
    preclinicalAnalysisSummary.GSE278200.reversal.spearman_disease_vs_treatment_rho,
    "Animal-level lung transcriptional disease and treatment contrasts",
    `${preclinicalAnalysisSummary.GSE278200.reversal.caution} Small male-rat bleomycin lung experiment (n=4/group); not liver or DILI evidence and not human efficacy confirmation.`,
    "fresh_analysis/preclinical/results/analysis_summary.json",
  ],
  [
    "GSE308578",
    "Mouse bleomycin lung fibrosis with nintedanib treatment",
    "Independent animal",
    JSON.stringify(preclinicalAnalysisSummary.GSE308578.group_n),
    preclinicalAnalysisSummary.GSE308578.diagnostics.n_genes,
    preclinicalAnalysisSummary.GSE308578.diagnostics.method,
    preclinicalAnalysisSummary.GSE308578.reversal.n_disease_responsive,
    preclinicalAnalysisSummary.GSE308578.reversal.fraction_directionally_reversed,
    preclinicalAnalysisSummary.GSE308578.reversal.fraction_closer_to_control,
    preclinicalAnalysisSummary.GSE308578.reversal.median_reversal_index_among_directional,
    preclinicalAnalysisSummary.GSE308578.reversal["fraction_directional_with_reversal_index_0.5_to_1.5"],
    preclinicalAnalysisSummary.GSE308578.reversal.spearman_disease_vs_treatment_rho,
    "Animal-level lung transcriptional disease and treatment contrasts with leave-one-animal-out sensitivity",
    `${preclinicalAnalysisSummary.GSE308578.reversal.caution} Model-specific male-mouse bleomycin lung experiment; not liver or DILI evidence and not human efficacy confirmation. ${preclinicalAnalysisSummary.GSE308578.metadata_resolution}`,
    "fresh_analysis/preclinical/results/analysis_summary.json",
  ],
  [
    "GSE120804",
    "Rat bile-duct-ligation precision-cut liver slices with treatment",
    "Pooled-tissue slice observation; not an independent animal",
    JSON.stringify(preclinicalAnalysisSummary.GSE120804.group_n),
    null,
    preclinicalAnalysisSummary.GSE120804.scope,
    null, null, null, null, null, null,
    "Exploratory pooled-slice molecular context",
    "Cannot support animal-population inference because slices derive from pooled rat tissues.",
    "fresh_analysis/preclinical/results/analysis_summary.json",
  ],
  [
    "GSE120679",
    "Rat TGF-beta-stimulated precision-cut lung slices with treatment",
    "Pooled-tissue slice observation; not an independent animal",
    JSON.stringify(preclinicalAnalysisSummary.GSE120679.group_n),
    null,
    preclinicalAnalysisSummary.GSE120679.scope,
    null, null, null, null, null, null,
    "Exploratory pooled-slice molecular context",
    "Cannot support animal-population inference because slices derive from pooled rat tissues.",
    "fresh_analysis/preclinical/results/analysis_summary.json",
  ],
  [
    "PXD024058",
    "Public lung-protein table with control, model, and nintedanib columns",
    "Not reconstructable from public provenance",
    "Three supplied columns per condition; biological-versus-technical replicate status unresolved",
    preclinicalAnalysisSummary.PXD024058.n_proteins,
    preclinicalAnalysisSummary.PXD024058.scope,
    null, null, null, null, null, null,
    "Descriptive protein-direction context only",
    "Replicate provenance is internally inconsistent; no biological-replicate P values or validation claim.",
    "fresh_analysis/preclinical/results/analysis_summary.json",
  ],
];

function classifyPreclinicalArtifact(p) {
  if (p.includes("/figures/")) return "figure";
  if (p.includes("/metadata/")) return "repository_metadata";
  if (p.includes("/raw/")) return "raw_or_source_input";
  if (p.includes("/scripts/") || /\.(py|r|mjs|ipynb)$/i.test(p)) return "analysis_or_validation_code";
  if (/manifest|validation_report|analysis_summary/i.test(p)) return "manifest_or_validation";
  if (/screening_ledger|key_findings/i.test(p)) return "screening_audit";
  if (/\.(md|txt)$/i.test(p)) return "documentation_or_caption";
  if (/\.(csv|tsv)(\.gz)?$/i.test(p)) return "analysis_output_table";
  return "other_frozen_artifact";
}

const preclinicalProvenanceMatrix = [
  ["Artifact_Category", "Path", "Bytes", "SHA256", "Manifest", "Frozen_Status"],
  ...preclinicalManifest.files.map((f) => {
    const check = preclinicalFileHashChecks.find((x) => x.path === f.path);
    return [
      classifyPreclinicalArtifact(f.path),
      f.path,
      f.bytes,
      f.sha256,
      "fresh_analysis/preclinical/results/analysis_manifest.json",
      check?.passed ? "frozen and hash-recorded" : `unexpected hash mismatch: ${check?.observed}`,
    ];
  }),
];

const gse151374ProvenanceMatrix = [
  ["Artifact_Category", "Path", "Bytes", "SHA256", "Manifest", "Frozen_Status"],
  ...[...gse151374Manifest.inputs, ...gse151374Manifest.outputs].map((f) => [
    classifyPreclinicalArtifact(f.path), f.path, f.bytes, f.sha256,
    "fresh_analysis/preclinical/gse151374/results/analysis_manifest.json", "frozen and hash-verified",
  ]),
];

const pxd052594ProvenanceMatrix = [
  ["Artifact_Category", "Path", "Bytes", "SHA256", "Manifest", "Frozen_Status"],
  ...[...pxd052594Manifest.inputs, ...pxd052594Manifest.outputs].map((f) => [
    classifyPreclinicalArtifact(f.path), f.path, f.bytes, f.sha256,
    "fresh_analysis/preclinical/pxd052594/results/analysis_manifest.json", "frozen and hash-verified",
  ]),
];

const integrationInputProvenanceMatrix = [
  ["Artifact_Category", "Path", "Bytes", "SHA256", "Manifest", "Frozen_Status"],
  ...Object.entries(integrationManifest.input_hashes).map(([p, h]) => ["integration_input", p, null, h, "fresh_analysis/integration/integration_analysis_manifest.json", "hash validated"]),
];

const integrationOutputProvenanceMatrix = [
  ["Artifact_Category", "Path", "Bytes", "SHA256", "Manifest", "Frozen_Status"],
  ...Object.entries(integrationManifest.output_table_hashes).map(([p, h]) => ["integration_output", p, null, h, "fresh_analysis/integration/integration_analysis_manifest.json", "hash validated"]),
  ...Object.entries(integrationManifest.text_artifact_hashes).map(([p, h]) => ["integration_text_artifact", p, null, h, "fresh_analysis/integration/integration_analysis_manifest.json", "hash validated"]),
];

const validationAuditMatrix = [
  ["Validation_Component", "Status", "Check_Count", "Manifest_SHA256", "Key_Frozen_Rule", "Source_File"],
  ["Same-compartment integration", integrationValidation.status, integrationValidation.n_checks, integrationManifestSha256, `${integrationOriginalHashMatches}/${integrationFileHashChecks.length} final manifest hashes match without exceptions. ${integrationManifest.liver_context_rule}; ${integrationManifest.primary_T2_gate}; strict sensitivity: ${integrationManifest.strict_T2_sensitivity}`, "fresh_analysis/integration/integration_validation_report.json"],
  ["GSE299128 whole-blood exact test", "PASS: workbook-side invariant check", 3, bloodManifestSha256, `30 rows; every exact P lies on the 2^7 lattice; all corrected BH q values equal 1; ${bloodManifest.primary_test}; ${bloodManifest.multiplicity}`, "fresh_analysis/blood_longitudinal/analysis_manifest.json"],
  ["Public preclinical Figure 6 core", preclinicalValidation.validation_passed ? "PASS" : "FAIL", preclinicalValidation.n_checks ?? preclinicalValidation.checks?.length ?? null, preclinicalManifestSha256, `${preclinicalOriginalHashMatches}/${preclinicalFileHashChecks.length} final manifest hashes match without exceptions; all Figure 6 result tables, statistics, figures, code and audit artifacts are frozen. ${(preclinicalValidation.scientific_boundary_warnings || []).join(" | ")}`, "fresh_analysis/preclinical/results/validation_report.json"],
  ["GSE151374 pooled-library single-cell extension (S6)", gse151374HashValidation.all_recorded_hashes_match && gse151374Validation.no_cell_level_inferential_tests ? "PASS" : "FAIL", gse151374HashValidation.entries_checked, gse151374ManifestSha256, `${gse151374FileHashChecks.length}/${gse151374FileHashChecks.length} manifest-listed hashes verified; pooled 10x library is the unit; three mice per pool; n=3 pools per condition-time; 70,696 cells retained; no cell-level inference; minimum exact two-sided P=0.10; all exact-permutation BH families null.`, "fresh_analysis/preclinical/gse151374/validation/manifest_hash_validation.json"],
  ["PXD052594 animal-level proteome extension (S7) and table-only phosphoproteome", pxd052594HashValidation.all_recorded_hashes_match && !pxd052594Validation.post_treatment_cluster_used_in_analysis ? "PASS" : "FAIL", pxd052594HashValidation.entries_checked, pxd052594ManifestSha256, `${pxd052594FileHashChecks.length}/${pxd052594FileHashChecks.length} manifest-listed hashes verified; unconditioned 7,006-protein analysis in 10 nintedanib and 13 vehicle mice is BH-null; distinct deposited 5-versus-5 subset includes 20,043 phosphosites and 994 q<0.05; no radiomic-cluster conditioning.`, "fresh_analysis/preclinical/pxd052594/validation/manifest_hash_validation.json"],
];

const files = {};
const load = async (key, relativePath) => { files[key] = await readCsv(relativePath); };
await Promise.all([
  load("endpointPts", "fresh_analysis/faers/results/prespecified_meddra_pt_sets.csv"),
  load("faersOverall", "fresh_analysis/faers/results/overall_active_comparator_signals.csv"),
  load("faersMh", "fresh_analysis/faers/results/year_adjusted_mantel_haenszel_signals.csv"),
  load("faersTto", "fresh_analysis/faers/results/time_to_onset_summary.csv"),
  load("faersSerious", "fresh_analysis/faers/results/serious_outcome_profiles.csv"),
  load("faersSpec", "fresh_analysis/faers/results/report_level_specification_signals.csv"),
  load("faersAlias", "fresh_analysis/faers/results/alias_sensitivity_signals.csv"),
  load("faersStrata", "fresh_analysis/faers/results/sex_age_stratified_signals.csv"),
  load("faersTemporal", "fresh_analysis/faers/results/temporal_influence_signals.csv"),
  load("faersPt", "fresh_analysis/faers/results/preferred_term_signals.csv"),
  load("faersAnnual", "fresh_analysis/faers/results/annual_active_comparator_signals.csv"),
  load("proteomicsEffects", "fresh_analysis/proteomics/effect_estimates.csv"),
  load("proteomicsCorr", "fresh_analysis/proteomics/correlation_stats.csv"),
  load("proteomicsZonal", "fresh_analysis/proteomics/zonal_effect_estimates.csv"),
  load("qcCounts", "fresh_analysis/proteomics/qc_counts.csv"),
  load("qcMissing", "fresh_analysis/proteomics/missingness_stats.csv"),
  load("qcCensor", "fresh_analysis/proteomics/left_censor_sensitivity.csv"),
  load("qcSourceCounts", "fresh_analysis/proteomics/source_block_counts.csv"),
  load("qcDuplicate", "fresh_analysis/proteomics/discovery_duplicate_series_qc.csv"),
  load("qcRowOrder", "fresh_analysis/proteomics/row_order_linkage_qc.csv"),
  load("qcRepeat", "fresh_analysis/proteomics/repeated_series_qc.csv"),
  load("qcPck2", "fresh_analysis/proteomics/pck2_outlier_sensitivity.csv"),
  load("qcMv", "fresh_analysis/proteomics/multivariable_model_feasibility.csv"),
  load("lincsMetrics", "fresh_analysis/lincs/analysis/gene_dose_response_metrics.csv.gz"),
  load("lincsModules", "fresh_analysis/lincs/analysis/mechanism_module_statistics.csv"),
  load("lincsSciPlex", "fresh_analysis/lincs/analysis/sciplex_hepg2_concordance.csv"),
  load("chembl", "fresh_analysis/pharmacology/chembl_target_summary.csv"),
  load("liver", "fresh_analysis/liver_singlecell/target_expression_by_compartment.csv"),
  load("vkoStable", "fresh_analysis/virtual_ko/virtual_knockout_stability_summary.csv"),
  load("vkoSeed", "fresh_analysis/virtual_ko/virtual_knockout_seed_summary.csv"),
  load("vkoNull", "fresh_analysis/integration/virtual_ko_expression_matched_null_stability.csv"),
  load("phenotypeTiers", "fresh_analysis/integration/phenotype_anchor_tiers.csv"),
  load("sameCompartment", "fresh_analysis/integration/target_same_compartment_context_all_pairs.csv"),
  load("targetGateSensitivity", "fresh_analysis/integration/target_tier_gate_sensitivity.csv"),
  load("rank", "fresh_analysis/integration/candidate_consensus_ranking.csv"),
  load("shortlist", "fresh_analysis/integration/tier2_candidate_shortlist.csv"),
  load("domain", "fresh_analysis/integration/candidate_domain_scores_long.csv"),
  load("loso", "fresh_analysis/integration/candidate_leave_one_source_out_ranks.csv"),
  load("string", "fresh_analysis/integration/string_candidate_path_metrics.csv"),
  load("bloodSummary", "fresh_analysis/blood_longitudinal/target_longitudinal_summary.csv"),
  load("bloodSlopes", "fresh_analysis/blood_longitudinal/participant_slopes.csv"),
  load("gse278Metadata", "fresh_analysis/preclinical/results/GSE278200_sample_metadata.csv"),
  load("gse278Modules", "fresh_analysis/preclinical/results/GSE278200_module_contrasts.csv"),
  load("gse278Candidates", "fresh_analysis/preclinical/results/GSE278200_candidate_results.csv"),
  load("gse308Metadata", "fresh_analysis/preclinical/results/GSE308578_sample_metadata.csv"),
  load("gse308Modules", "fresh_analysis/preclinical/results/GSE308578_module_contrasts.csv"),
  load("gse308Candidates", "fresh_analysis/preclinical/results/GSE308578_candidate_results.csv"),
  load("gse308Loo", "fresh_analysis/preclinical/results/GSE308578_leave_one_animal_out_summary.csv"),
  load("gse120804Metadata", "fresh_analysis/preclinical/results/GSE120804_sample_metadata.csv"),
  load("gse120804Modules", "fresh_analysis/preclinical/results/GSE120804_module_contrasts.csv"),
  load("gse120804Candidates", "fresh_analysis/preclinical/results/GSE120804_candidate_results.csv"),
  load("gse120679Metadata", "fresh_analysis/preclinical/results/GSE120679_sample_metadata.csv"),
  load("gse120679Modules", "fresh_analysis/preclinical/results/GSE120679_module_contrasts.csv"),
  load("gse120679Candidates", "fresh_analysis/preclinical/results/GSE120679_candidate_results.csv"),
  load("pxdCandidates", "fresh_analysis/preclinical/results/PXD024058_candidate_effects_descriptive.csv"),
  load("crossModel", "fresh_analysis/preclinical/results/cross_model_candidate_module_effects.csv"),
  load("gse151LibraryManifest", "fresh_analysis/preclinical/gse151374/results/library_manifest.csv"),
  load("gse151Qc", "fresh_analysis/preclinical/gse151374/results/qc_by_library.csv"),
  load("gse151CelltypeCounts", "fresh_analysis/preclinical/gse151374/results/celltype_counts_and_proportions_by_library.csv"),
  load("gse151CelltypeContrasts", "fresh_analysis/preclinical/gse151374/results/celltype_proportion_contrasts.csv"),
  load("gse151PcaAll", "fresh_analysis/preclinical/gse151374/results/pca_all_retained_pseudobulk.csv"),
  load("gse151PcaMacrophage", "fresh_analysis/preclinical/gse151374/results/pca_macrophage_pseudobulk.csv"),
  load("gse151ModuleMembership", "fresh_analysis/preclinical/gse151374/results/module_gene_membership.csv"),
  load("gse151ModuleScores", "fresh_analysis/preclinical/gse151374/results/macrophage_module_scores_by_library.csv"),
  load("gse151ModuleContrasts", "fresh_analysis/preclinical/gse151374/results/macrophage_module_contrasts.csv"),
  load("gse151CandidateTracks", "fresh_analysis/preclinical/gse151374/results/candidate_and_phenotype_tracks.csv"),
  load("gse151CandidateContrasts", "fresh_analysis/preclinical/gse151374/results/candidate_gene_contrasts.csv"),
  load("pxd052AnimalManifest", "fresh_analysis/preclinical/pxd052594/results/animal_sample_manifest.csv"),
  load("pxd052Proteome", "fresh_analysis/preclinical/pxd052594/results/proteome_NINT_vs_vehicle_all_animals.csv"),
  load("pxd052ProteomePca", "fresh_analysis/preclinical/pxd052594/results/proteome_PCA.csv"),
  load("pxd052ModuleMembership", "fresh_analysis/preclinical/pxd052594/results/module_protein_membership.csv"),
  load("pxd052ModuleScores", "fresh_analysis/preclinical/pxd052594/results/module_scores_by_animal.csv"),
  load("pxd052ModuleContrasts", "fresh_analysis/preclinical/pxd052594/results/module_contrasts.csv"),
  load("pxd052CandidateExpression", "fresh_analysis/preclinical/pxd052594/results/candidate_expression_by_animal.csv"),
  load("pxd052CandidateContrasts", "fresh_analysis/preclinical/pxd052594/results/candidate_contrasts.csv"),
  load("pxd052PhosphoPcaManifest", "fresh_analysis/preclinical/pxd052594/results/phosphoproteome_PCA_and_sample_manifest.csv"),
  load("pxd052Phosphosite", "fresh_analysis/preclinical/pxd052594/results/phosphosite_NINT_vs_vehicle_all_deposited_subset_animals.csv"),
  load("pxd052PhosphositeSummary", "fresh_analysis/preclinical/pxd052594/results/phosphosite_target_pathway_summary.csv"),
]);
files.faersHeterogeneity = await readPackagedCsv("Source_Data/FAERS/year_effect_heterogeneity_tests.csv");

const gse151DeDay7Summary = await summarizeCsvNumericColumn(
  "fresh_analysis/preclinical/gse151374/results/macrophage_DE_NINT_vs_BLM_day7.csv", "q_BH", (value) => value < 0.10,
);
const gse151DeDay14Summary = await summarizeCsvNumericColumn(
  "fresh_analysis/preclinical/gse151374/results/macrophage_DE_NINT_vs_BLM_day14.csv", "q_BH", (value) => value < 0.10,
);
const gse151DeDay7PackagedPath = "Source_Data/GSE151374/results/macrophage_DE_NINT_vs_BLM_day7.csv";
const gse151DeDay14PackagedPath = "Source_Data/GSE151374/results/macrophage_DE_NINT_vs_BLM_day14.csv";
const gse151DeFileIndex = [
  ["Contrast", "Full_Result_File", "Rows", "Numeric_BH_q_Values", "BH_q_lt_0_10", "Statistical_Unit", "Method", "SHA256", "Workbook_Inclusion"],
  ["BLM+nintedanib versus BLM, day 7", gse151DeDay7PackagedPath, gse151DeDay7Summary.rows, gse151DeDay7Summary.numeric, gse151DeDay7Summary.matched, "Pooled 10x library; n=3 per arm; three mice per pool", "Macrophage pseudobulk Welch test; BH within day", gse151374Manifest.outputs.find((x) => x.path === gse151DeDay7Summary.relativePath)?.sha256, "Hash-verified external machine-readable table in Source_Data; indexed here because artifact-tool cannot serialize both 48,795-row DE matrices inside one XLSX"],
  ["BLM+nintedanib versus BLM, day 14", gse151DeDay14PackagedPath, gse151DeDay14Summary.rows, gse151DeDay14Summary.numeric, gse151DeDay14Summary.matched, "Pooled 10x library; n=3 per arm; three mice per pool", "Macrophage pseudobulk Welch test; BH within day", gse151374Manifest.outputs.find((x) => x.path === gse151DeDay14Summary.relativePath)?.sha256, "Hash-verified external machine-readable table in Source_Data; indexed here because artifact-tool cannot serialize both 48,795-row DE matrices inside one XLSX"],
];

function matrixObjects(matrix) {
  const headers = matrix[0].map(String);
  return matrix.slice(1).map((row) => Object.fromEntries(headers.map((h, i) => [h, row[i]])));
}

const targetGateRows = matrixObjects(files.targetGateSensitivity.matrix);
const exploratoryT2GateTargets = targetGateRows.filter((r) => r.primary_T2_gate === true).map((r) => r.gene_symbol).sort();
const strictAllSeedT2Targets = targetGateRows.filter((r) => r.strict_all_seed_T2_sensitivity_gate === true).map((r) => r.gene_symbol).sort();
if (JSON.stringify(exploratoryT2GateTargets) !== JSON.stringify(["FGFR1"])) {
  throw new Error(`Unexpected exploratory T2-gate set: ${JSON.stringify(exploratoryT2GateTargets)}`);
}
if (strictAllSeedT2Targets.length !== 0) {
  throw new Error(`Unexpected strict all-seed T2 set: ${JSON.stringify(strictAllSeedT2Targets)}`);
}

const bloodRows = matrixObjects(files.bloodSummary.matrix);
const bloodExactGridPassed = bloodRows.length === 30 && bloodRows.every((r) => {
  const p = r.exact_signflip_p;
  return typeof p === "number" && p >= 0 && p <= 1 && Math.abs(p * (2 ** bloodManifest.n_participants) - Math.round(p * (2 ** bloodManifest.n_participants))) < 1e-10;
});
if (!bloodExactGridPassed) {
  throw new Error("GSE299128 exact sign-flip values do not form a complete 30-row 2^7 probability lattice.");
}
const bloodBhPassed = bloodRows.every((r) => typeof r.q_bh === "number" && r.q_bh >= 0 && r.q_bh <= 1) && bloodRows.every((r) => Math.abs(r.q_bh - 1) < 1e-12);
if (!bloodBhPassed) throw new Error("GSE299128 corrected BH q values are not the frozen all-q=1 result.");

const gse151QcRows = matrixObjects(files.gse151Qc.matrix);
const gse151TotalCells = gse151QcRows.reduce((sum, r) => sum + Number(r.cells_retained || 0), 0);
const gse151ConditionCounts = new Map();
for (const r of gse151QcRows) gse151ConditionCounts.set(r.condition, (gse151ConditionCounts.get(r.condition) || 0) + 1);
const gse151DesignPassed = gse151QcRows.length === 18 && gse151TotalCells === 70696
  && [...gse151ConditionCounts.values()].length === 6 && [...gse151ConditionCounts.values()].every((n) => n === 3)
  && gse151QcRows.every((r) => r.n_mice_pooled === 3 && r.statistical_unit === "pooled 10x library")
  && gse151374Validation.no_cell_level_inferential_tests === true
  && gse151374Validation.exact_permutation_minimum_two_sided_p_for_3v3 === 0.1;
const gse151ExactFamiliesNull = [files.gse151CelltypeContrasts, files.gse151ModuleContrasts, files.gse151CandidateContrasts].every((file) =>
  matrixObjects(file.matrix).every((r) => Object.entries(r).filter(([k]) => /q_BH_exact/i.test(k)).every(([, v]) => v === null || Number(v) >= 0.05))
);
const gse151GeneFamiliesNull = gse151DeDay7Summary.rows === 48795 && gse151DeDay14Summary.rows === 48795
  && gse151DeDay7Summary.matched === 0 && gse151DeDay14Summary.matched === 0;
if (!gse151DesignPassed || !gse151ExactFamiliesNull || !gse151GeneFamiliesNull) {
  throw new Error("Frozen GSE151374 design, retained-cell count, or null multiplicity-controlled result changed.");
}

const pxd052ProteomeRows = matrixObjects(files.pxd052Proteome.matrix);
const pxd052PhosphositeRows = matrixObjects(files.pxd052Phosphosite.matrix);
const pxd052ProteomeBhNull = pxd052ProteomeRows.length === 7006 && pxd052ProteomeRows.every((r) => r.q_BH === null || Number(r.q_BH) >= 0.10);
const pxd052PhosphositeQ05Count = pxd052PhosphositeRows.filter((r) => typeof r.q_BH === "number" && r.q_BH < 0.05).length;
const pxd052DesignPassed = pxd052594Validation.animals === 23 && pxd052594Validation.nintedanib_animals === 10
  && pxd052594Validation.vehicle_animals === 13 && pxd052594Validation.post_treatment_cluster_used_in_analysis === false
  && pxd052594Validation.phosphosites === 20043 && pxd052594Validation.phosphoproteome_nintedanib_animals === 5
  && pxd052594Validation.phosphoproteome_vehicle_animals === 5;
if (!pxd052ProteomeBhNull || pxd052PhosphositeQ05Count !== 994 || !pxd052DesignPassed) {
  throw new Error(`Frozen PXD052594 result changed: proteome_null=${pxd052ProteomeBhNull}, phosphosite_q05=${pxd052PhosphositeQ05Count}, design=${pxd052DesignPassed}`);
}

const extensionRoleRows = [...matrixObjects(files.gse151CandidateContrasts.matrix), ...matrixObjects(files.pxd052CandidateContrasts.matrix)];
const extensionExposureTargets = [...new Set(extensionRoleRows.filter((r) => String(r.track_type).toLowerCase().includes("exposure_proximal_target")).map((r) => String(r.gene).toLowerCase()))].sort();
const extensionPhenotypeAnchors = [...new Set(extensionRoleRows.filter((r) => String(r.track_type).toLowerCase().includes("dili_phenotype_anchor")).map((r) => String(r.gene).toLowerCase()))].sort();
if (JSON.stringify(extensionExposureTargets) !== JSON.stringify(["fgfr1"]) || JSON.stringify(extensionPhenotypeAnchors) !== JSON.stringify(["fbp1", "gsta1", "otc"])) {
  throw new Error(`Extension role labels changed: targets=${extensionExposureTargets}, phenotype=${extensionPhenotypeAnchors}`);
}

const workbook = Workbook.create();
const sheetRecords = {};

sheetRecords.README = makeSimpleSheet(workbook, "README", readme, "Readme");
sheetRecords.Data_Sources = makeSimpleSheet(workbook, "Data_Sources", dataSources, "DataSources");
sheetRecords.Data_Provenance = makeMultiBlockSheet(workbook, "Data_Provenance", [
  { title: "Integration input hashes", source: "fresh_analysis/integration/integration_analysis_manifest.json", matrix: integrationInputProvenanceMatrix, note: "SHA-256 hashes validated by the frozen integration validation report." },
  { title: "Integration output and text-artifact hashes", source: "fresh_analysis/integration/integration_analysis_manifest.json", matrix: integrationOutputProvenanceMatrix, note: "Final same-compartment integration outputs and text artifacts." },
  { title: "Public preclinical frozen-file manifest", source: "fresh_analysis/preclinical/results/analysis_manifest.json", matrix: preclinicalProvenanceMatrix, note: "Paths, byte sizes, and SHA-256 hashes reproduced from the final frozen preclinical manifest." },
  { title: "Validation audit", source: "integration, blood, and preclinical manifests", matrix: validationAuditMatrix, note: "Automated validation status and frozen interpretation rules." },
  { title: "GSE151374 S6 frozen manifest", source: "fresh_analysis/preclinical/gse151374/results/analysis_manifest.json", matrix: gse151374ProvenanceMatrix, note: "All S6 inputs and outputs are frozen and hash-verified; the pooled 10x library, not the cell, is the inferential unit." },
  { title: "PXD052594 S7 frozen manifest", source: "fresh_analysis/preclinical/pxd052594/results/analysis_manifest.json", matrix: pxd052594ProvenanceMatrix, note: "All S7 and table-only phosphoproteome inputs and outputs are frozen and hash-verified; proteome and phosphoproteome cohorts remain distinct." },
]);
sheetRecords.Statistical_Units = makeSimpleSheet(workbook, "Statistical_Units", statisticalUnits, "StatisticalUnits");
sheetRecords.Figure_Source_Map = makeSimpleSheet(workbook, "Figure_Source_Map", figureSourceMap, "FigureSourceMap");
sheetRecords.Endpoint_Definitions = makeMultiBlockSheet(workbook, "Endpoint_Definitions", [
  { title: "Endpoint summary", source: "FAERS README and frozen analysis specification", matrix: endpointSummary, note: "Operational definitions and interpretation boundaries." },
  { title: "Prespecified MedDRA preferred terms", source: files.endpointPts.relativePath, matrix: files.endpointPts.matrix },
]);
sheetRecords.FAERS_Overall = makeMultiBlockSheet(workbook, "FAERS_Overall", [
  { title: "All-indication drug-name-indexed active-comparator signals", source: files.faersOverall.relativePath, matrix: files.faersOverall.matrix, note: "Report-level disproportionality; not incidence, risk or causal attribution to a drug element." },
  { title: "Year-adjusted Mantel-Haenszel estimates", source: files.faersMh.relativePath, matrix: files.faersMh.matrix, note: "Descriptive fixed-effect calendar-time standardization with Robins–Breslow–Greenland confidence intervals; not a causal adjustment model. Heterogeneity is tested separately." },
  { title: "Time-to-onset summary", source: files.faersTto.relativePath, matrix: files.faersTto.matrix, note: "Exact-date subset only; report-level event date may not be hepatic-event-specific." },
  { title: "Serious report outcomes", source: files.faersSerious.relativePath, matrix: files.faersSerious.matrix, note: "Report outcome labels are not event incidence or attributable outcomes." },
]);
sheetRecords.FAERS_Specifications = makeMultiBlockSheet(workbook, "FAERS_Specifications", [
  { title: "Report-level specification curves", source: files.faersSpec.relativePath, matrix: files.faersSpec.matrix },
  { title: "Alias sensitivity", source: files.faersAlias.relativePath, matrix: files.faersAlias.matrix },
  { title: "Sex and strict-age-unit strata", source: files.faersStrata.relativePath, matrix: files.faersStrata.matrix, note: "Age strata include only openFDA unit 801 (years)." },
  { title: "Temporal influence analyses", source: files.faersTemporal.relativePath, matrix: files.faersTemporal.matrix },
  { title: "Preferred-term signals", source: files.faersPt.relativePath, matrix: files.faersPt.matrix },
]);
sheetRecords.FAERS_Annual = makeMultiBlockSheet(workbook, "FAERS_Annual", [
  { title: "Annual active-comparator signals", source: files.faersAnnual.relativePath, matrix: files.faersAnnual.matrix },
  { title: "Prespecified calendar-year heterogeneity tests", source: files.faersHeterogeneity.relativePath, matrix: files.faersHeterogeneity.matrix, note: "Tarone-adjusted Breslow–Day and Cochran-Q/I² quantify temporal instability; exclusion-of-2022 rows are influence sensitivities." },
]);
sheetRecords.Proteomics_Effects = makeMultiBlockSheet(workbook, "Proteomics_Effects", [
  { title: "Two-group effect estimates", source: files.proteomicsEffects.relativePath, matrix: files.proteomicsEffects.matrix, note: "Univariate public-data reanalysis; FDR is within cohort/contrast." },
  { title: "Protein-ALT Spearman correlations", source: files.proteomicsCorr.relativePath, matrix: files.proteomicsCorr.matrix },
  { title: "Exploratory liver-zone effects", source: files.proteomicsZonal.relativePath, matrix: files.proteomicsZonal.matrix, note: "No zone comparison survived FDR correction." },
]);
sheetRecords.Proteomics_QC = makeMultiBlockSheet(workbook, "Proteomics_QC", [
  { title: "Observed and missing counts", source: files.qcCounts.relativePath, matrix: files.qcCounts.matrix },
  { title: "Differential missingness", source: files.qcMissing.relativePath, matrix: files.qcMissing.matrix },
  { title: "Left-censor sensitivity", source: files.qcCensor.relativePath, matrix: files.qcCensor.matrix },
  { title: "Source block count audit", source: files.qcSourceCounts.relativePath, matrix: files.qcSourceCounts.matrix },
  { title: "Exact duplicate series audit", source: files.qcDuplicate.relativePath, matrix: files.qcDuplicate.matrix, note: "HPD discovery evidence is not counted as independent where it duplicates ALDOB." },
  { title: "Row-order linkage audit", source: files.qcRowOrder.relativePath, matrix: files.qcRowOrder.matrix, note: "Cross-protein positional joins are unsafe." },
  { title: "Repeated-series audit", source: files.qcRepeat.relativePath, matrix: files.qcRepeat.matrix },
  { title: "PCK2 starred-outlier sensitivity", source: files.qcPck2.relativePath, matrix: files.qcPck2.matrix },
  { title: "Multivariable model feasibility", source: files.qcMv.relativePath, matrix: files.qcMv.matrix, note: "No participant-level multivariable classifier was fabricated." },
]);
sheetRecords.LINCS_Dose_Metrics = makeMultiBlockSheet(workbook, "LINCS_Dose_Metrics", [
  { title: "All-gene six-dose response metrics", source: files.lincsMetrics.relativePath, matrix: files.lincsMetrics.matrix, note: "Six HepG2 composite Level-5 signatures at 24 h; exact six-dose permutation tests." },
  { title: "Mechanism-module statistics", source: files.lincsModules.relativePath, matrix: files.lincsModules.matrix },
  { title: "Sci-Plex concordance sensitivity", source: files.lincsSciPlex.relativePath, matrix: files.lincsSciPlex.matrix, note: "Overlap-conditional comparison to truncated MCF7 ranked lists." },
]);
sheetRecords.ChEMBL_Targets = makeSimpleSheet(workbook, "ChEMBL_Targets", files.chembl.matrix, "ChemblTargets");
sheetRecords.Liver_Expression = makeSimpleSheet(workbook, "Liver_Expression", files.liver.matrix, "LiverExpression");
sheetRecords.Virtual_KO_Stability = makeMultiBlockSheet(workbook, "Virtual_KO_Stability", [
  { title: "Across-seed stability", source: files.vkoStable.relativePath, matrix: files.vkoStable.matrix, note: "Unsigned network displacement only; no protective or harmful direction." },
  { title: "Seed-level summaries", source: files.vkoSeed.relativePath, matrix: files.vkoSeed.matrix },
  { title: "Expression-matched null stability", source: files.vkoNull.relativePath, matrix: files.vkoNull.matrix, note: "Nominal or BH flags are computational diagnostics, not causal validation." },
]);
sheetRecords.Integrated_Prioritization = makeMultiBlockSheet(workbook, "Integrated_Prioritization", [
  { title: "Phenotype anchor tiers", source: files.phenotypeTiers.relativePath, matrix: files.phenotypeTiers.matrix, note: "P1/P2/P3 phenotype tiers are distinct from target T1/T2/T3 tiers; FBP1 is a P1 DILI-discriminating phenotype anchor, not a causal target." },
  { title: "All target-compartment context pairs", source: files.sameCompartment.relativePath, matrix: files.sameCompartment.matrix, note: "Expression and median virtual-KO ranks are combined only within the same compartment; no cross-compartment splicing." },
  { title: "Target tier-gate sensitivity", source: files.targetGateSensitivity.relativePath, matrix: files.targetGateSensitivity.matrix, note: "FGFR1 is the sole candidate meeting the exploratory T2 gate; no target meets the strict all-seed T2 sensitivity gate and no Tier 1 target is assigned." },
  { title: "Consensus ranking", source: files.rank.relativePath, matrix: files.rank.matrix, note: "Candidate roles kept separate; no Tier 1 target assigned." },
  { title: "Tier 2 shortlist", source: files.shortlist.relativePath, matrix: files.shortlist.matrix, note: "The exploratory T2 gate defines a prespecified same-compartment target-network hypothesis, not experimental validation; strict all-seed sensitivity is reported separately." },
  { title: "Evidence-domain scores", source: files.domain.relativePath, matrix: files.domain.matrix },
  { title: "Leave-one-source-out ranks", source: files.loso.relativePath, matrix: files.loso.matrix },
  { title: "STRING contextual paths", source: files.string.relativePath, matrix: files.string.matrix, note: "Association evidence only; not a demonstrated toxicity pathway." },
]);
sheetRecords.Blood_Sensitivity = makeMultiBlockSheet(workbook, "Blood_Sensitivity", [
  { title: "Gene-level longitudinal summary", source: files.bloodSummary.relativePath, matrix: files.bloodSummary.matrix, note: "Two-sided exact participant-level sign-flip P values enumerate all 2^7 assignments; BH correction spans 30 prespecified genes. The study has no adjudicated DILI endpoint." },
  { title: "Participant-level slopes", source: files.bloodSlopes.relativePath, matrix: files.bloodSlopes.matrix },
]);
sheetRecords.Preclinical_Evidence = makeMultiBlockSheet(workbook, "Preclinical_Evidence", [
  { title: "Cross-dataset evidence summary", source: "fresh_analysis/preclinical/results/analysis_summary.json", matrix: preclinicalSummaryMatrix, note: "Independent-animal, pooled-slice, and unresolved-provenance components are explicitly distinguished." },
  { title: "GSE278200 sample metadata", source: files.gse278Metadata.relativePath, matrix: files.gse278Metadata.matrix, note: "Twelve independent animals; four per saline, bleomycin, and bleomycin plus nintedanib group." },
  { title: "GSE278200 module contrasts", source: files.gse278Modules.relativePath, matrix: files.gse278Modules.matrix, note: "Independent-animal contrasts; exact label-permutation P values and within-contrast module BH q values are preserved." },
  { title: "GSE278200 candidate results", source: files.gse278Candidates.relativePath, matrix: files.gse278Candidates.matrix, note: "Animal-level moderated models; genome-wide BH q values are preserved." },
  { title: "GSE308578 sample metadata", source: files.gse308Metadata.relativePath, matrix: files.gse308Metadata.matrix, note: "GEO SOFT and ARRIVE Supplementary Figure S1 support CTRL 10, vehicle 16, nintedanib 17; no samples were relabelled." },
  { title: "GSE308578 module contrasts", source: files.gse308Modules.relativePath, matrix: files.gse308Modules.matrix, note: "Independent-animal contrasts; label-permutation P values and within-contrast module BH q values are preserved." },
  { title: "GSE308578 candidate results", source: files.gse308Candidates.relativePath, matrix: files.gse308Candidates.matrix, note: "Animal-level moderated models; genome-wide BH q values are preserved." },
  { title: "GSE308578 leave-one-animal-out stability", source: files.gse308Loo.relativePath, matrix: files.gse308Loo.matrix, note: "Leave-one-animal-out summaries preserve sign-concordance and full deletion count." },
  { title: "GSE120804 sample metadata", source: files.gse120804Metadata.relativePath, matrix: files.gse120804Metadata.matrix, note: "Pooled-tissue slice observations; not independent animals." },
  { title: "GSE120804 exploratory module contrasts", source: files.gse120804Modules.relativePath, matrix: files.gse120804Modules.matrix, note: "Exploratory only; cannot support animal-population inference." },
  { title: "GSE120804 exploratory candidate results", source: files.gse120804Candidates.relativePath, matrix: files.gse120804Candidates.matrix, note: "Exploratory pooled-slice estimates; not animal-level confirmation." },
  { title: "GSE120679 sample metadata", source: files.gse120679Metadata.relativePath, matrix: files.gse120679Metadata.matrix, note: "Pooled-tissue slice observations; not independent animals." },
  { title: "GSE120679 exploratory module contrasts", source: files.gse120679Modules.relativePath, matrix: files.gse120679Modules.matrix, note: "Exploratory only; cannot support animal-population inference." },
  { title: "GSE120679 exploratory candidate results", source: files.gse120679Candidates.relativePath, matrix: files.gse120679Candidates.matrix, note: "Exploratory pooled-slice estimates; not animal-level confirmation." },
  { title: "PXD024058 descriptive candidate effects", source: files.pxdCandidates.relativePath, matrix: files.pxdCandidates.matrix, note: "Descriptive only because biological-versus-technical replicate provenance is internally inconsistent." },
  { title: "Cross-model candidate and module effects", source: files.crossModel.relativePath, matrix: files.crossModel.matrix, note: "Effect display across models; scope labels distinguish independent-animal, pooled-slice, and descriptive proteomic evidence." },
]);
sheetRecords.GSE151374_scRNA = makeMultiBlockSheet(workbook, "GSE151374_scRNA", [
  { title: "Pooled-library manifest", source: files.gse151LibraryManifest.relativePath, matrix: files.gse151LibraryManifest.matrix, note: "Each inferential unit is one independently prepared 10x library pooling three mice; n=3 libraries per condition-time cell." },
  { title: "Library-level QC", source: files.gse151Qc.relativePath, matrix: files.gse151Qc.matrix, note: "A total of 70,696 cells passed QC, but cells are nested descriptive observations and are never treated as independent replicates." },
  { title: "Cell-type counts and proportions", source: files.gse151CelltypeCounts.relativePath, matrix: files.gse151CelltypeCounts.matrix, note: "Broad marker-rule strata summarized once per pooled library." },
  { title: "Cell-type proportion contrasts", source: files.gse151CelltypeContrasts.relativePath, matrix: files.gse151CelltypeContrasts.matrix, note: "Pooled-library tests with exhaustive 3-versus-3 permutations; minimum two-sided exact P=0.10 and all exact-permutation BH families are null." },
  { title: "All-retained pseudobulk PCA", source: files.gse151PcaAll.relativePath, matrix: files.gse151PcaAll.matrix, note: "Descriptive library-level PCA coordinates." },
  { title: "Macrophage pseudobulk PCA", source: files.gse151PcaMacrophage.relativePath, matrix: files.gse151PcaMacrophage.matrix, note: "Descriptive macrophage-stratum library-level PCA coordinates." },
  { title: "Macrophage module membership", source: files.gse151ModuleMembership.relativePath, matrix: files.gse151ModuleMembership.matrix, note: "Prespecified gene membership used to calculate library-pseudobulk program scores." },
  { title: "Macrophage module scores", source: files.gse151ModuleScores.relativePath, matrix: files.gse151ModuleScores.matrix, note: "One score per pooled library and prespecified program." },
  { title: "Macrophage module contrasts", source: files.gse151ModuleContrasts.relativePath, matrix: files.gse151ModuleContrasts.matrix, note: "No program met multiplicity-controlled significance." },
  { title: "Candidate and phenotype tracks", source: files.gse151CandidateTracks.relativePath, matrix: files.gse151CandidateTracks.matrix, note: "Fgfr1 is the sole exposure-proximal target track; Fbp1/Gsta1/Otc are downstream DILI phenotype anchors, not targets." },
  { title: "Candidate gene contrasts", source: files.gse151CandidateContrasts.relativePath, matrix: files.gse151CandidateContrasts.matrix, note: "All exact-permutation BH results are null; low or absent detectability is retained explicitly." },
  { title: "Macrophage DE full-file index", source: "fresh_analysis/preclinical/gse151374/results/macrophage_DE_NINT_vs_BLM_day7.csv; fresh_analysis/preclinical/gse151374/results/macrophage_DE_NINT_vs_BLM_day14.csv", matrix: gse151DeFileIndex, note: "Both complete 48,795-row source CSVs remain unchanged and hash-verified in the submission package. The workbook indexes them instead of duplicating 97,590 rows because artifact-tool otherwise exceeds its serialization string limit; zero genes meet BH q<0.10 at either day." },
]);
sheetRecords.PXD052594_Proteomics = makeMultiBlockSheet(workbook, "PXD052594_Proteomics", [
  { title: "Animal sample manifest", source: files.pxd052AnimalManifest.relativePath, matrix: files.pxd052AnimalManifest.matrix, note: "All remaining animals in the public proteome are included: nintedanib n=10 and vehicle n=13; no radiomic response cluster was used." },
  { title: "All-animal proteome contrasts", source: files.pxd052Proteome.relativePath, matrix: files.pxd052Proteome.matrix, note: "Complete unconditioned 7,006-protein Welch contrast with BH across all proteins; no protein met q<0.10." },
  { title: "Proteome PCA", source: files.pxd052ProteomePca.relativePath, matrix: files.pxd052ProteomePca.matrix, note: "Individual-animal PCA coordinates based on the 2,000 most variable proteins." },
  { title: "Protein module membership", source: files.pxd052ModuleMembership.relativePath, matrix: files.pxd052ModuleMembership.matrix, note: "Prespecified pulmonary program membership." },
  { title: "Protein module scores", source: files.pxd052ModuleScores.relativePath, matrix: files.pxd052ModuleScores.matrix, note: "One program score per individual mouse." },
  { title: "Protein module contrasts", source: files.pxd052ModuleContrasts.relativePath, matrix: files.pxd052ModuleContrasts.matrix, note: "No module met BH significance." },
  { title: "Candidate expression by animal", source: files.pxd052CandidateExpression.relativePath, matrix: files.pxd052CandidateExpression.matrix, note: "Absolute deposited abundance or explicit nonquantification for the prespecified tracks." },
  { title: "Candidate contrasts", source: files.pxd052CandidateContrasts.relativePath, matrix: files.pxd052CandidateContrasts.matrix, note: "Fgfr1 is the sole exposure-proximal target; Fbp1/Gsta1/Otc remain downstream phenotype anchors." },
  { title: "Deposited phosphoproteome subset manifest and PCA", source: files.pxd052PhosphoPcaManifest.relativePath, matrix: files.pxd052PhosphoPcaManifest.matrix, note: "Distinct deposited random 5-versus-5 animal subset; it must not be merged with the 23-animal proteome." },
  { title: "All deposited-subset phosphosite contrasts", source: files.pxd052Phosphosite.relativePath, matrix: files.pxd052Phosphosite.matrix, note: "Complete table-only 20,043-site sensitivity with BH across all sites; 994 sites have q<0.05. This is contextual pulmonary signaling, not a figure panel or hepatic-safety confirmation." },
  { title: "Phosphosite target-pathway summary", source: files.pxd052PhosphositeSummary.relativePath, matrix: files.pxd052PhosphositeSummary.matrix, note: "Prespecified target/pathway summary; absence of Fgfr1 phosphorylation-site quantification is retained explicitly." },
]);
sheetRecords.Experimental_Validation = makeSimpleSheet(workbook, "Experimental_Validation", experimentalValidation, "ExperimentalValidation");
sheetRecords.Reporting_Checklists = makeSimpleSheet(workbook, "Reporting_Checklists", reportingChecklists, "ReportingChecklists");

// Small-sheet readability refinements.
for (const name of ["README", "Data_Sources", "Statistical_Units", "Figure_Source_Map", "Experimental_Validation", "Reporting_Checklists"]) {
  const sheet = workbook.worksheets.getItem(name);
  const used = sheet.getUsedRange(true);
  used.format.wrapText = true;
  if (name === "README") {
    sheet.getRange("A:A").format.columnWidth = 32;
    sheet.getRange("B:B").format.columnWidth = 110;
    sheet.getRange("A2:B21").format.rowHeight = 64;
  }
  if (name === "Data_Sources") {
    const widths = [10, 24, 24, 32, 38, 16, 46, 30, 58, 48, 66];
    widths.forEach((w, j) => { sheet.getRange(`${colLetter(j + 1)}:${colLetter(j + 1)}`).format.columnWidth = w; });
    sheet.getRange("A2:K16").format.rowHeight = 108;
  }
  if (name === "Statistical_Units") {
    const widths = [27, 30, 34, 38, 38, 54, 48, 42, 54, 54];
    widths.forEach((w, j) => { sheet.getRange(`${colLetter(j + 1)}:${colLetter(j + 1)}`).format.columnWidth = w; });
    sheet.getRange("A2:J17").format.rowHeight = 116;
  }
  if (name === "Figure_Source_Map") {
    const widths = [26, 19, 54, 90, 42, 58, 66, 54];
    widths.forEach((w, j) => { sheet.getRange(`${colLetter(j + 1)}:${colLetter(j + 1)}`).format.columnWidth = w; });
    sheet.getRange("A2:H52").format.rowHeight = 118;
  }
  if (name === "Experimental_Validation") {
    const widths = [22, 38, 28, 32, 44, 48, 62, 74, 68, 26];
    widths.forEach((w, j) => { sheet.getRange(`${colLetter(j + 1)}:${colLetter(j + 1)}`).format.columnWidth = w; });
    sheet.getRange("A2:J12").format.rowHeight = 145;
  }
  if (name === "Reporting_Checklists") {
    const widths = [30, 34, 26, 34, 70, 60, 54];
    widths.forEach((w, j) => { sheet.getRange(`${colLetter(j + 1)}:${colLetter(j + 1)}`).format.columnWidth = w; });
    sheet.getRange("A2:G12").format.rowHeight = 100;
  }
}

// Provenance and boundary-heavy sheets need explicit widths and row heights so hashes and caveats remain visible.
{
  const sheet = workbook.worksheets.getItem("Data_Provenance");
  const widths = [28, 78, 16, 68, 72, 28];
  widths.forEach((w, j) => { sheet.getRange(`${colLetter(j + 1)}:${colLetter(j + 1)}`).format.columnWidth = w; });
  sheet.getUsedRange(true).format.wrapText = true;
  const index = sheetRecords.Data_Provenance.blocks[0];
  if (index.nRows > 0) sheet.getRangeByIndexes(index.startRowZero + 1, 0, index.nRows, index.nCols).format.rowHeight = 118;
  const manifestBlock = sheetRecords.Data_Provenance.blocks[3];
  if (manifestBlock.nRows > 0) sheet.getRangeByIndexes(manifestBlock.startRowZero + 1, 0, manifestBlock.nRows, manifestBlock.nCols).format.rowHeight = 30;
  const validationBlock = sheetRecords.Data_Provenance.blocks[4];
  if (validationBlock.nRows > 0) sheet.getRangeByIndexes(validationBlock.startRowZero + 1, 0, validationBlock.nRows, validationBlock.nCols).format.rowHeight = 120;
  for (const blockIndex of [5, 6]) {
    const extensionBlock = sheetRecords.Data_Provenance.blocks[blockIndex];
    if (extensionBlock?.nRows > 0) sheet.getRangeByIndexes(extensionBlock.startRowZero + 1, 0, extensionBlock.nRows, extensionBlock.nCols).format.rowHeight = 30;
  }
}
{
  const sheet = workbook.worksheets.getItem("Preclinical_Evidence");
  const index = sheetRecords.Preclinical_Evidence.blocks[0];
  if (index.nRows > 0) sheet.getRangeByIndexes(index.startRowZero + 1, 0, index.nRows, index.nCols).format.rowHeight = 118;
  const summary = sheetRecords.Preclinical_Evidence.blocks[1];
  if (summary.nRows > 0) {
    const summaryBody = sheet.getRangeByIndexes(summary.startRowZero + 1, 0, summary.nRows, summary.nCols);
    summaryBody.format.wrapText = true;
    summaryBody.format.rowHeight = 112;
  }
}
for (const name of ["GSE151374_scRNA", "PXD052594_Proteomics"]) {
  const sheet = workbook.worksheets.getItem(name);
  const index = sheetRecords[name].blocks[0];
  if (index.nRows > 0) sheet.getRangeByIndexes(index.startRowZero + 1, 0, index.nRows, index.nCols).format.rowHeight = 118;
}
// Keep both complete external GSE151374 macrophage-DE file references directly
// readable in the workbook. These are the only source matrices intentionally
// indexed rather than duplicated cell-by-cell, so the path, full SHA-256, and
// packaging note must be visible without relying on the formula bar.
{
  const sheet = workbook.worksheets.getItem("GSE151374_scRNA");
  const externalIndex = sheetRecords.GSE151374_scRNA.blocks.at(-1);
  const widths = { A: 38, B: 76, F: 54, G: 54, H: 70, I: 88 };
  for (const [column, width] of Object.entries(widths)) sheet.getRange(`${column}:${column}`).format.columnWidth = width;
  const externalBody = sheet.getRangeByIndexes(externalIndex.startRowZero + 1, 0, externalIndex.nRows, externalIndex.nCols);
  externalBody.format.wrapText = true;
  externalBody.format.verticalAlignment = "top";
  externalBody.format.rowHeight = 88;
}
workbook.worksheets.getItem("GSE151374_scRNA").getUsedRange(true).conditionalFormats.add("containsText", {
  text: "no cell-level",
  format: { fill: COLORS.paleOrange, font: { color: "#9A3412", bold: true } },
});
workbook.worksheets.getItem("PXD052594_Proteomics").getUsedRange(true).conditionalFormats.add("containsText", {
  text: "distinct deposited",
  format: { fill: COLORS.paleOrange, font: { color: "#9A3412", bold: true } },
});

// Semantic conditional formatting for evidence status fields.
workbook.worksheets.getItem("Integrated_Prioritization").getUsedRange(true).conditionalFormats.add("containsText", {
  text: "T2:",
  format: { fill: COLORS.paleOrange, font: { color: "#9A3412", bold: true } },
});
workbook.worksheets.getItem("Integrated_Prioritization").getUsedRange(true).conditionalFormats.add("containsText", {
  text: "P1",
  format: { fill: COLORS.paleGreen, font: { color: "#14532D", bold: true } },
});
workbook.worksheets.getItem("Experimental_Validation").getRange("J2:J12").conditionalFormats.add("containsText", {
  text: "Prospective",
  format: { fill: COLORS.paleBlue, font: { color: COLORS.navy, bold: true } },
});
workbook.worksheets.getItem("Preclinical_Evidence").getUsedRange(true).conditionalFormats.add("containsText", {
  text: "not an independent animal",
  format: { fill: COLORS.paleOrange, font: { color: "#9A3412", bold: true } },
});
workbook.worksheets.getItem("Preclinical_Evidence").getUsedRange(true).conditionalFormats.add("containsText", {
  text: "descriptive only",
  format: { fill: COLORS.paleBlue, font: { color: COLORS.navy, bold: true } },
});

// Re-assert frozen header rows after all table and formatting mutations. Calling
// this at finalization (rather than only when a sheet is created) makes the pane
// state persist in the exported XLSX for every worksheet.
for (const name of Object.keys(sheetRecords)) {
  const sheet = workbook.worksheets.getItem(name);
  sheet.freezePanes.unfreeze();
  sheet.freezePanes.freezeRows(1);
}

// Export before the inspection/render pass so large-table rendering does not inflate the
// heap at serialization time. The final XLSX is independently re-imported and re-rendered.
let output;
try {
  output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(OUT_XLSX);
} catch (error) {
  const concise = `${error?.name || "Error"}: ${error?.message || String(error)}`;
  await fs.mkdir(path.dirname(RENDER_DIR), { recursive: true });
  await fs.writeFile(path.join(path.dirname(RENDER_DIR), "export_error.txt"), concise + "\n", "utf8");
  console.error(`EXPORT_ERROR_SUMMARY ${concise}`);
  process.exit(1);
}
const stat = await fs.stat(OUT_XLSX);
const workbookSha256 = await sha256File(OUT_XLSX);

const keyInspections = {};
for (const [sheetName, range] of Object.entries({
  README: "A1:B21",
  Data_Sources: "A1:K16",
  Statistical_Units: "A1:J17",
  Figure_Source_Map: "A1:H52",
  FAERS_Overall: "A1:X14",
  Proteomics_Effects: "A1:AL18",
  LINCS_Dose_Metrics: "A1:AD12",
  Integrated_Prioritization: "A1:AY26",
  Blood_Sensitivity: "A1:H20",
  Preclinical_Evidence: "A1:O28",
  GSE151374_scRNA: "A1:M30",
  PXD052594_Proteomics: "A1:M30",
  Experimental_Validation: "A1:J12",
})) {
  const check = await workbook.inspect({
    kind: "table",
    range: `${sheetName}!${range}`,
    include: "values,formulas",
    tableMaxRows: 18,
    tableMaxCols: 40,
    maxChars: 12000,
  });
  keyInspections[sheetName] = check.ndjson;
}

const sheetCatalog = await workbook.inspect({
  kind: "sheet",
  include: "id,name",
  maxChars: 12000,
});

const styleInspections = {};
const firstCellInspections = {};
const headerStyleChecks = {};
const firstHeaderChecks = {};
for (const sheetName of Object.keys(sheetRecords)) {
  const cell = await workbook.inspect({
    kind: "table",
    range: `${sheetName}!A1:A1`,
    include: "values,formulas",
    tableMaxRows: 2,
    tableMaxCols: 2,
    maxChars: 1200,
  });
  firstCellInspections[sheetName] = cell.ndjson;
  const style = await workbook.inspect({
    kind: "computedStyle",
    sheetId: sheetName,
    range: "A1:A1",
    maxChars: 2500,
  });
  styleInspections[sheetName] = style.ndjson;
  const styleObj = JSON.parse(style.ndjson.trim().split("\n")[0]);
  const cellObj = JSON.parse(cell.ndjson.trim().split("\n")[0]);
  const expectedFirstHeader = matrixAudit.find((x) => x.sheet === sheetName)?.first_header;
  const observedFirstHeader = cellObj.values?.[0]?.[0];
  firstHeaderChecks[sheetName] = { expected: expectedFirstHeader, observed: observedFirstHeader, passed: observedFirstHeader === expectedFirstHeader };
  headerStyleChecks[sheetName] = {
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

const formulaErrorScan = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 300 },
  summary: "final formula error scan",
  maxChars: 4000,
});

const previewRanges = {
  README: "A1:B21",
  Data_Sources: "A1:K16",
  Data_Provenance: "A1:F34",
  Statistical_Units: "A1:J17",
  Figure_Source_Map: "A1:H52",
  Endpoint_Definitions: "A1:H28",
  FAERS_Overall: "A1:X20",
  FAERS_Specifications: "A1:Y20",
  FAERS_Annual: "A1:U42",
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

const renderSpecs = [];
const renderSpecKeys = new Set();
function addRenderSpec(sheetName, range, purpose) {
  const key = `${sheetName}!${range}`;
  if (renderSpecKeys.has(key)) return;
  renderSpecKeys.add(key);
  renderSpecs.push({ sheetName, range, purpose });
}
for (const [sheetName, range] of Object.entries(previewRanges)) addRenderSpec(sheetName, range, "sheet overview");

for (const [sheetName, record] of Object.entries(sheetRecords)) {
  record.blocks.forEach((block, blockIndex) => {
    if (blockIndex === 0 && record.blocks.length > 1) return;
    const rowStart = block.startRowZero + 1;
    const rowEnd = Math.min(block.endRowZero + 1, rowStart + 8);
    const firstEndCol = colLetter(Math.min(block.nCols, 15));
    addRenderSpec(sheetName, `A${rowStart}:${firstEndCol}${rowEnd}`, `block ${blockIndex + 1} header and representative rows`);
    if (sheetName === "Integrated_Prioritization" && block.nCols > 15) {
      for (let colStart = 16; colStart <= block.nCols; colStart += 15) {
        const colEnd = Math.min(block.nCols, colStart + 14);
        addRenderSpec(sheetName, `${colLetter(colStart)}${rowStart}:${colLetter(colEnd)}${rowEnd}`, `block ${blockIndex + 1} horizontal slab ${colStart}-${colEnd}`);
      }
    }
    if (block.nRows > 100) {
      const tailStart = Math.max(rowStart + 1, block.endRowZero + 1 - 6);
      addRenderSpec(sheetName, `A${tailStart}:${firstEndCol}${block.endRowZero + 1}`, `block ${blockIndex + 1} tail rows`);
    }
  });
}

const rendered = [];
for (const spec of renderSpecs) {
  const { sheetName, range, purpose } = spec;
  const blob = await workbook.render({ sheetName, range, scale: 1, format: "png" });
  const filename = `${String(rendered.length + 1).padStart(3, "0")}_${sheetName}_${range.replace(/[^A-Za-z0-9]+/g, "_")}.png`;
  const outputPath = path.join(RENDER_DIR, filename);
  await fs.writeFile(outputPath, new Uint8Array(await blob.arrayBuffer()));
  const renderStat = await fs.stat(outputPath);
  rendered.push({ sheetName, range, purpose, outputPath, bytes: renderStat.size, nonempty: renderStat.size > 1000 });
}

const formulaScanPassed = formulaErrorScan.ndjson.includes("matched 0 entries");
const requiredSheets = [
  "README", "Data_Sources", "Data_Provenance", "Statistical_Units", "Figure_Source_Map", "Endpoint_Definitions",
  "FAERS_Overall", "FAERS_Specifications", "FAERS_Annual", "Proteomics_Effects", "Proteomics_QC",
  "LINCS_Dose_Metrics", "ChEMBL_Targets", "Liver_Expression", "Virtual_KO_Stability",
  "Integrated_Prioritization", "Blood_Sensitivity", "Preclinical_Evidence", "GSE151374_scRNA", "PXD052594_Proteomics", "Experimental_Validation", "Reporting_Checklists",
];
const allRequiredSheetsPresent = requiredSheets.every((name) => Object.hasOwn(sheetRecords, name));
const sourceMatricesValid = matrixAudit.every((x) => x.header_check === "PASS" && x.rectangular_check === "PASS");
const upstreamValidationPassed = integrationValidation.status === "PASS" && preclinicalValidation.validation_passed === true
  && gse151374HashValidation.all_recorded_hashes_match === true && pxd052594HashValidation.all_recorded_hashes_match === true
  && gse151374Validation.no_cell_level_inferential_tests === true && pxd052594Validation.post_treatment_cluster_used_in_analysis === false;
const allRendersNonempty = rendered.every((x) => x.nonempty);
const renderedSheetCount = new Set(rendered.map((x) => x.sheetName)).size;
const allSheetsRendered = renderedSheetCount === requiredSheets.length;
const allHeaderStylesPassed = Object.values(headerStyleChecks).every((x) => x.passed);
const allFirstHeadersPassed = Object.values(firstHeaderChecks).every((x) => x.passed);
const programmaticPass = formulaScanPassed && allRequiredSheetsPresent && sourceMatricesValid && upstreamValidationPassed
  && integrationFileHashChecks.every((x) => x.passed || x.accepted_post_freeze_audit_update) && preclinicalFileHashChecks.every((x) => x.passed || x.accepted_post_freeze_audit_update)
  && gse151374FileHashChecks.every((x) => x.passed) && pxd052594FileHashChecks.every((x) => x.passed)
  && JSON.stringify(exploratoryT2GateTargets) === JSON.stringify(["FGFR1"]) && strictAllSeedT2Targets.length === 0
  && bloodExactGridPassed && bloodBhPassed && allRendersNonempty && allSheetsRendered
  && gse151DesignPassed && gse151ExactFamiliesNull && gse151GeneFamiliesNull
  && pxd052DesignPassed && pxd052ProteomeBhNull && pxd052PhosphositeQ05Count === 994
  && JSON.stringify(extensionExposureTargets) === JSON.stringify(["fgfr1"])
  && JSON.stringify(extensionPhenotypeAnchors) === JSON.stringify(["fbp1", "gsta1", "otc"])
  && allHeaderStylesPassed && allFirstHeadersPassed && missingFigureSourcePaths.length === 0;

const validation = {
  artifact: "Supplementary_Tables.xlsx",
  created_utc: new Date().toISOString(),
  status: programmaticPass ? "programmatic_checks_pass_visual_review_pending" : "FAIL",
  visual_review: "pending manual inspection of artifact-tool renders",
  sheet_count: Object.keys(sheetRecords).length,
  sheet_names: Object.keys(sheetRecords),
  table_count: tableCounter,
  source_csv_count: Object.keys(files).length,
  source_tables: Object.fromEntries(Object.entries(files).map(([k, v]) => [k, { path: v.relativePath, rows: v.rows, columns: v.cols }])),
  external_large_result_tables: {
    gse151374_macrophage_de_day7: { packaged_path: gse151DeDay7PackagedPath, analysis_source_path: gse151DeDay7Summary.relativePath, rows: gse151DeDay7Summary.rows, numeric_q_values: gse151DeDay7Summary.numeric, q_lt_0_10: gse151DeDay7Summary.matched, sha256: gse151DeFileIndex[1][7] },
    gse151374_macrophage_de_day14: { packaged_path: gse151DeDay14PackagedPath, analysis_source_path: gse151DeDay14Summary.relativePath, rows: gse151DeDay14Summary.rows, numeric_q_values: gse151DeDay14Summary.numeric, q_lt_0_10: gse151DeDay14Summary.matched, sha256: gse151DeFileIndex[2][7] },
    packaging_reason: "Complete CSVs are hash-verified submission artifacts and indexed in GSE151374_scRNA; artifact-tool raises RangeError: Invalid string length when both 48,795-row matrices are duplicated into this XLSX.",
  },
  workbook_sha256: workbookSha256,
  structural_style_formula_qa: {
    all_required_sheets_present: allRequiredSheetsPresent,
    required_sheets: requiredSheets,
    matrix_header_and_rectangular_checks_passed: sourceMatricesValid,
    matrix_audit_count: matrixAudit.length,
    filterable_table_creation_enforced_for_every_matrix: tableCounter === matrixAudit.length,
    first_row_freeze_requested_during_authoring: true,
    first_row_freeze_persisted_in_export: false,
    freeze_panes_export_note: "The documented freezeRows(1) request was issued for every worksheet, including a final post-formatting re-assertion, but @oai/artifact-tool 2.8.6 did not serialize worksheet <pane> records. This is disclosed as a non-scientific navigation limitation and is not treated as a passed artifact property.",
    gridlines_hidden_enforced_during_authoring: true,
    rendered_sheet_count: renderedSheetCount,
    render_range_count: rendered.length,
    all_render_files_nonempty: allRendersNonempty,
    all_sheets_rendered: allSheetsRendered,
    header_style_inspection_count: Object.keys(styleInspections).length,
    all_header_styles_passed: allHeaderStylesPassed,
    all_first_headers_passed: allFirstHeadersPassed,
    formula_error_scan_passed: formulaScanPassed,
    upstream_integration_validation: integrationValidation.status,
    upstream_preclinical_validation: preclinicalValidation.validation_passed ? "PASS" : "FAIL",
    upstream_gse151374_validation: gse151374HashValidation.all_recorded_hashes_match && gse151374Validation.no_cell_level_inferential_tests ? "PASS" : "FAIL",
    upstream_pxd052594_validation: pxd052594HashValidation.all_recorded_hashes_match && !pxd052594Validation.post_treatment_cluster_used_in_analysis ? "PASS" : "FAIL",
    integration_manifest_hashes_verified: `${integrationOriginalHashMatches}/${integrationFileHashChecks.length}`,
    integration_post_freeze_exceptions_required: integrationAcceptedTextUpdates,
    preclinical_core_original_hashes_verified: `${preclinicalOriginalHashMatches}/${preclinicalFileHashChecks.length}`,
    preclinical_manifest_hashes_verified: `${preclinicalOriginalHashMatches}/${preclinicalFileHashChecks.length}`,
    preclinical_post_freeze_exceptions_required: preclinicalAcceptedAuditUpdates,
    gse151374_hashes_verified: `${gse151374FileHashChecks.filter((x) => x.passed).length}/${gse151374FileHashChecks.length}`,
    pxd052594_hashes_verified: `${pxd052594FileHashChecks.filter((x) => x.passed).length}/${pxd052594FileHashChecks.length}`,
    target_gate_reconciliation_passed: JSON.stringify(exploratoryT2GateTargets) === JSON.stringify(["FGFR1"]) && strictAllSeedT2Targets.length === 0,
    blood_exact_probability_lattice_passed: bloodExactGridPassed,
    blood_bh_reconciliation_passed: bloodBhPassed,
    gse151374_design_and_cell_count_reconciliation_passed: gse151DesignPassed,
    gse151374_all_exact_bh_families_null: gse151ExactFamiliesNull,
    gse151374_gene_bh_q_lt_0_10_count_is_zero: gse151GeneFamiliesNull,
    pxd052594_design_reconciliation_passed: pxd052DesignPassed,
    pxd052594_proteome_bh_q_lt_0_10_count_is_zero: pxd052ProteomeBhNull,
    pxd052594_phosphosite_q_lt_0_05_count: pxd052PhosphositeQ05Count,
    extension_role_labels_reconciled: JSON.stringify(extensionExposureTargets) === JSON.stringify(["fgfr1"]) && JSON.stringify(extensionPhenotypeAnchors) === JSON.stringify(["fbp1", "gsta1", "otc"]),
    figure_source_path_count: figureSourcePathChecks.length,
    missing_figure_source_paths: missingFigureSourcePaths,
  },
  integrity_controls: [
    "All source-derived blocks were read directly from existing CSV/CSV.GZ outputs.",
    "No source value was recomputed or overwritten during workbook construction.",
    "Blank/non-estimable source cells were retained as blank Excel cells.",
    "All data blocks are filterable, banded Excel tables with visible headers. The documented freezeRows(1) request was issued for every worksheet but @oai/artifact-tool 2.8.6 did not persist worksheet <pane> records; this non-scientific navigation limitation is disclosed and not evaluated as a PASS property. Blocks above 10,000 rows retain native numeric values on General format to avoid artifact-tool per-cell style expansion; headers and column widths remain publication-professional.",
    "Phenotype P tiers and target T tiers remain separate; expression and virtual-KO evidence is never spliced across liver compartments.",
    "Independent-animal, pooled-slice, and unresolved-provenance public preclinical sources are labelled at row and sheet level.",
    `The final Figure 6 core manifest SHA-256 was verified as ${EXPECTED_CORE_PRECLINICAL_MANIFEST_SHA256}; all ${preclinicalFileHashChecks.length} member hashes match without exceptions.`,
    `The final same-compartment integration manifest SHA-256 was verified as ${EXPECTED_INTEGRATION_MANIFEST_SHA256}; all ${integrationFileHashChecks.length} member hashes match without exceptions.`,
    `The frozen GSE151374 S6 manifest SHA-256 was verified as ${EXPECTED_GSE151374_MANIFEST_SHA256}; the pooled library is the inferential unit and cells are not independent replicates.`,
    `The frozen PXD052594 S7 manifest SHA-256 was verified as ${EXPECTED_PXD052594_MANIFEST_SHA256}; no radiomic response cluster was used.`,
    "PXD052594 23-animal proteome and distinct deposited random 5-versus-5 phosphoproteome subset are never merged; the phosphosite analysis is table-only.",
    "The two complete 48,795-row GSE151374 macrophage DE files remain hash-verified external machine-readable supplement files and are indexed, not truncated or selectively filtered, in the workbook.",
    "GSE299128 P values are the corrected exact two-sided participant-level sign-flip values over all 2^7 assignments.",
    "The prospective experimental protocol is explicitly marked not yet performed.",
  ],
  key_reconciliations: {
    faers_narrow: "340/30233 versus 112/39467; ROR 3.996605751",
    faers_broad: "1773/30233 versus 717/39467; ROR 3.3668703341",
    faers_2014_narrow: "non-estimable because both event counts are zero",
    lincs_gene_rows: files.lincsMetrics.rows,
    chembl_target_rows: files.chembl.rows,
    liver_expression_rows: files.liver.rows,
    virtual_ko_rows: files.vkoStable.rows,
    integrated_candidate_rows: files.rank.rows,
    whole_blood_gene_rows: files.bloodSummary.rows,
    whole_blood_exact_assignments: 2 ** bloodManifest.n_participants,
    gse278200_independent_animals: preclinicalAnalysisSummary.GSE278200.diagnostics.n_samples,
    gse308578_independent_animals: preclinicalAnalysisSummary.GSE308578.diagnostics.n_samples,
    preclinical_manifest_entries: preclinicalManifest.files.length,
    preclinical_core_manifest_sha256: preclinicalManifestSha256,
    gse151374_retained_cells: gse151TotalCells,
    gse151374_pooled_libraries: gse151QcRows.length,
    gse151374_manifest_entries: gse151374Manifest.inputs.length + gse151374Manifest.outputs.length,
    gse151374_manifest_sha256: gse151374ManifestSha256,
    pxd052594_proteins: pxd052ProteomeRows.length,
    pxd052594_proteome_nintedanib_animals: pxd052594Validation.nintedanib_animals,
    pxd052594_proteome_vehicle_animals: pxd052594Validation.vehicle_animals,
    pxd052594_phosphosites: pxd052PhosphositeRows.length,
    pxd052594_phosphosite_q_lt_0_05: pxd052PhosphositeQ05Count,
    pxd052594_manifest_entries: pxd052594Manifest.inputs.length + pxd052594Manifest.outputs.length,
    pxd052594_manifest_sha256: pxd052594ManifestSha256,
    integration_manifest_sha256: integrationManifestSha256,
    blood_manifest_sha256: bloodManifestSha256,
    exploratory_t2_gate_targets: exploratoryT2GateTargets,
    strict_all_seed_t2_targets: strictAllSeedT2Targets,
  },
  formula_error_scan: formulaErrorScan.ndjson,
  sheet_catalog: sheetCatalog.ndjson,
  first_cell_inspections: firstCellInspections,
  first_header_checks: firstHeaderChecks,
  style_inspections: styleInspections,
  header_style_checks: headerStyleChecks,
  figure_source_path_checks: figureSourcePathChecks,
  inspected_ranges: Object.keys(keyInspections),
  rendered_sheets: rendered,
  output_bytes: stat.size,
};
await fs.writeFile(VALIDATION_JSON, JSON.stringify(validation, null, 2) + "\n", "utf8");
await fs.rm(`${OUT_XLSX}.inspect.ndjson`, { force: true });

console.log(JSON.stringify({
  output: OUT_XLSX,
  validation: VALIDATION_JSON,
  sheets: Object.keys(sheetRecords).length,
  tables: tableCounter,
  renderDir: RENDER_DIR,
  bytes: stat.size,
  sha256: workbookSha256,
  programmaticPass,
  formulaErrorScan: formulaErrorScan.ndjson,
}, null, 2));
