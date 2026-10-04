import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const root=process.argv[2] || '/workspace/scratch/c5634a2d588f';
const out=path.join(root,'work/tables');
await fs.mkdir(path.join(out,'renders'),{recursive:true});
const defs=[
 ['N1_Binding_Summary','Binding assay-stratified summary','binding_analysis/results/Binding_Assay_Stratified_Summary.csv','ChEMBL biochemical activity records; IC50, Kd and Ki retain their endpoint labels.'],
 ['N2_Document_Contrasts','Within-document target contrasts','binding_analysis/results/Binding_Within_Document_Contrasts.csv','Matched document and endpoint summaries; differences are descriptive.'],
 ['N3_Binding_LODO','Binding leave-one-document-out analysis','binding_analysis/results/Binding_Leave_One_Document_Out.csv','Independent-document sensitivity; omitted document recorded for every estimate.'],
 ['N4_Physicochemical','Inductive Bio physicochemical predictions','binding_analysis/results/Inductive_Physicochemical_Predictions.csv','Connector predictions for exact nintedanib parent SMILES; returned bounds have unspecified coverage.'],
 ['N5_Binding_Audit','Binding record audit counts','binding_analysis/results/Binding_Record_Audit_Counts.csv','All target activity records assigned to explicit audit categories.'],
 ['N6_Binding_Records','Complete audited binding records','binding_analysis/results/Binding_All_Target_Records_Audited.csv','All 26 PDGFRA/FLT4 activity records retained, including excluded and non-quantitative records.'],
 ['N7_Primary_Calibration','Virtual intervention primary calibration','new_virtual/primary_calibration.csv','Nine tests; BH q values retained for all targets and compartments.'],
 ['N8_Intervention_Runs','All virtual intervention runs','new_virtual/perturbation_metrics.csv','Full baseline, sham, intervention, restoration and donor-omission computational outputs.'],
 ['N9_Nonadditivity','Paired intervention nonadditivity calibration','new_virtual/nonadditivity_calibration.csv','Three compartment tests. Network nonadditivity does not establish biological synergy.'],
 ['N10_Nonadditivity_LODO','Nonadditivity under donor omission','new_virtual/nonadditivity_by_donor_omission.csv','Donor-omission sensitivity of the same prespecified joint-intervention contrast.'],
 ['N11_Target_Expression','Donor-level target expression','new_virtual/donor_target_expression.csv','Measured atlas transcript detection and abundance, stratified by donor and compartment.'],
 ['N12_Attenuation','Incident-edge attenuation sensitivity','new_virtual/edge_attenuation_sensitivity.csv','Graded network-edge attenuation. This fraction is not a drug dose or receptor-inhibition percentage.'],
 ['N13_Matched_Controls','Matched control gene pools','new_virtual/matched_control_pools.csv','Expression, detection and graph-degree matched genes used for empirical null calibration.'],
 ['N14_Target_Readouts','Target-to-target network redistribution','new_virtual/target_to_target_redistribution.csv','Model score redistribution after graph deletion; not a directed receptor-regulation estimate.'],
 ['N15_Covariance_QC','Donor covariance diagnostics','new_virtual/donor_covariance_qc.csv','Sample counts, nonconstant feature counts and shrinkage coefficients.'],
 ['N16_Omission_Ranges','Donor-omission effect ranges','new_virtual/donor_omission_ranges.csv','Minimum and maximum across donor omissions; these ranges are not confidence intervals.'],
 ['N17_Input_Provenance','Virtual intervention input provenance','new_virtual/input_manifest.csv','Input matrix filenames, SHA-256 hashes, cells, features and biological donors.'],
];
const wb=Workbook.create();
const index=wb.worksheets.add('Table_Index');
const manifest=[];
const qa={expectedSheets:18,checks:[],dataCellsChecked:0,renderedSheets:[],formulaErrors:null};
const col=i=>{let s='';for(let n=i+1;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
function typed(v,h){
 if(v===null||v==='') return null;
 const s=String(v);
 if(s.startsWith('='))return "'"+s;
 if(s==='True'||s==='False')return s;
 if(/^(activity_id|record_id|src_id|toid|target_tax_id|document_year)$/.test(h))return s;
 if(/^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/.test(s)&&Number.isFinite(Number(s)))return Number(s);
 return s;
}
function preferredColumns(file,headers){
 let first=[];
 if(file.includes('All_Target_Records'))first=['gene_symbol','activity_id','document_chembl_id','document_year','standard_type','standard_value','standard_units','pchembl_value','audit_category','in_frozen_analysis','included_format_sensitivity'];
 if(file.includes('Physicochemical'))first=['compound','model_id','predicted_value','returned_lower_bound','returned_upper_bound','units','out_of_domain_flag','low_confidence_flag','interval_coverage'];
 if(file.includes('primary_calibration'))first=['compartment','intervention','n_donors','n_cells','clinical_weighted_relative_change','empirical_p','q_bh_9','omitted_donor'];
 if(file==='nonadditivity_calibration.csv')first=['compartment','clinical_nonadditivity_relative_L1','nonadditivity_fraction_of_pair_change','empirical_p','q_bh_3','n_null'];
 return [...first,...headers.filter(h=>!first.includes(h))];
}
function setup(sheet,title,source,headers,data,tableName){
 sheet.showGridLines=false;
 const last=col(headers.length-1),end=data.length+5;
 sheet.getRange(`A1:${last}${end}`).format.font={name:'Arial',size:10,color:'#17212B'};
 sheet.getRange(`A1:${last}${end}`).format.verticalAlignment='center';
 sheet.getRange('A2').values=[[title]];sheet.getRange('A2').format.font={name:'Arial',size:14,bold:true,color:'#17324D'};
 sheet.getRange('A3').values=[[source]];sheet.getRange('A3').format.font={name:'Arial',size:10,italic:true,color:'#536171'};
 sheet.getRange(`A5:${last}5`).values=[headers];
 if(data.length)sheet.getRange(`A6:${last}${end}`).values=data;
 const table=sheet.tables.add(`A5:${last}${end}`,true,tableName);table.showFilterButton=true;table.style='TableStyleLight1';
 sheet.getRange(`A5:${last}5`).format={fill:'#244662',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:false,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:30};
 if(data.length)sheet.getRange(`A6:${last}${end}`).format.rowHeight=23;
 headers.forEach((h,j)=>{
  let width=Math.max(16,h.length+3,...data.slice(0,100).map(r=>typeof r[j]==='string'?Math.min(r[j].length*.88+2,50):14));
  if(h==='compound')width=36;
  if(h==='file'||h==='sha256'||h==='canonical_smiles'||h==='exact_smiles'||h==='assay_description')width=65;
  sheet.getRange(`${col(j)}5:${col(j)}${end}`).format.columnWidth=width;
  if(['file','assay_description','canonical_smiles','exact_smiles','ligand_efficiency','interval_coverage','p_value_note'].includes(h)){
   sheet.getRange(`${col(j)}6:${col(j)}${end}`).format.wrapText=true;
   data.forEach((r,i)=>{if(typeof r[j]==='string'&&r[j].length>width)sheet.getRange(`A${i+6}:${last}${i+6}`).format.rowHeight=Math.max(38,Math.ceil(r[j].length/(width*.96))*14);});
  }
  const nums=data.map(r=>r[j]).filter(v=>typeof v==='number');
  if(nums.length){
   const isCount=/^n_|documents|^remaining_documents$|^Additional assay format$|^Duplicate flag$|^Frozen eligible$|^Kinetics without value$|^Single-concentration inhibition$/.test(h);
   const tiny=nums.some(v=>v!==0&&Math.abs(v)<.0001);
   const form=isCount?'#,##0':tiny?'0.000E+00':/fraction|detection|percentile/.test(h)?'0.0000':/pchembl|_nM|predicted|bound|mean_expression/.test(h)?'0.0000':'0.000000';
   sheet.getRange(`${col(j)}6:${col(j)}${end}`).setNumberFormat(form);
   sheet.getRange(`${col(j)}6:${col(j)}${end}`).format.horizontalAlignment='right';
  }
 });
 sheet.freezePanes.freezeRows(5);sheet.freezePanes.freezeColumns(Math.min(2,headers.length));
 return {last,end};
}

for(let i=0;i<defs.length;i++){
 const [name,title,rel,note]=defs[i],filename=path.basename(rel);
 const csv=await fs.readFile(path.join(root,'work',rel),'utf8');
 // Exact-relation ChEMBL fields contain literal '='. Escape before CSV import
 // so the importer cannot mistake an assay relation for a spreadsheet formula.
 const safeCsv=csv.replace(/(^|,)(=)(?=,|\r?\n)/gm,"$1'=");
 const tmp=await Workbook.fromCSV(safeCsv,{sheetName:'CSV'});
 const vals=tmp.worksheets.getItemAt(0).getUsedRange().values;
 const rawHeaders=vals[0].map(String),headers=preferredColumns(filename,rawHeaders);
 const body=vals.slice(1).filter(row=>row.some(v=>v!==''&&v!==null));
 const data=body.map(row=>headers.map(h=>typed(row[rawHeaders.indexOf(h)],h)));
 const sh=wb.worksheets.add(name);
 setup(sh,`Supplementary Table N${i+1}. ${title}`,`Source: ${filename}`,headers,data,`Supplementary_N${i+1}`);
 const expected=data.map(r=>r.map(v=>typeof v==='string'&&v.startsWith("'=")?v.slice(1):v)),actual=sh.getRange(`A6:${col(headers.length-1)}${data.length+5}`).values;
 if(JSON.stringify(actual)!==JSON.stringify(expected))throw Error(`Data mismatch in ${name}`);
 qa.dataCellsChecked+=data.length*headers.length;
 qa.checks.push({sheet:name,rows:data.length,columns:headers.length,exactTypedCsvMatch:true});
 manifest.push({number:`N${i+1}`,sheet:name,filename,rows:data.length,columns:headers.length,sourceSha256:crypto.createHash('sha256').update(csv).digest('hex'),note,headerOrder:headers});
}

const idxHeaders=['Table','Title','Worksheet or archive location','CSV source','Records','Evidence type','Source reference'];
const idxData=manifest.map((m,i)=>[m.number,defs[i][1],m.sheet,m.filename,m.rows,defs[i][3],i<6?(i===3?'Inductive Bio connector; models and versions in N4':'ChEMBL CHEMBL502835; https://www.ebi.ac.uk/chembl/explore/compound/CHEMBL502835'):'GSE115469; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469']);
for(const [num,name,filename] of [['N18','Complete gene-level perturbation responses','gene_level_perturbations.csv.gz'],['N19','Complete matched-null distributions','matched_null_metrics.csv.gz']]){
 const zlib=await import('node:zlib');const raw=await fs.readFile(path.join(root,'work/new_virtual',filename));const txt=zlib.gunzipSync(raw).toString();const rows=txt.trimEnd().split('\n').length-1;
 idxData.push([num,name,'Reproducibility archive, compressed CSV',filename,rows,'Complete computational output retained; no statistical-significance filtering.','GSE115469; same frozen virtual-intervention specification']);
 manifest.push({number:num,filename,rows,external:true,sourceSha256:crypto.createHash('sha256').update(raw).digest('hex')});
}
setup(index,'Targeted binding and virtual-intervention tables','Exact source-derived results; predictions and network interventions are identified separately.',idxHeaders,idxData,'Table_Index');
index.getRange('A5:A24').format.columnWidth=10;
index.getRange('B5:B24').format.columnWidth=44;
index.getRange('C5:C24').format.columnWidth=40;
index.getRange('D5:D24').format.columnWidth=48;
index.getRange('E5:E24').format.columnWidth=12;
index.getRange('F5:F24').format.columnWidth=68;
index.getRange('G5:G24').format.columnWidth=82;
index.getRange('A6:G24').format.wrapText=true;
index.getRange('A6:G24').format.rowHeight=44;
index.getRange('E6:E24').setNumberFormat('#,##0');
index.getRange('A27').values=[['Notes']];index.getRange('A27').format.font={name:'Arial',size:10,bold:true};
const notes=[
 'All original Supplementary Tables remain in Additional file 3. This workbook contains supplementary Tables N1–N19.',
 'Blank fields retain unavailable source values and are not zero. Numeric results are stored as numbers without rounding the underlying values.',
 'pChEMBL is a potency transformation, not a statistical P value. IC50, Kd and Ki are different assay endpoints.',
 'Virtual intervention and null outputs were computed from measured liver-atlas expression; they are model-conditional network sensitivities.',
 'Donor omission ranges are descriptive. Network deletion, restoration and nonadditivity do not establish receptor inhibition, biological rescue or synergy.',
 'Inductive Bio outputs are physicochemical predictions, not target-affinity predictions or molecular-dynamics simulations.',
 'Table N18 and N19 remain complete compressed CSVs in the reproducibility archive; no rows are suppressed by significance.',
];
notes.forEach((t,i)=>index.getRange(`A${28+i}`).values=[[t]]);
index.getRange('A27:G34').format.font={name:'Arial',size:10,color:'#17212B'};
index.tabColor='#244662';
wb.recalculate();
for(const name of ['Table_Index',...defs.map(d=>d[0])]){
 const sh=wb.worksheets.getItem(name);
 const h=name==='Table_Index'?idxHeaders:manifest.find(m=>m.sheet===name).headerOrder;
 const range=`A1:${col(Math.min(7,h.length-1))}${name==='Table_Index'?10:Math.min(12,(manifest.find(m=>m.sheet===name)?.rows||5)+5)}`;
 const img=await wb.render({sheetName:name,range,scale:1.4,format:'png'});
 await fs.writeFile(path.join(out,'renders',`${name}.png`),new Uint8Array(await img.arrayBuffer()));
 qa.renderedSheets.push({sheet:name,range});
}
const inspected=await wb.inspect({kind:'table',range:'N7_Primary_Calibration!A5:H14',include:'values,formulas',tableMaxRows:10,tableMaxCols:8,maxChars:6000});
await fs.writeFile(path.join(out,'primary_inspect.ndjson'),inspected.ndjson);
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:50},summary:'Formula error scan',maxChars:3000});
qa.formulaErrors=errors.ndjson;
const file=await SpreadsheetFile.exportXlsx(wb);await file.save(path.join(out,'Additional_file_5_Targeted_Analysis_Tables.xlsx'));
await fs.writeFile(path.join(out,'table_manifest.json'),JSON.stringify(manifest,null,2));
await fs.writeFile(path.join(out,'workbook_qa.json'),JSON.stringify(qa,null,2));
console.log(JSON.stringify({workbook:path.join(out,'Additional_file_5_Targeted_Analysis_Tables.xlsx'),sheets:18,dataCellsChecked:qa.dataCellsChecked,renders:qa.renderedSheets.length}));
