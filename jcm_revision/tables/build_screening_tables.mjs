import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';
const root=process.argv[2];
if(!root)throw Error('Supply revision root');
const out=path.join(root,'tables/output'); await fs.mkdir(path.join(out,'renders'),{recursive:true});
const defs=[
 ['N20_Library','Compound identity and clinical context','candidate_library/Candidate_Library_Evidence.csv','Verified parent structures, approved indications, mechanisms and hepatic warnings.'],
 ['N21_Boltz','AI screening of both receptors','screening/Boltz_All32_Dispositions.csv','All 32 input dispositions; filtered compounds have no score. Scores are target-specific model outputs.'],
 ['N22_Docking_Runs','Individual docking searches','docking/Docking_All_Completed_Runs.csv','All 96 primary searches; three computational seeds per compound and target.'],
 ['N23_Docking_Summary','Docking score and pose consistency','docking/Docking_Summary.csv','Within-target median and observed range; seeds are not biological replicates.'],
 ['N24_Properties','Inductive molecular properties','screening/Inductive_Properties_All16.csv','All 48 predictions with versions, structure identity and provider uncertainty flags.'],
 ['N25_Assay_Disposition','Complete experimental activity records','candidate_library/All_Target_Activity_Disposition.csv','All 191 retrieved records, including censored, variant, duplicate and ineligible observations.'],
 ['N26_Activity_Summary','Endpoint-specific experimental summary','candidate_library/Endpoint_Separated_Experimental_Summary.csv','Kd, Ki and IC50 remain separate; no cross-assay inferential pooling.'],
 ['N27_Matched_Kd','Source-matched kinase-binding panel','candidate_library/Davis2011_Matched_Kd_Benchmark.csv','Reported inequalities and potential-duplicate flags are retained; not independent replicate measurements.']
];
const wb=Workbook.create();const index=wb.worksheets.add('Table_Index');
const col=i=>{let s='';for(let n=i+1;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
const typed=(v,h)=>{if(v===null||v==='')return null;const s=String(v);if(s.startsWith('='))return "'"+s;if(/id$|identifier|version|year|seed|smiles|key|path|sha256/.test(h))return s;if(/^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/.test(s)&&Number.isFinite(Number(s)))return Number(s);return s;};
const manifest=[];const qa=[];
function setup(sh,title,source,headers,data,tabname){
 const last=col(headers.length-1),end=data.length+5;
 sh.showGridLines=false;sh.getRange(`A1:${last}${end}`).format.font={name:'Arial',size:10,color:'#17212B'};
 sh.getRange(`A1:${last}${end}`).format.verticalAlignment='center';
 sh.getRange('A2').values=[[title]];sh.getRange('A2').format.font={name:'Arial',size:14,bold:true,color:'#17324D'};
 sh.getRange('A3').values=[[source]];sh.getRange('A3').format.font={name:'Arial',size:10,italic:true,color:'#536171'};
 sh.getRange(`A5:${last}5`).values=[headers];if(data.length)sh.getRange(`A6:${last}${end}`).values=data;
 const t=sh.tables.add(`A5:${last}${end}`,true,tabname);t.style='TableStyleLight1';t.showFilterButton=true;
 sh.getRange(`A5:${last}5`).format={fill:'#244662',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',rowHeight:44};
 if(data.length){sh.getRange(`A6:${last}${end}`).format.rowHeight=60;sh.getRange(`A6:${last}${end}`).format.wrapText=true;}
 const widths=[];
 for(let j=0;j<headers.length;j++){
  const h=headers[j];const range=sh.getRange(`${col(j)}5:${col(j)}${end}`);
  range.format.columnWidth=Math.min(55,Math.max(18,Math.ceil(h.length*.82)+2));
  if(/name|compound/.test(h))range.format.columnWidth=24;
  if(/role|mechanism|status|context|resolution/.test(h))range.format.columnWidth=42;
  if(/smiles|key|config_version/.test(h))range.format.columnWidth=68;
  if(/smiles|description|warning|source|indication|note|url|reason|provenance|identity_resolution|safety/.test(h))range.format.columnWidth=60;
  widths.push(range.format.columnWidth);
  if(data.length)sh.getRange(`${col(j)}6:${col(j)}${end}`).format.horizontalAlignment='left';
  const nums=data.map(r=>r[j]).filter(v=>typeof v==='number');
  if(nums.length){const format=/^n_|count|rank_|records|searches|duplicate$/i.test(h)?'#,##0':nums.some(x=>x!==0&&Math.abs(x)<.0001)?'0.000E+00':'0.0000';sh.getRange(`${col(j)}6:${col(j)}${end}`).setNumberFormat(format);sh.getRange(`${col(j)}6:${col(j)}${end}`).format.horizontalAlignment='center';}
 }
 for(let i=0;i<data.length;i++){const lines=Math.max(...data[i].map((v,j)=>typeof v==='string'?Math.ceil(v.length/(widths[j]*.78)):1));sh.getRange(`A${i+6}:${last}${i+6}`).format.rowHeight=Math.min(400,Math.max(32,lines*13+15));}
 sh.freezePanes.freezeRows(5);sh.freezePanes.freezeColumns(Math.min(2,headers.length));
}
for(let i=0;i<defs.length;i++){
 const [name,title,rel,note]=defs[i];const csv=await fs.readFile(path.join(root,rel),'utf8');
 const safe=csv.replace(/(^|,)(=)(?=,|\r?\n)/gm,"$1'=");
 const tmp=await Workbook.fromCSV(safe,{sheetName:'CSV'});const vals=tmp.worksheets.getItemAt(0).getUsedRange().values;
 const originalHeaders=vals[0].map(String);
 const priority={
  N20_Library:['compound','candidate_role','regulatory_status','target_mechanism','hepatic_safety_summary','hepatic_warning_category'],
  N22_Docking_Runs:['name','target','seed','vina_score_kcal_mol','elapsed_seconds','role'],
  N23_Docking_Summary:['name','target','completed_searches','median_vina_score_kcal_mol','min_vina_score_kcal_mol','max_vina_score_kcal_mol','max_pairwise_top_pose_RMSD_A','receptor_type'],
  N24_Properties:['compound','model_id','value','model_lower_bound','model_upper_bound','out_of_domain_flag','low_confidence_flag','bounds_coverage','model_version','config_version'],
  N25_Assay_Disposition:['compound','gene_symbol','standard_type','standard_relation','standard_value','standard_units','eligible_exact_human_mapped_nonvariant','potential_duplicate','disposition'],
  N27_Matched_Kd:['compound','gene_symbol','standard_relation','standard_value','standard_units','potential_duplicate','document_doi']
 }[name]||[];
 const headers=[...priority,...originalHeaders.filter(h=>!priority.includes(h))];
 const data=vals.slice(1).filter(r=>r.some(x=>x!==null&&x!=='')).map(r=>headers.map(h=>typed(r[originalHeaders.indexOf(h)],h)));
 if(name==='N22_Docking_Runs'&&data.length!==96)throw Error('Primary docking incomplete');
 if(name==='N21_Boltz'&&data.length!==32)throw Error('Boltz dispositions incomplete');
 const sh=wb.worksheets.add(name);setup(sh,`Supplementary Table N${20+i}. ${title}`,`Source: ${path.basename(rel)}`,headers,data,`Table_N${20+i}`);
 const actual=sh.getRange(`A6:${col(headers.length-1)}${data.length+5}`).values;
 const expected=data.map(r=>r.map(v=>typeof v==='string'&&v.startsWith("'=")?v.slice(1):v));
 if(JSON.stringify(actual)!==JSON.stringify(expected))throw Error(`Source-value mismatch: ${name}`);
 manifest.push({table:`N${20+i}`,sheet:name,title,file:rel,rows:data.length,columns:headers.length,note,sha256:crypto.createHash('sha256').update(csv).digest('hex')});
 qa.push({sheet:name,exact_source_values:true,rows:data.length,columns:headers.length});
}
const ih=['Table','Content','Worksheet','Records','Interpretation'];
const id=manifest.map(m=>[m.table,m.title,m.sheet,m.rows,m.note]);
setup(index,'Receptor screening and molecular-property tables','Supplementary Tables N20–N27; source files retain the complete precision.',ih,id,'Screening_Table_Index');
index.getRange('B5:B13').format.columnWidth=48;index.getRange('E5:E13').format.columnWidth=80;
index.getRange('A6:E13').format.wrapText=true;index.getRange('A6:E13').format.rowHeight=60;index.getRange('D6:D13').setNumberFormat('0');index.tabColor='#244662';
wb.recalculate();
for(const name of ['Table_Index',...defs.map(x=>x[0])]){
 const range=name==='Table_Index'?'A1:E13':`A1:${col(Math.min(5,manifest.find(m=>m.sheet===name).columns-1))}12`;
 const blob=await wb.render({sheetName:name,range,scale:1.2,format:'png'});await fs.writeFile(path.join(out,'renders',name+'.png'),new Uint8Array(await blob.arrayBuffer()));
}
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:50},maxChars:3000});
const f=await SpreadsheetFile.exportXlsx(wb);await f.save(path.join(out,'Additional_file_8_Screening_Tables.xlsx'));
await fs.writeFile(path.join(out,'Table_Manifest.json'),JSON.stringify(manifest,null,2));
await fs.writeFile(path.join(out,'Workbook_QA.json'),JSON.stringify({sheets:9,checks:qa,formula_errors:errors.ndjson},null,2));
console.log(JSON.stringify({file:path.join(out,'Additional_file_8_Screening_Tables.xlsx'),sheets:9,rows:manifest.map(m=>[m.table,m.rows])}));
