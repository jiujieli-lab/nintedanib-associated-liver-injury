#!/usr/bin/env python3
"""Read-only integrity and static source checks. Does not execute analysis code."""
from pathlib import Path
import ast
import csv
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent

def main():
    errors = []
    checked = 0
    listed = set()
    for line in (ROOT / 'SHA256SUMS.txt').read_text(encoding='utf8').splitlines():
        if not line.strip():
            continue
        digest, rel = line.split('  ', 1)
        path = (ROOT / rel).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            errors.append('Missing or unsafe path: ' + rel)
            continue
        listed.add(rel)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            errors.append('Hash mismatch: ' + rel)
        checked += 1
    inventory = list(csv.DictReader((ROOT / 'documentation/Code_Inventory.csv').open(encoding='utf8', newline='')))
    code_hashes = set()
    syntax_count = 0
    for row in inventory:
        path = ROOT / row['package_path']
        if not path.is_file():
            errors.append('Missing original script: ' + row['package_path'])
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != row['sha256']:
            errors.append('Original source changed: ' + row['package_path'])
        code_hashes.add(digest)
        if path.suffix == '.py':
            try:
                ast.parse(path.read_text(encoding='utf8'), filename=row['package_path'])
                syntax_count += 1
            except SyntaxError as exc:
                errors.append('Python syntax: ' + row['package_path'] + ': ' + str(exc))
    summary = json.loads((ROOT / 'documentation/Consolidation_Summary.json').read_text(encoding='utf8'))
    if len(code_hashes) != summary['unique_original_source_code_files']:
        errors.append('Original code count differs from consolidation summary.')
    source_map = list(csv.DictReader((ROOT / 'documentation/Source_File_Map.csv').open(encoding='utf8', newline='')))
    for row in source_map:
        if row['package_path']:
            path = ROOT / row['package_path']
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
                errors.append('Source mapping mismatch: ' + row['source_archive'] + '!/' + row['source_path'])
    counts = {
        'jcm_revision/screening/Boltz_All32_Dispositions.csv': 32,
        'jcm_revision/screening/Inductive_Properties_All16.csv': 48,
        'jcm_revision/docking/Docking_All_Completed_Runs.csv': 96,
        'jcm_revision/docking/Docking_Summary.csv': 32,
        'jcm_revision/candidate_library/All_Target_Activity_Disposition.csv': 191,
        'jcm_revision/candidate_library/Davis2011_Matched_Kd_Benchmark.csv': 12,
    }
    observed = {}
    for rel, expected in counts.items():
        path = ROOT / rel
        with path.open(encoding='utf8', newline='') as handle:
            observed[rel] = sum(1 for _ in csv.DictReader(handle))
        if observed[rel] != expected:
            errors.append('Archived row count differs: ' + rel)
    # Extra runtime files are not included in the frozen integrity claim.
    unlisted = sorted(p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
                      if p.is_file() and p.name != 'SHA256SUMS.txt'
                      and p.relative_to(ROOT).as_posix() not in listed)
    report = {'status': 'PASS' if not errors else 'FAIL', 'hash_checked_files': checked,
              'unique_original_scripts': len(code_hashes), 'python_syntax_pass': syntax_count,
              'source_member_records': len(source_map), 'archived_counts': observed,
              'unlisted_extra_files': unlisted, 'errors': errors,
              'analyses_executed': False, 'network_calls': False}
    print(json.dumps(report, indent=2))
    return 1 if errors else 0

if __name__ == '__main__':
    sys.exit(main())
