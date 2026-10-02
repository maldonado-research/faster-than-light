#!/usr/bin/env python3
"""Verify strengthened registered checks reject in-memory helper sign/formula defects."""
import argparse
import sys
sys.dont_write_bytecode = True
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--candidate',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args()
candidate_path=args.candidate.resolve()
source_hash=hashlib.sha256(candidate_path.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('candidate_gate_review',candidate_path)
candidate=importlib.util.module_from_spec(spec); spec.loader.exec_module(candidate)
baseline=candidate.run_checks()
checks=[]

def substitute(name,make_defect,description):
    original=getattr(candidate,name)
    try:
        setattr(candidate,name,make_defect(original))
        result=candidate.run_checks()
    finally:
        setattr(candidate,name,original)
    checks.append({'defect':description,'rejected':not result['passed'],
                   'failure_count':len(result['failures']),'failure_examples':result['failures'][:5]})

substitute('frequency_squared',lambda original:lambda *parameters:abs(original(*parameters)),
           'Erase negative dispersion values; healthy frequencies unchanged')
substitute('free_energy_density',lambda original:lambda *parameters:abs(original(*parameters)),
           'Erase negative Hamiltonian density; healthy energies unchanged')
substitute('principal_coefficients',lambda original:lambda *parameters:(abs(original(*parameters)[0]),*original(*parameters)[1:]),
           'Erase negative boosted time coefficient')
substitute('omega',lambda original:lambda *parameters:original(*parameters)*(1+1e-8),
           'Rescale mode frequency by 1+1e-8')
result={'run_utc':datetime.now(timezone.utc).isoformat(),'candidate_sha256':source_hash,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'method':'In-memory substitutions only; tracked source remains unchanged',
        'baseline_passed':baseline['passed'],'registered_result_hash_matches':
             json.loads((candidate_path.parent/'RESULTS.json').read_text())['script_sha256']==source_hash,
        'candidate_unchanged_during_run':hashlib.sha256(candidate_path.read_bytes()).hexdigest()==source_hash,
        'defects':checks,'passed':baseline['passed'] and all(check['rejected'] for check in checks)}
args.out.parent.mkdir(parents=True,exist_ok=True)
args.out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
raise SystemExit(0 if result['passed'] and result['candidate_unchanged_during_run'] else 1)
