#!/usr/bin/env python3
"""Independent geometric and coefficient checks of the registered preferred toy."""
from datetime import datetime, timezone
from fractions import Fraction
import argparse
import sys
sys.dont_write_bytecode = True
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--candidate', type=Path, required=True, help='Path to the registered reproduce.py')
parser.add_argument('--out', type=Path, required=True, help='JSON result path')
args = parser.parse_args()
CANDIDATE = args.candidate.resolve()
source_hash = hashlib.sha256(CANDIDATE.read_bytes()).hexdigest()
spec = importlib.util.spec_from_file_location('registered_preferred_candidate', CANDIDATE)
candidate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(candidate)
rng = random.Random(202610041)
failures = []

def dot(v,w): return math.fsum(x*y for x,y in zip(v,w))
def norm(v): return math.sqrt(dot(v,v))

def intersect_time(speed, separation, target_velocity):
    """First future sphere/moving-target intersection by bracketed root search."""
    size=norm(separation)
    if size == 0: return 0.0
    lower, upper = 0.0, size/(speed-norm(target_velocity))
    for _ in range(85):
        mid=(lower+upper)/2
        displacement=tuple(r+v*mid for r,v in zip(separation,target_velocity))
        if norm(displacement) > speed*mid: lower=mid
        else: upper=mid
    return (lower+upper)/2

max_relay_error=max_emitter_error=max_lab_error=max_ray_error=0.0
min_preferred_ratio=math.inf
for i in range(2400):
    a=1+10**rng.uniform(-3,1)
    v=rng.uniform(0,.999)
    time=10**rng.uniform(-3,2)
    delay=0 if i%4==0 else 10**rng.uniform(-3,2)
    gamma=1/math.sqrt(1-v*v)
    reception=time+intersect_time(a,(v*time,0,0),(v,0,0))
    departure=reception+gamma*delay*time
    reception_alice=departure+intersect_time(a,(-v*departure,0,0),(0,0,0))
    measured=reception_alice/time
    expected=(a+v)/(a-v)+gamma*delay*(1+v/a)
    max_relay_error=max(max_relay_error,abs(measured-expected)/expected)
    min_preferred_ratio=min(min_preferred_ratio,measured)
    if measured<1 or not math.isclose(measured,expected,rel_tol=1e-10,abs_tol=1e-10):
        failures.append('Preferred relay '+str(i))

    # Distinct emitter-frame rule: second future cone is chosen in Bob's frame.
    bob_time=gamma*(reception-v*v*reception)+delay*time
    alice_bob_time=bob_time+intersect_time(a,(-v*bob_time,0,0),(-v,0,0))
    emitter_ratio=alice_bob_time/(gamma*time)
    emitter_expected=a*a*(1-v*v)/((a-v)*(a-v))+delay*a/(gamma*(a-v))
    max_emitter_error=max(max_emitter_error,abs(emitter_ratio-emitter_expected)/max(1,abs(emitter_expected)))
    if not math.isclose(emitter_ratio,emitter_expected,rel_tol=1e-10,abs_tol=1e-10):
        failures.append('Emitter-rule comparison '+str(i))

def lab_intersections(a,v,length,angle):
    gamma=1/math.sqrt(1-v*v)
    rod=(length*math.cos(angle)/gamma,length*math.sin(angle),0.0)
    motion=(v,0,0)
    outward=intersect_time(a,rod,motion)
    # At reception, source and detector have translated together: relative rod unchanged.
    returning=intersect_time(a,tuple(-x for x in rod),motion)
    for delta,initial in ((outward,rod),(returning,tuple(-x for x in rod))):
        displacement=tuple(x+u*delta for x,u in zip(initial,motion))
        yield_residual=abs(norm(displacement)-a*delta)/max(1,a*delta)
        global max_ray_error
        max_ray_error=max(max_ray_error,yield_residual)
    return (outward+returning)/gamma

for i in range(1600):
    a=1+10**rng.uniform(-3,1)
    v=rng.uniform(0,.999)
    length=10**rng.uniform(-2,2)
    angle=rng.uniform(0,math.pi)
    gamma=1/math.sqrt(1-v*v)
    rod=(length*math.cos(angle)/gamma,length*math.sin(angle),0)
    measured=lab_intersections(a,v,length,angle)
    expected=2*math.sqrt(dot((v,0,0),rod)**2+(a*a-v*v)*dot(rod,rod))/(gamma*(a*a-v*v))
    max_lab_error=max(max_lab_error,abs(measured-expected)/expected)
    if measured<=0 or not math.isclose(measured,expected,rel_tol=1e-10,abs_tol=1e-10):
        failures.append('Arbitrary-angle lab '+str(i))

for a in (1.0,1.001,1.1,2.0,10.0):
    for v in (0,.1,.5,.9,.99):
        comparison=candidate.laboratory_times(a,v)
        measured=(lab_intersections(a,v,1,0),lab_intersections(a,v,1,math.pi/2))
        for x,y in zip(measured,comparison):
            max_lab_error=max(max_lab_error,abs(x-y)/y)
            if not math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-10):
                failures.append('Registered lab formula')

def matrix_principal(a,v):
    b=((Fraction(1),-v),(-v,Fraction(1)))
    p=((Fraction(1),Fraction(0)),(Fraction(0),-a*a))
    return tuple(tuple(sum(b[i][k]*p[k][l]*b[j][l] for k in range(2) for l in range(2))/(1-v*v)
                       for j in range(2)) for i in range(2))

matrix_count=0
for a in map(Fraction,('3/2','2','3')):
    for v in map(Fraction,('0','1/4','1/3','1/2','2/3','3/4','9/10')):
        p=matrix_principal(a,v); matrix_count+=1
        if candidate.principal_coefficients(a,v) != (p[0][0],p[0][1],p[1][1]):
            failures.append('Registered coefficients disagree with tensor transformation')
        if p[0][0]*p[1][1]-p[0][1]*p[1][0] != -a*a:
            failures.append('Exact tensor transform determinant')
        if (p[0][0]>0)!=(a*abs(v)<1) or (p[0][0]==0)!=(a*abs(v)==1):
            failures.append('Exact tensor transform slicing classification')
a,v,kx,kperp,m=Fraction(2),Fraction(3,4),Fraction(0),Fraction(1),Fraction(0)
p=matrix_principal(a,v)
quarter_discriminant=p[0][1]**2*kx*kx-p[0][0]*(p[1][1]*kx*kx-a*a*kperp*kperp-m*m)
if quarter_discriminant != Fraction(-80,7): failures.append('Computed transverse discriminant')

result={'run_utc':datetime.now(timezone.utc).isoformat(),'candidate_source_sha256':source_hash,
        'candidate_unchanged_during_run':hashlib.sha256(CANDIDATE.read_bytes()).hexdigest()==source_hash,
        'check_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'status':'Independent conditional synthetic geometry and exact tensor-transform checks; no physical observation',
        'passed':not failures,'failures':failures,
        'counts':{'preferred_relays':2400,'distinct_emitter_rule_events':2400,'arbitrary_angle_labs':1600,
                  'registered_parallel_transverse_labs_including_luminal':25,'exact_tensor_transforms':matrix_count},
        'residuals':{'max_preferred_relay_relative':max_relay_error,'min_preferred_arrival_ratio':min_preferred_ratio,
                     'max_emitter_rule_scaled':max_emitter_error,'max_lab_relative':max_lab_error,
                     'max_ray_support_equation_scaled':max_ray_error},
        'computed_transverse_discriminant_over4':str(quarter_discriminant)}
args.out.parent.mkdir(parents=True,exist_ok=True)
args.out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
raise SystemExit(0 if result['passed'] and result['candidate_unchanged_during_run'] else 1)
