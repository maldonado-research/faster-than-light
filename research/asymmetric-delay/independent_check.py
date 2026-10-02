#!/usr/bin/env python3
"""Independent unequal-leg conditional SR checks; Python >=3.10, stdlib only.

Inputs are deterministic mathematical parameters, not physical observations.
Run: python3 reproduce.py --out results.json
The companion AUDIT.md records assumptions and exact derivations.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import random


def factor(speed, beta):
    return speed * math.sqrt((1-beta)*(1+beta)) / (speed-beta)


def ratio(u, v, beta, delay):
    return factor(v, beta) * (factor(u, beta)+delay)


def direct_event(u, v, beta, delay):
    # Independently construct the Bob departure event in Alice coordinates.
    gamma = 1/math.sqrt((1-beta)*(1+beta))
    reception_t = u/(u-beta)
    departure_t = reception_t + gamma*delay
    departure_x = beta*departure_t
    # Bob-future elapsed coordinate time along the reply in Bob's frame.
    bob_elapsed = departure_x/(gamma*(v-beta))
    return departure_t + gamma*bob_elapsed*(1-beta*v)


def allowance(u, v, beta):
    return 1/factor(v, beta)-factor(u, beta)


def decimal_gap(u_text, v_text, delay_text):
    """Root in q=exp(-rapidity), using a polynomial independent of gap route.

    q=sqrt((1-beta)/(1+beta)). The polynomial arises from multiplying
    positive denominators, not squaring a radical equation.
    """
    with localcontext() as ctx:
        ctx.prec=110
        u,v,d=map(Decimal,(u_text,v_text,delay_text))
        qc=((u-1)*(v-1)/((u+1)*(v+1))).sqrt()
        if not d:
            return +(2*qc*qc/(1+qc*qc))
        lo,hi=Decimal(0),qc
        for _ in range(380):
            q=(lo+hi)/2
            polynomial=((u+1)*(v+1)*q**4 - 2*d*v*(u+1)*q**3
                        - 2*(u*v+1)*q*q - 2*d*v*(u-1)*q + (u-1)*(v-1))
            if polynomial>0:
                lo=q
            else:
                hi=q
        q=(lo+hi)/2
        return +(2*q*q/(1+q*q))


def run(root_result_path):
    failures=[]
    rng=random.Random(20261002)
    max_event_relative_error=0.
    for i in range(20000):
        u=1+10**rng.uniform(-3,2)
        v=1+10**rng.uniform(-3,2)
        beta=rng.random()*.999999
        delay=0. if i%5==0 else 10**rng.uniform(-3,2)
        expected=ratio(u,v,beta,delay)
        actual=direct_event(u,v,beta,delay)
        max_event_relative_error=max(max_event_relative_error,abs(actual-expected)/expected)
        if not math.isclose(actual,expected,rel_tol=2e-10,abs_tol=2e-10):
            failures.append(f"Event mismatch {i}")

    exact_cases=0
    for beta,s in [(Fraction(3,5),Fraction(4,5)),(Fraction(4,5),Fraction(3,5)),
                   (Fraction(5,13),Fraction(12,13)),(Fraction(12,13),Fraction(5,13))]:
        for u in map(Fraction,("3/2","2","3","10")):
            for v in map(Fraction,("3/2","2","3","10")):
                for d in map(Fraction,("0","1/5","1/2","2")):
                    departure_t=u/(u-beta)+d/s
                    q=beta*departure_t*s/(v-beta)
                    actual=departure_t+q*(1-beta*v)/s
                    a=u*s/(u-beta);b=v*s/(v-beta)
                    exact_cases+=1
                    if actual!=b*(a+d):
                        failures.append("Exact event mismatch")
                    boundary=1/b-a
                    if b*(a+boundary)!=1:
                        failures.append("Exact allowance boundary mismatch")

    threshold_cases=side_cases=0
    max_threshold_residual=Decimal(0)
    threshold_rows=[]
    with localcontext() as ctx:
        ctx.prec=85
        for ut in ("1.001","1.1","2","10","100"):
            for vt in ("1.001","1.1","2","10","100"):
                for dt in ("0","0.1","0.5","1","10","100"):
                    u,v,d=map(Decimal,(ut,vt,dt))
                    e=decimal_gap(ut,vt,dt)
                    def at_gap(gap):
                        s=(gap*(2-gap)).sqrt()
                        return v*s/(v-1+gap)*(u*s/(u-1+gap)+d)
                    residual=abs(at_gap(e)-1)
                    max_threshold_residual=max(max_threshold_residual,residual)
                    threshold_cases+=1
                    if residual>Decimal("1e-60"):
                        failures.append(f"Threshold residual {ut},{vt},{dt}")
                    for m,before in ((Decimal(".999"),True),(Decimal("1.001"),False)):
                        side_cases+=1
                        if (at_gap(e*m)<1)!=before:
                            failures.append(f"Threshold side {ut},{vt},{dt}")
                    if dt in ("0","1","100") and ut in ("1.1","2") and vt in ("2","10"):
                        threshold_rows.append({"U":ut,"V":vt,"d":dt,"one_minus_beta_threshold":str(e)})

    cap_checks=0
    for u in (1.01,1.1,2,10,100):
        for v in (1.01,1.1,2,10,100):
            bc=(u+v)/(u*v+1)
            if not (bc>1/u and bc>1/v and bc<1):
                failures.append("Immediate threshold location")
            for mix in (.25,.5,.75):
                cap=bc+(1-bc)*mix
                bound=allowance(u,v,cap)
                for i in range(501):
                    beta=cap*i/500
                    if ratio(u,v,beta,bound)<1-2e-10:
                        failures.append("Cap boundary failure")
                    if ratio(u,v,beta,bound+1)<1-2e-10:
                        failures.append("Cap above boundary failure")
                    cap_checks+=2
                if not ratio(u,v,cap,.999*bound)<1:
                    failures.append("Cap strict lower side failure")
                cap_checks+=1

    luminal_checks=0
    for u,v in ((1.,2.),(2.,1.),(1.,1.),(1.,100.),(100.,1.)):
        for beta in (0.,.1,.5,.9,.99,.999999):
            for d in (0.,.1,1.,100.):
                if ratio(u,v,beta,d)<1-1e-12:
                    failures.append("Luminal leg created modeled return-before-send")
                luminal_checks+=1

    asymptotic_checks=0
    asymptotic_rows=[]
    for u,v in ((1.1,2.),(2.,3.),(100.,1.01)):
        for d in (0.,.1,2.):
            for e in (1e-8,1e-10,1e-12):
                s=math.sqrt(e*(2-e))
                actual=v*s/(v-1+e)*(u*s/(u-1+e)+d)
                approximation=d*math.sqrt(2)*v/(v-1)*math.sqrt(e)+2*u*v/((u-1)*(v-1))*e
                relative=abs(actual-approximation)/actual
                if relative>2e-6:
                    failures.append("Near-beta-one asymptotic")
                asymptotic_checks+=1
                if e==1e-12:
                    asymptotic_rows.append({"U":u,"V":v,"d":d,"gap":e,"relative_error_two_terms":relative})

    large_delay_rows=[]
    for ut,vt in (("1.1","2"),("2","3"),("3","2"),("100","1.01")):
        v=float(vt)
        coefficient=(v-1)**2/(2*v*v)
        for dt in ("100","1000","10000"):
            e=float(decimal_gap(ut,vt,dt));d=float(dt)
            relative=abs(d*d*e-coefficient)/coefficient
            if dt=="10000" and relative>1e-6:
                failures.append("Large-delay threshold asymptotic")
            large_delay_rows.append({"U":ut,"V":vt,"d":dt,"d2_times_threshold_gap":d*d*e,
                                     "predicted_limit":coefficient,"relative_error":relative})

    beta=Fraction(4,5);s=Fraction(3,5);d=Fraction(1,5)
    def exact_ratio(u,v):
        return v*s/(v-beta)*(u*s/(u-beta)+d)
    fast_outgoing=exact_ratio(Fraction(3),Fraction(2))
    fast_return=exact_ratio(Fraction(2),Fraction(3))
    if fast_outgoing!=Fraction(56,55) or fast_return!=Fraction(54,55):
        failures.append("Exact orientation counterexample")

    root_result=json.loads(root_result_path.read_text())
    comparisons=[]
    max_root_gap_relative_difference=Decimal(0)
    with localcontext() as ctx:
        ctx.prec=110
        for row in root_result["thresholds"]:
            independent=decimal_gap(row["U"],row["V"],row["d"])
            reported=Decimal(row["gap_1_minus_beta_threshold"])
            relative=abs(independent-reported)/independent
            max_root_gap_relative_difference=max(max_root_gap_relative_difference,relative)
            if relative>Decimal("1e-60"):
                failures.append("Independent q-polynomial comparison failed")
            comparisons.append({"U":row["U"],"V":row["V"],"d":row["d"],
                                "independent_gap":str(independent),"relative_difference":str(relative)})
    root_script=root_result_path.parent/"reproduce.py"
    root_hash_matches=(root_result["script_sha256"]==hashlib.sha256(root_script.read_bytes()).hexdigest())
    if not root_hash_matches:
        failures.append("Root script hash differs from root result")

    return {"run_utc":datetime.now(timezone.utc).isoformat(),
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "scope":"Internal conditional mathematical audit; deterministic parameters, no observations, physical FTL channel or novelty claim",
            "passed":not failures,"failures":failures,
            "checks":{"seed":20261002,"floating_event_cases":20000,"max_event_relative_error":max_event_relative_error,
                      "exact_event_and_allowance_cases":exact_cases,"decimal_precision":85,
                      "decimal_threshold_cases":threshold_cases,"threshold_side_cases":side_cases,
                      "max_threshold_residual":str(max_threshold_residual),"beta_cap_checks":cap_checks,
                      "luminal_leg_checks":luminal_checks,"near_beta_one_asymptotic_checks":asymptotic_checks,
                      "independent_threshold_solver_precision":110,"root_threshold_comparison_cases":len(comparisons),
                      "max_root_gap_relative_difference":str(max_root_gap_relative_difference)},
            "root_result_provenance":{"file":root_result_path.name,"sha256":hashlib.sha256(root_result_path.read_bytes()).hexdigest(),
                                      "root_script_hash_matches":root_hash_matches},
            "root_gap_comparisons":comparisons,
            "threshold_examples":threshold_rows,"near_beta_one_examples":asymptotic_rows,
            "large_delay_threshold_examples":large_delay_rows,
            "exact_orientation_counterexample":{"beta":"4/5","d":"1/5","R_U3_V2":str(fast_outgoing),"R_U2_V3":str(fast_return)}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--root-result",type=Path,
                        default=Path(__file__).with_name("RESULTS.json"))
    args=parser.parse_args()
    result=run(args.root_result)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"passed":result["passed"],"failures":result["failures"],"checks":result["checks"]}))
    raise SystemExit(0 if result["passed"] else 1)


if __name__=="__main__":
    main()
