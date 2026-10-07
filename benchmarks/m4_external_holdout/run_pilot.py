#!/usr/bin/env python3
"""Run only the Glomus chinense external-reference holdout."""
import argparse
import json
from pathlib import Path
import sys
import zipfile
import external_holdout as h
from verify_code import verify


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--threads',type=int,default=8)
    ap.add_argument('--mafft',default='mafft')
    ap.add_argument('--raxml-ng',default='raxml-ng')
    ap.add_argument('--epa-ng',default='epa-ng')
    ap.add_argument('--check-only',action='store_true')
    args=ap.parse_args()
    if args.threads<1:ap.error('threads must be positive')
    verify();frozen_hash=h.verify_frozen();tools=h.tool_info(args)
    plan=h.load_plan('species',['Glomus chinense'])[0]
    base_key=h.digest({'frozen':frozen_hash,'protocol':h.sha(h.ROOT/'PROTOCOL.md'),'plan':h.sha(h.ROOT/'holdout_plan.json'),
                       'algorithm':'fresh-reference-fixed-topology-v1','versions':{k:v['version'] for k,v in tools.items()}})
    signature=h.digest({'base':base_key,'plan':plan})
    out=h.ROOT/'runs';fold=out/plan['key']
    h.valid_stage(fold/'01_reference',signature+':reference')
    h.valid_stage(fold/'02_model',signature+':model')
    if args.check_only:
        print('Pilot preflight passed: frozen inputs, tool versions and saved reference/model checkpoints match.')
        return
    h.write_json(out/'pilot_invocation.json',{'argv':sys.argv,'script_sha256':h.sha(__file__),
                   'engine_sha256':h.sha(h.ROOT/'external_holdout.py'),'base_key':base_key,'tools':tools,'threads':args.threads})
    h.run_fold(plan,out,tools,args.threads,base_key)
    m=json.loads((fold/'06_evaluation/metrics.json').read_text())
    summary={'pilot_complete':True,'target':'Glomus chinense','reference_tips_removed':len(plan['removed_ref_ids']),
             'full_benchmark_complete':False,'main_copies':m['external_holdout']['held_n'],
             'full_reference_false_known':m['full_reference']['false_known'],
             'external_holdout_false_known':m['external_holdout']['false_known'],
             'metrics':m}
    h.write_json(out/'PILOT_SUMMARY.json',summary)
    destination=h.ROOT/'amf_m4_step1_results.zip'
    h.collect(out,destination)
    with zipfile.ZipFile(destination,'a',zipfile.ZIP_DEFLATED) as z:
        for name in ('PILOT_SUMMARY.json','pilot_invocation.json'):
            z.write(out/name,'results/'+name)
    print(f"Pilot complete. False-known calls: {summary['full_reference_false_known']}/{summary['main_copies']} with full V18 -> {summary['external_holdout_false_known']}/{summary['main_copies']} after removal.")
    print(f'Upload {destination.name}. We will inspect this one result before expanding.')


if __name__=='__main__':main()
