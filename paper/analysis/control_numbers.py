#!/usr/bin/env python3
"""Macros for the controls of the SSU-discordance flag (stage3/artifact_control.py; added after review, point 1.6).

  python3 analysis/control_numbers.py --control ../benchmarks/copy_checks/artifact_control --paper .
Writes numbers_artcontrol.tex.
"""
import argparse
import csv
import os


def rd(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--control', required=True)
    ap.add_argument('--paper', default='.')
    a = ap.parse_args()
    s = {r['kind']: r for r in rd(os.path.join(a.control, 'artifact_control.summary.tsv'))}
    q = rd(os.path.join(a.control, 'artifact_control.queries.tsv'))
    m = {}

    def pct(r):
        return f"{float(r['flagged_pct']):.0f}\\%" if float(r['flagged_pct']) >= 10 else f"{float(r['flagged_pct']):.1f}\\%"
    m['AcRealN'] = s['real']['queries']
    m['AcRealFlag'] = s['real']['flagged']
    m['AcRealSSUOne'] = sum(r['kind'] == 'real' and r['ssu_pct'] != '' and float(r['ssu_pct']) >= 1 for r in q)
    m['AcErrN'] = s['errors_2']['queries']
    for rate, tag in (('0.5', 'Half'), ('1', 'One'), ('2', 'Two'), ('3', 'Three'), ('5', 'Five'), ('8', 'Eight')):
        r = s[f'errors_{rate}']
        m[f'AcErr{tag}'] = pct(r)
        m[f'AcErr{tag}CI'] = f"{float(r['ci95_low']):.0f}--{float(r['ci95_high']):.0f}\\%"
    jn = [s['chim_junction (other genus)'], s['chim_junction (same genus)']]
    sn = [s['chim_ssu (other genus)'], s['chim_ssu (same genus)']]
    m['AcJunctionN'] = sum(int(r['queries']) for r in jn)
    m['AcJunctionFlag'] = sum(int(r['flagged']) for r in jn)
    m['AcSSUChimN'] = sum(int(r['queries']) for r in sn)
    m['AcSSUChimFlag'] = sum(int(r['flagged']) for r in sn)
    m['AcSSUChimFloor'] = sum(r['kind'].startswith('chim_ssu') and r['ssu_pct'] != '' and float(r['ssu_pct']) >= 1 for r in q)
    m['AcSSUChimOtherFlag'] = s['chim_ssu (other genus)']['flagged']
    with open(os.path.join(a.paper, 'numbers_artcontrol.tex'), 'w') as fh:
        fh.write('% written by analysis/control_numbers.py from benchmarks/copy_checks/artifact_control\n')
        for k, v in m.items():
            fh.write(f'\\newcommand{{\\{k}}}{{{v}}}\n')
    print(m)


if __name__ == '__main__':
    main()
