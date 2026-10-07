"""Preflight, no external tools: add the 280 test ASVs (with silent evidence and empty placements) to the
frozen data, rescore all 37 M4 species folds with the test cultures excluded, and require the stored
M4 v2 results back exactly. Usage: python3 tests/preflight_parity.py ../bench/longglodb"""
import csv, json, sys, time
from pathlib import Path
M4 = Path(__file__).resolve().parent.parent; LGDB = Path(sys.argv[1]); TMP = Path(sys.argv[2]) if len(sys.argv) > 2 else M4 / "runs_lgdb_preflight"
sys.path.insert(0, str(M4))
import lgdb_test as L
inputs = L.Inputs(LGDB)
data_ref = L.v1.Data(L.FROZEN / 'inputs/benchmark')
tc, tcp = L.test_set(inputs, set(data_ref.species.values()))
ev = {'VT': {}, 'SH': {}, 'VT_spanned': {}}
for x in tcp:
    for k in ev:
        ev[k][x] = {'unit_pident': {}, 'anchor_pident': None, 'best': ('', 0.)}
data = L.augmented_data(tc, tcp, inputs.win, ev, 'VT')
cols = list(next(iter(L.read_tsv(L.FROZEN / 'placement/summary/lsu_placements.copies.tsv'))).keys())
rows = []
for x in tcp:
    r = {k: '' for k in cols}
    r.update(accession=x, bridge_alignment_eligible='0', credible_mass_target='0.95', credible_species_clade_units='')
    rows.append(r)
TMP.mkdir(parents=True, exist_ok=True)
L.write_tsv(TMP / 'mock_test_placements.tsv', rows)
pa = L.placements_for(data, L.FROZEN / 'placement/summary/lsu_placements.copies.tsv', TMP / 'mock_test_placements.tsv')
t0 = time.time(); rule = L.ItsRule(data, L.v2.QUANTILE); print(f'ItsRule on {len(data.cop)} copies: {time.time()-t0:.0f}s', flush=True)
t0 = time.time()
par = L.parity_m4(data, rule, frozenset(tc), pa, M4 / 'runs')
print(f"parity: {par['copies']} copies in {par['folds']} folds reproduced exactly ({time.time()-t0:.0f}s)", flush=True)
before, after = L.score(data, rule, sorted(tcp), frozenset(tc), pa)
from collections import Counter
print('mock test scoring:', len(after), 'rows;', dict(Counter((tcp[r['accession']]['role'], r['status']) for r in after)))
print('VT states:', dict(Counter(r['VT_state'] for r in after)), '| groups:', dict(Counter(r['group'] for r in after)))
assert all(r['group'] == ('named, species elsewhere' if tcp[r['accession']]['role'] == 'known species' else "unnamed ('sp.')" if tcp[r['accession']]['role'] == 'genus only' else 'named, species held out') for r in after)
print('truth groups as expected for every test role')
