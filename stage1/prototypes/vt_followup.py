#!/usr/bin/env python3
"""Follow-up checks on the Stage 1 VT results. Run from amfsplit_stage1, next to vt_check.py.

1. Paraglomus brasilianum BR105_2435A: are its VTX00281 copies P. laccatum (mixed culture)?
2. Septoglomus viscosum VTX00409 copies: one variant class across two cultures, or a contaminant?
3. PX215117.1 (Gigaspora rosea): chimera, or a uniformly divergent SSU?
"""
import csv, edlib
from vt_check import read_fasta
S = '../bench/stefani'; P = f'{S}/out/windows/stefani'
T = {r['accession']: r for r in csv.DictReader(open(f'{S}/vt_check.copies.tsv'), delimiter='\t')}
seq = {k: read_fasta(f'{P}.{k}.fasta') for k in ('its1', 'its2', 'ssu_flank')}

def d(a, b):  # % edit distance, shorter aligned inside longer (end gaps free)
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    return 100 * edlib.align(q, t, mode='HW', task='distance')['editDistance'] / len(q)

def near(a, pool, k, s=None):
    s = s or seq[k].get(a)
    ds = [d(s, seq[k][b]) for b in pool if b != a and b in seq[k]] if s else []
    return f'{min(ds):.1f}' if ds else '-'

pick = lambda f: [a for a, r in T.items() if f(r)]
grp = lambda g: pick(lambda r: r['group'] == g)

print('1. Paraglomus brasilianum BR105_2435A: % to nearest copy in each set (ITS2 / ITS1)')
br = pick(lambda r: r['group'].endswith('BR105_2435A'))
lac = pick(lambda r: r['organism'] == 'Paraglomus laccatum')
bra = pick(lambda r: r['organism'] == 'Paraglomus brasilianum' and not r['group'].endswith('BR105_2435A'))
print(f"   {'copy':12s}{'VT':10s}{'P. laccatum':>14s}{'other P. brasil.':>18s}{'BR105, other VT':>18s}")
for a in sorted(br, key=lambda a: T[a]['vt']):
    oth = [b for b in br if T[b]['vt'] != T[a]['vt']]
    c = [f"{near(a, p, 'its2')} / {near(a, p, 'its1')}" for p in (lac, bra, oth)]
    print(f"   {a:12s}{T[a]['vt']:10s}{c[0]:>14s}{c[1]:>18s}{c[2]:>18s}")

print('\n2. Septoglomus viscosum VTX00409 copies: % distance to each other, and to the nearest main-type (VTX00063) copy of their own culture')
sv = pick(lambda r: r['vt'] == 'VTX00409')
for k in ('ssu_flank', 'its1', 'its2'):
    pair = '  '.join(f'{a[:8]}-{b[:8]} {d(seq[k][a], seq[k][b]):.1f}' for i, a in enumerate(sv) for b in sv[i + 1:])
    own = '  '.join(f"{a[:8]} {near(a, [b for b in grp(T[a]['group']) if T[b]['vt'] == 'VTX00063'], k)}" for a in sv)
    print(f'   {k:10s} each other: {pair} | own main type: {own}')
out = pick(lambda r: r['organism'] != 'Septoglomus viscosum')
for a in sv:
    b = min((b for b in out if b in seq['its2']), key=lambda b: d(seq['its2'][a], seq['its2'][b]))
    print(f"   {a} nearest ITS2 outside S. viscosum: {T[b]['organism']} {d(seq['its2'][a], seq['its2'][b]):.1f}")

print('\n3. PX215117.1 (Gigaspora rosea): % to the nearest other G. rosea copy, by segment')
a = 'PX215117.1'; s = seq['ssu_flank'][a]
gr = pick(lambda r: r['organism'] == 'Gigaspora rosea' and r['accession'] != a)
n = 4; w = len(s) // n
segs = [(i * w, len(s) if i == n - 1 else (i + 1) * w) for i in range(n)]
print('   SSU flank ' + '  '.join(f'{x + 1}-{y}: {near(a, gr, "ssu_flank", s[x:y])}' for x, y in segs))
print(f"   ITS1 {near(a, gr, 'its1')}   ITS2 {near(a, gr, 'its2')}")
