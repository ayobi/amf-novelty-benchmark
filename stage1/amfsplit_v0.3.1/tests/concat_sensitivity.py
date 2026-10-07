#!/usr/bin/env python3
"""Concatemer screen on synthetic reads with ONT-like errors: false-positive and detection rates.
Run from the amfsplit folder: python3 tests/concat_sensitivity.py"""
import random, sys, time
sys.path.insert(0, 'tests'); sys.path.insert(0, '.')
import test_v03 as T
import amfsplit as A

rng = random.Random(11)

def errify(s, rate):
    out = []
    for ch in s:
        if rng.random() < rate:
            r = rng.random()
            if r < 2 / 3: out.append(rng.choice([b for b in 'ACGT' if b != ch]))
            elif r < 5 / 6: out.append(ch + rng.choice('ACGT'))
        else:
            out.append(ch)
    return ''.join(out)

by_site = {}
for name, fwd, fs, rev, rs, lo, hi in [(w[0], w[1], w[2], w[3], w[4], w[5], w[6]) for w in T.PANEL]:
    by_site.setdefault(fs, fwd); by_site.setdefault(T.rc(rs), rev)
primers = {n: s for s, n in by_site.items()}
print(f"{len(primers)} distinct primers; {'error':>6s} {'single-copy flagged':>21s} {'tandem flagged':>16s} {'inverted flagged':>18s}  time/read")
for rate in (0.03, 0.08, 0.12, 0.15):
    frac = 0.20 if rate >= 0.12 else 0.15
    n = 150
    res = {'single': 0, 'tandem': 0, 'inverted': 0}
    t0 = time.time()
    for _ in range(n):
        ssu, its, lsu = T.make_copy(); full = ssu + its + lsu
        other = T.mutate(full, int(0.02 * len(full)))   # a second copy of the same molecule, ~2% different
        for kind, seq in (('single', full), ('tandem', full + other), ('inverted', full + T.rc(other))):
            res[kind] += A.concatemer_screen(errify(seq, rate), primers, frac, 150, 2, 0.05)[0]
    dt = (time.time() - t0) / (3 * n)
    print(f"        {'':>0s}{100 * rate:>5.0f}% (k frac {frac}) {100 * res['single'] / n:>14.1f}% {100 * res['tandem'] / n:>15.1f}% {100 * res['inverted'] / n:>17.1f}%  {1000 * dt:.0f} ms")
