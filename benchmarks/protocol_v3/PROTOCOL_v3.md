#!/usr/bin/env bash
# Protocol v3, procedure step 6: the protocol lock.
# Puts PROTOCOL_v3.md (the doc's own Markdown export), the eligibility rules, the per-stratum predictions as rates,
# the name-check edits, every Protocol v3 input, the arm code and the v0.1 commit into protocol_hashes.txt;
# commits, and tags both repositories protocol-v3-lock.
set -euo pipefail
B=~/projects/amf-novelty-benchmark
A=~/projects/amfsplit
P3=$B/benchmarks/protocol_v3
MD=$(ls -t ~/Downloads/*rotocol*v3*.md 2>/dev/null | head -1 || true)
[ -n "$MD" ] && [ -f "$MD" ] || { echo "Download the protocol doc as Markdown to ~/Downloads first."; exit 1; }
echo "protocol export: $MD"

# both repositories current, with nothing uncommitted
for r in "$B" "$A"; do
  git -C "$r" pull -q origin main
  [ -z "$(git -C "$r" status --porcelain --untracked-files=no)" ] || { echo "$r has uncommitted changes"; exit 1; }
done
[ "$(git -C "$A" rev-parse HEAD)" = "fa3404f3c372116e4fb5267726d6a9e28716d09f" ] || echo "note: amfsplit HEAD is not fa3404f; the lock records the HEAD it finds"

cp "$MD" "$P3/PROTOCOL_v3.md"
python3 - "$P3" <<'PY'
import sys, csv
p3 = sys.argv[1]
md = open(f'{p3}/PROTOCOL_v3.md', encoding='utf-8').read()
start = md.index('**Confirmation set, by eligibility rule.**')
end = md.index('**Measurement E', start)
open(f'{p3}/eligibility_rules.md', 'w', encoding='utf-8').write(
    '# Protocol v3: eligibility rules (extracted from PROTOCOL_v3.md)\n\n' + md[start:end].rstrip() + '\n')
rows = list(csv.DictReader(open(f'{p3}/rehearsal/adopt_scores.tsv'), delimiter='\t'))
with open(f'{p3}/predictions_v3.tsv', 'w') as f:
    f.write('arm\toutcome\tdevelopment_n\tfrozen_rate\tarm_rate\tdiff_low_points\tdiff_high_points\n')
    for r in rows:
        n = int(r['n'])
        fr, ar = (int(r['frozen']) / n, int(r['arm']) / n) if n else (0.0, 0.0)
        f.write(f"{r['design']}\t{r['outcome']}\t{n}\t{fr:.6f}\t{ar:.6f}\t{r['diff_low']}\t{r['diff_high']}\n")
print(f'eligibility rules ({end - start} characters) and {len(rows)} predictions written')
PY

cd "$P3"
{
  echo "# Protocol v3 lock, $(date -u +%Y-%m-%d) (procedure step 6). SHA-256 of each locked file, then the code commits."
  sha256sum PROTOCOL_v3.md eligibility_rules.md predictions_v3.tsv r4_name_checks.tsv E_PREDICTION.md E_prediction.lock \
            r6_triage_positives.tsv r6_triage_positives.summary.json r6_negatives_triage.tsv.gz \
            r6_negatives_triage.summary.json r6_negatives_named_rows.tsv.gz
  find withheld_k rehearsal -type f | sort | xargs sha256sum
  echo "# arm code: amfsplit $(git -C "$A" rev-parse HEAD) ($(git -C "$A" log -1 --format=%s | cut -c1-90))"
  ( cd "$A" && git ls-files crates/amfsplit-naming/src crates/amfsplit/src | sort | xargs sha256sum | sed 's/^/arm  /' )
  echo "# v0.1: amfsplit $(git -C "$A" rev-parse 'v0.1.0-rc.1^{commit}') (tag v0.1.0-rc.1)"
} > protocol_hashes.txt
echo "protocol_hashes.txt: $(wc -l < protocol_hashes.txt) lines"

printf '| `PROTOCOL_v3.md`, `eligibility_rules.md`, `predictions_v3.tsv`, `protocol_hashes.txt` | %s | The protocol lock (procedure step 6): the protocol as exported from its document, its eligibility rules, the per-stratum predictions as rates (from the dress rehearsal), and the SHA-256 of every locked file with the arm code and v0.1 commits. Tagged protocol-v3-lock in both repositories |\n' "$(date -u +%Y-%m-%d)" >> README.md

cd "$B"
git add benchmarks/protocol_v3
git commit -q -m "Protocol v3: the protocol lock (procedure step 6): PROTOCOL_v3.md, eligibility rules, per-stratum predictions and protocol_hashes.txt"
git tag -a protocol-v3-lock -m "Protocol v3 locked $(date -u +%Y-%m-%d); arm code amfsplit $(git -C "$A" rev-parse --short HEAD)"
git push -q origin main protocol-v3-lock
git -C "$A" tag -a protocol-v3-lock -m "Arm code of Protocol v3, locked $(date -u +%Y-%m-%d); hashes in amf-novelty-benchmark benchmarks/protocol_v3/protocol_hashes.txt"
git -C "$A" push -q origin protocol-v3-lock
echo "locked: benchmark $(git -C "$B" rev-parse --short HEAD), amfsplit $(git -C "$A" rev-parse --short HEAD)"
echo "protocol_hashes.txt: $(sha256sum "$P3/protocol_hashes.txt" | cut -c1-64)"
echo "PROTOCOL_v3.md:      $(sha256sum "$P3/PROTOCOL_v3.md" | cut -c1-64)"
