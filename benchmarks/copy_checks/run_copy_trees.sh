#!/usr/bin/env bash
# Per-marker trees and monophyly (PROTOCOL_COPY_CHECKS.md, section 1). Run from the repository root:
#   bash benchmarks/copy_checks/run_copy_trees.sh            # THREADS defaults to nproc, JOBS to 4
#   THREADS=32 JOBS=4 bash benchmarks/copy_checks/run_copy_trees.sh
# Steps: prep -> MAFFT L-INS-i (12 alignments) -> concat -> IQ-TREE (16 trees) -> score.
# A finished IQ-TREE run is skipped on a rerun and an interrupted one resumes from its checkpoint,
# so the script can be restarted.
set -euo pipefail

REF=${REF:-reference/culture_ref_v4}
OUT=${OUT:-benchmarks/copy_checks/trees}
COMPLEXES=${COMPLEXES:-paper/analysis/v4.complexes.tsv}
THREADS=${THREADS:-$(nproc)}
JOBS=${JOBS:-4}
PER=$(( THREADS / JOBS )); [ "$PER" -ge 1 ] || PER=1
SEED=20261011
IQTREE=${IQTREE:-$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)}
[ -n "$IQTREE" ] || { echo "IQ-TREE not found (iqtree3, iqtree2 or iqtree)"; exit 1; }
command -v mafft >/dev/null || { echo "mafft not found"; exit 1; }

mkdir -p "$OUT"
{
  printf 'mafft\t%s\n' "$(mafft --version 2>&1 | head -1)"
  printf 'iqtree\t%s\n' "$("$IQTREE" --version 2>&1 | head -1)"
  printf 'python\t%s\n' "$(python3 --version 2>&1)"
  printf 'started\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
  printf 'host\t%s\n' "$(hostname)"
} > "$OUT/versions.txt"
cat "$OUT/versions.txt"

python3 stage3/copy_trees.py prep --ref "$REF" --out "$OUT"

# 1. MAFFT L-INS-i on the unique sequences of each order and marker
ls "$OUT"/*/*.unique.fasta | xargs -P "$JOBS" -I{} bash -c '
  f={}; o=${f%.unique.fasta}.aln.fasta
  if [ ! -s "$o" ]; then
    mafft --localpair --maxiterate 1000 --thread '"$PER"' --quiet "$f" > "$o.tmp" && mv "$o.tmp" "$o"
    echo "aligned $o"
  fi'

python3 stage3/copy_trees.py concat --out "$OUT"

# 2. IQ-TREE: single markers, then the partitioned linked alignment
run_iq() {
  local aln=$1 prefix=$2 part=${3:-}
  # IQ-TREE writes .treefile during the search too; only the final log line marks a finished run
  if grep -qs '^Total wall-clock time used' "$prefix.log"; then echo "done $prefix"; return; fi
  if [ -n "$part" ]; then
    "$IQTREE" -s "$aln" -p "$part" -m MFP -B 1000 --alrt 1000 -czb --seed "$SEED" -T "$PER" --prefix "$prefix" --quiet
  else
    "$IQTREE" -s "$aln" -m MFP -B 1000 --alrt 1000 -czb --seed "$SEED" -T "$PER" --prefix "$prefix" --quiet
  fi
  echo "tree $prefix"
}
export -f run_iq; export IQTREE SEED PER
{
  for d in "$OUT"/*/; do
    for m in ssu_flank full_its lsu_flank; do echo "${d}${m}.aln.fasta ${d}${m}"; done
    echo "${d}linked.aln.fasta ${d}linked ${d}linked.partitions.nex"
  done
} | xargs -P "$JOBS" -L 1 bash -c 'run_iq "$@"' _

# 3. Monophyly of cultures, named species and complexes
python3 stage3/copy_trees.py score --ref "$REF" --out "$OUT" --complexes "$COMPLEXES"
printf 'finished\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" >> "$OUT/versions.txt"
