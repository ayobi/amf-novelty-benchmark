#!/usr/bin/env bash
# Copy the files behind the paper from the AMF_test working folder into this repository.
#
#   bash stage_repo.sh ~/Downloads/AMF_test           dry run: lists what would be copied, changes nothing
#   bash stage_repo.sh ~/Downloads/AMF_test --copy    copies, removes LongGloDB sequences, writes MANIFEST.tsv
#
# Rules:
#   - code, edits files, protocols, locks, logs and result tables are copied;
#   - sequence files (FASTA/FASTQ/gz/zip) are left out everywhere except the culture reference v4,
#     which is derived from public GenBank records (Stefani et al. 2025);
#   - per-fold run folders (runs/) and files over 50 MB are left out and listed;
#   - LongGloDB sequences are never copied: any file under a LongGloDB path that holds a sequence
#     of 200 nt or more is removed after copying and listed;
#   - third-party references (MaarjAM, UNITE, Delavaux V18, LongGloDB) are not copied; their
#     checksums go to external_refs/CHECKSUMS.sha256.
set -euo pipefail

SRC=${1:?usage: bash stage_repo.sh <AMF_test folder> [--copy]}
MODE=${2:-dry}
SRC=$(cd "$SRC" && pwd)
DEST=$(cd "$(dirname "$0")" && pwd)
command -v rsync >/dev/null || { echo "rsync is needed: sudo apt install rsync"; exit 1; }
COPY=0; [ "$MODE" = "--copy" ] && COPY=1
[ $COPY = 1 ] && echo "MODE: copy" || echo "MODE: dry run (nothing is changed; add --copy to copy)"
TMPD=$(mktemp -d); trap 'rm -rf "$TMPD"' EXIT
echo "FROM: $SRC"; echo "TO:   $DEST"; echo

BASE=(rsync -a --prune-empty-dirs --max-size=50m
      --exclude='__pycache__/' --exclude='.ipynb_checkpoints/' --exclude='*.pyc' --exclude='.git/'
      --exclude='runs/' --exclude='*.zip' --exclude='*.bam')
[ $COPY = 1 ] || BASE+=(--dry-run)
NOSEQ=(--exclude='*.fasta' --exclude='*.fa' --exclude='*.fas' --exclude='*.fna'
       --exclude='*.fastq' --exclude='*.fq' --exclude='*.gz')
CODE=(--include='*/' --include='*.py' --include='*.sh' --include='*.R' --include='*.tsv'
      --include='*.md' --include='*.txt' --include='*.json' --include='*.yml' --include='*.yaml'
      --include='*.toml' --include='*.cfg' --exclude='*')

MAPPED=()
copy() {  # copy <source subfolder> <destination subfolder> [rsync filters...]
  local rel=$1 s="$SRC/$1" d="$DEST/$2"; shift 2
  MAPPED+=("$rel")
  if [ ! -d "$s" ]; then echo "== MISSING $s"; echo; return 0; fi
  echo "== $rel/  ->  $(basename "$DEST")/${d#$DEST/}/"
  find "$s" -type f -size +50M -not -path '*/runs/*' -printf '   SKIPPED (over 50 MB): %P\n' 2>/dev/null || true
  [ $COPY = 1 ] && mkdir -p "$d"
  if [ $COPY = 1 ]; then
    "${BASE[@]}" "$@" --out-format='   %n' "$s/" "$d/" | grep -v '/$' || true
  else
    mkdir -p "$TMPD/$rel"; "${BASE[@]}" "$@" --out-format='   %n' "$s/" "$TMPD/$rel/" | grep -v '/$' || true
  fi
  echo
}

copy amfsplit_stage3                   stage3                           "${CODE[@]}"
copy amfsplit_stage1                   stage1/prototypes                "${CODE[@]}"
copy amfsplit_v03/amfsplit_stage1_v03  stage1/amfsplit_v0.3.1           "${CODE[@]}"
copy ref/culture_ref/v4                reference/culture_ref_v4
copy bench/stefani                     benchmarks/stefani               "${NOSEQ[@]}"
copy amf_m4_step1                      benchmarks/m4_external_holdout   "${NOSEQ[@]}"
copy bench/longglodb                   benchmarks/longglodb             "${NOSEQ[@]}"

echo "== Folders in AMF_test not copied (tell Claude if any belong in the repo):"
for d in "$SRC"/*/ "$SRC"/bench/*/ "$SRC"/ref/*/; do
  [ -d "$d" ] || continue
  rel=${d#$SRC/}; rel=${rel%/}
  hit=0
  for m in "${MAPPED[@]}"; do
    case "$m" in "$rel"|"$rel"/*) hit=1;; esac
    case "$rel" in "$m"/*) hit=1;; esac
  done
  case "$rel" in bench|ref|ref/*) hit=1;; esac   # ref/* are listed with the third-party references
  [ $hit = 0 ] && printf '   %s  (%s)\n' "$rel" "$(du -sh "$d" 2>/dev/null | cut -f1)"
done
echo

echo "== Third-party references: checksums only (external_refs/CHECKSUMS.sha256)"
mapfile -t EXT < <(find "$SRC/ref" -type f \( -name '*.fasta' -o -name '*.fasta.gz' -o -name '*.fa' \
                    -o -name '*.fas' -o -name '*.fa.gz' \) -not -path '*/culture_ref/*' -not -path '*/.git/*' | sort)
for f in "${EXT[@]}"; do echo "   ${f#$SRC/}"; done
if [ $COPY = 1 ]; then
  mkdir -p "$DEST/external_refs"
  ( cd "$SRC" && for f in "${EXT[@]}"; do sha256sum "${f#$SRC/}"; done ) > "$DEST/external_refs/CHECKSUMS.sha256"
fi
echo

if [ $COPY = 1 ]; then
  echo "== Removing any LongGloDB-derived file that holds a sequence (200 nt or more)"
  n=0
  while IFS= read -r f; do
    case "$f" in *longglodb*|*lgdb*|*LongGloDB*|*glodb*) echo "   removed ${f#$DEST/}"; rm -f "$f"; n=$((n+1));;
                 *) echo "   KEPT, check: ${f#$DEST/} holds sequences (fine if from GenBank/Stefani)";;
    esac
  done < <(grep -rlE '[ACGTNacgtn]{200,}' "$DEST" --exclude-dir=.git --exclude-dir=reference \
           --exclude-dir=paper --exclude='*.sh' --exclude='MANIFEST.tsv' 2>/dev/null || true)
  echo "   $n removed"; echo

  echo "== Conda environments (env/)"
  if command -v conda >/dev/null; then
    mkdir -p "$DEST/env"
    for e in amfsplit amf-m4; do
      conda env export -n "$e" --from-history > "$DEST/env/$e.yml" 2>/dev/null \
        && echo "   env/$e.yml" || { rm -f "$DEST/env/$e.yml"; echo "   could not export $e"; }
    done
  else
    echo "   conda not found; export the environments by hand"
  fi
  echo

  echo "== Writing MANIFEST.tsv"
  ( cd "$DEST" && printf 'path\tbytes\tsha256\n' && \
    find . -type f -not -path './.git/*' -not -name MANIFEST.tsv -printf '%P\n' | LC_ALL=C sort | \
    while IFS= read -r p; do printf '%s\t%s\t%s\n' "$p" "$(stat -c %s "$p")" "$(sha256sum "$p" | cut -d' ' -f1)"; done
  ) > "$DEST/MANIFEST.tsv"
  echo "   $(($(wc -l < "$DEST/MANIFEST.tsv")-1)) files, $(du -sh --exclude=.git "$DEST" | cut -f1) in total"
fi
