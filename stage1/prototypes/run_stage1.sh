#!/usr/bin/env bash
# Stage 1 driver (v0.2) for one sample:
#   ITSxRust anchors -> amfsplit windows -> LSU validation -> AMF triage -> AMF-only subsets
#
#   run_stage1.sh <reads.fq.gz> <sample_name> <hifi|ont> [outdir]
#
# Needs on PATH: itsxrust (>=0.3.0), nhmmer (HMMER 3), vsearch, python3 with edlib.
# Env:
#   HMM           ITSx fungal profile (F.hmm)                                   [required]
#   DELAVAUX_DB   LROR-FLR2 reference (triage via its 4 outgroups if no TRIAGE_DB)
#   TRIAGE_DB     broader LSU reference with lineage in headers (e.g. EUKARYOME LSU)
#   TRIAGE_PATTERN regex marking AMF targets in TRIAGE_DB headers   [Glomeromycota]
#   MAARJAM_DB    MaarjAM reference; checked against the AMF-only VT windows
#   ANCHOR_MODE   tolerant (default) | strict. Anchors only delimit regions; taxonomy is
#                 triage's job, so tolerant thresholds keep divergent lineages anchored.
#   THREADS [8], MAX_ERR_FRAC [0.15 hifi / 0.20 ont], MIN_ID [90 hifi / 85 ont]
set -euo pipefail

READS=${1:?reads}; SAMPLE=${2:?sample name}; PRESET=${3:?hifi or ont}; OUT=${4:-stage1_out}
HERE=$(cd "$(dirname "$0")" && pwd)
HMM=${HMM:?set HMM to the ITSx fungal F.hmm}
PANEL=${PANEL:-$HERE/panels/amf_legacy_windows.tsv}
THREADS=${THREADS:-8}
ANCHOR_MODE=${ANCHOR_MODE:-tolerant}
if [[ "$PRESET" == "ont" ]]; then
  MAX_ERR_FRAC=${MAX_ERR_FRAC:-0.20}; MIN_ID=${MIN_ID:-85}
else
  MAX_ERR_FRAC=${MAX_ERR_FRAC:-0.15}; MIN_ID=${MIN_ID:-90}
fi
ITSX_THRESH=()
if [[ "$ANCHOR_MODE" == "tolerant" ]]; then
  ITSX_THRESH=(--inc-e 1e-3 --min-anchor-score 15 --max-anchor-evalue 1e-3)
fi

mkdir -p "$OUT/itsx" "$OUT/windows"
W="$OUT/windows/$SAMPLE"

# 1) one nhmmer pass; ITS1/ITS2/full plus anchors, then SSU/LSU flanks from the saved tblout
itsxrust extract -i "$READS" --hmm "$HMM" --preset "$PRESET" "${ITSX_THRESH[@]}" \
  --hmmer-cpu "$THREADS" --region all -o "$OUT/itsx/$SAMPLE" --output-format fasta --plain-ids \
  --tblout "$OUT/itsx/$SAMPLE.tblout" --anchors-tsv "$OUT/itsx/$SAMPLE.anchors.tsv" \
  --qc-json "$OUT/itsx/$SAMPLE.qc.json" --write-skipped "$OUT/itsx/$SAMPLE.skipped.fasta"
for region in ssu lsu; do
  itsxrust extract -i "$READS" --tblout-existing "$OUT/itsx/$SAMPLE.tblout" --preset "$PRESET" \
    "${ITSX_THRESH[@]}" --region "$region" -o "$OUT/itsx/$SAMPLE.$region.fasta" \
    --output-format fasta --plain-ids
done

# 2) multi-marker records + legacy windows (open-ended where a primer site lies outside the amplicon)
python3 "$HERE/amfsplit.py" -i "$READS" -a "$OUT/itsx/$SAMPLE.anchors.tsv" -p "$PANEL" \
  -o "$W" --max-err-frac "$MAX_ERR_FRAC" --only-anchored --open-ends

# 3) triage on the LSU window: broad reference if given, else the Delavaux outgroups
LSU="$W.LSU_LROR-FLR2.fasta"
if [[ -n "${TRIAGE_DB:-}" && -s "$LSU" ]]; then
  python3 "$HERE/validate_windows.py" -w "$LSU" -r "$TRIAGE_DB" -o "$W.LSU_LROR-FLR2.vs_triage" \
    --threads "$THREADS" > /dev/null
  python3 "$HERE/triage.py" --hits "$W.LSU_LROR-FLR2.vs_triage.hits.tsv" --records "$W.records.tsv" \
    --fasta-prefix "$W" --amf-pattern "${TRIAGE_PATTERN:-Glomeromycota}" --min-id "$MIN_ID" --out "$W"
elif [[ -n "${DELAVAUX_DB:-}" && -s "$LSU" ]]; then
  python3 "$HERE/validate_windows.py" -w "$LSU" -r "$DELAVAUX_DB" -o "$W.LSU_LROR-FLR2.vs_delavaux" \
    --threads "$THREADS" > /dev/null
  python3 "$HERE/triage.py" --hits "$W.LSU_LROR-FLR2.vs_delavaux.hits.tsv" --records "$W.records.tsv" \
    --fasta-prefix "$W" --min-id "$MIN_ID" --out "$W"
fi

# 4) legacy-reference checks on the AMF-only subsets
if [[ -n "${DELAVAUX_DB:-}" && -s "$W.LSU_LROR-FLR2.amf.fasta" ]]; then
  python3 "$HERE/validate_windows.py" -w "$W.LSU_LROR-FLR2.amf.fasta" -r "$DELAVAUX_DB" \
    -o "$W.LSU_LROR-FLR2.amf.vs_delavaux" --threads "$THREADS" > /dev/null
fi
if [[ -n "${MAARJAM_DB:-}" && -s "$W.SSU_VT_NS31-AM1.amf.fasta" ]]; then
  python3 "$HERE/validate_windows.py" -w "$W.SSU_VT_NS31-AM1.amf.fasta" -r "$MAARJAM_DB" \
    -o "$W.SSU_VT_NS31-AM1.amf.vs_maarjam" --threads "$THREADS" > /dev/null
fi
echo "done: $W.records.tsv  (triage: $W.triage.tsv)"
