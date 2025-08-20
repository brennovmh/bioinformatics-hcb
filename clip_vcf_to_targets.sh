#!/usr/bin/env bash
set -euo pipefail

# Clip CNV VCF endpoints to the convex hull of on-target regions in a BED.
# Usage: ./clip_vcf_to_targets.sh INPUT.vcf PANEL.bed OUTPUT.vcf
# Requires: awk, sort, bedtools
# Optional: KEEP_OFFTARGET=1 to keep variants with zero on-target overlap (unchanged)

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 INPUT.vcf PANEL.bed OUTPUT.vcf" >&2
  exit 1
fi

IN_VCF="$1"
IN_BED="$2"
OUT_VCF="$3"
KEEP="${KEEP_OFFTARGET:-0}"   # 0 = drop 100% off-target SVs; 1 = keep them unchanged

workdir="$(dirname "$OUT_VCF")"
mkdir -p "$workdir"

SV_BED6="$workdir/.sv.bed6"
SV_BED6_SORT="$workdir/.sv.sorted.bed6"
BED3_CLEAN="$workdir/.panel.clean3.bed"
BED3_SORT="$workdir/.panel.sorted.bed"
MAP_OUT="$workdir/.sv.map"
CLIP_TSV="$workdir/.sv.clip.tsv"
AUDIT_TSV="${OUT_VCF%.vcf}.clip_audit.tsv"

# 1) Clean panel BED -> 3 columns, tabs, no headers
sed -E 's/\r$//; s/[[:space:]]+/\t/g; s/\t+$//' "$IN_BED" \
 | grep -v -E '^(track|browser|#)' \
 | awk 'BEGIN{OFS="\t"} NF>=3 && $3>$2 {print $1,$2,$3}' \
 > "$BED3_CLEAN"

# 2) Detect chrom style in VCF (chr10 vs 10) from first non-header line
VCF_CHR="$(awk '!/^#/ {print $1; exit}' "$IN_VCF" 2>/dev/null || true)"
if [[ -z "${VCF_CHR:-}" ]]; then
  echo "[WARN] VCF seems empty of variants; copying as-is." >&2
  cp "$IN_VCF" "$OUT_VCF"
  exit 0
fi
if [[ "$VCF_CHR" == chr* ]]; then
  awk 'BEGIN{OFS="\t"} {if($1 !~ /^chr/) $1="chr"$1; print}' "$BED3_CLEAN" > "$BED3_SORT.tmp"
else
  awk 'BEGIN{OFS="\t"} {sub(/^chr/,"",$1); print}' "$BED3_CLEAN" > "$BED3_SORT.tmp"
fi

# 3) Sort panel BED
LC_ALL=C sort -k1,1 -k2,2n "$BED3_SORT.tmp" > "$BED3_SORT"
rm -f "$BED3_SORT.tmp"

# 4) Build BED6 of SVs from VCF body: chr start0 end SVTYPE SVLEN rid
awk -F'\t' 'BEGIN{OFS="\t"; rid=0}
  /^#/ {next}
  {
    if (NF<8) next
    rid++
    # extract END, SVTYPE, SVLEN from INFO
    split($8,info,";"); end=""; svt=""; svlen=""
    for(i=1;i<=length(info);i++){
      split(info[i],kv,"=")
      if(kv[1]=="END") end=kv[2]
      else if(kv[1]=="SVTYPE") svt=kv[2]
      else if(kv[1]=="SVLEN") svlen=kv[2]
    }
    if(end=="") end=$2
    if($2!~/^[0-9]+$/ || end!~/^[0-9]+$/) next
    print $1, $2-1, end, svt, svlen, rid
  }' "$IN_VCF" > "$SV_BED6"

# 5) Sort SVs for bedtools
LC_ALL=C sort -k1,1 -k2,2n "$SV_BED6" > "$SV_BED6_SORT"

# 6) bedtools map: for each SV, get min(target.start), max(target.end)
bedtools map -a "$SV_BED6_SORT" -b "$BED3_SORT" -c 2,3 -o min,max > "$MAP_OUT"

# 7) Clip to on-target convex hull and compute new POS/END/SVLEN
#    Output: rid  newPOS  newEND  newSVLEN   oldPOS  oldEND  oldSVLEN  SVTYPE
awk 'BEGIN{OFS="\t"}
  {
    chr=$1; s0=$2; e=$3; svt=$4; svlen=$5; rid=$6; minb=$7; maxb=$8
    if(minb==-1 || maxb==-1) next   # no overlap -> drop here; optional re-add later
    s = s0 + 1                      # to 1-based
    cs = (s>minb ? s : minb)        # clipped start
    ce = (e<maxb ? e : maxb)        # clipped end
    if (ce<=cs) next
    newlen = ce - cs
    if     (svt=="DEL") newsvlen = -newlen
    else if(svt=="DUP") newsvlen =  newlen
    else                newsvlen =  newlen
    print rid, cs, ce, newsvlen, s, e, svlen, svt
  }' "$MAP_OUT" > "$CLIP_TSV"

# 8) Recompose VCF (headers preserved)
#    - if KEEP=0, SVs without on-target overlap are dropped
#    - if KEEP=1, they are kept unchanged
awk -F'\t' -v OFS="\t" -v KEEP="$KEEP" '
  NR==FNR {
    rid=$1; pos[rid]=$2; end[rid]=$3; svlen[rid]=$4; next
  }
  /^#/ {print; next}
  {
    rec=++k
    if(rec in pos){
      $2 = pos[rec]
      # update END
      if($8 ~ /(^|;)END=/) sub(/END=[^;]*/, "END=" end[rec], $8)
      else $8=$8";END="end[rec]
      # update SVLEN
      if($8 ~ /(^|;)SVLEN=/) sub(/SVLEN=[^;]*/, "SVLEN=" svlen[rec], $8)
      else $8=$8";SVLEN="svlen[rec]
      print
    } else {
      if (KEEP) print
    }
  }' "$CLIP_TSV" "$IN_VCF" > "$OUT_VCF"

# 9) Write audit table
{
  echo -e "RID\tSVTYPE\tOLD_POS\tOLD_END\tOLD_SVLEN\tNEW_POS\tNEW_END\tNEW_SVLEN"
  awk 'BEGIN{OFS="\t"} {print $1,$8,$5,$6,$7,$2,$3,$4}' "$CLIP_TSV" | LC_ALL=C sort -n
} > "$AUDIT_TSV"

echo "[OK] Wrote clipped VCF: $OUT_VCF"
echo "[OK] Audit TSV: $AUDIT_TSV"
