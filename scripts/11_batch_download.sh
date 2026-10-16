#!/usr/bin/env bash
# Phase 1 批量下载 PE scRNA-seq 数据集（大文件自动多线程）
set -u
BASE=/root/pe_signature/data
mkdir -p $BASE/gse
ACC="GSE282038 GSE298119 GSE298602 GSE240922 GSE290578 GSE329173"

for GSE in $ACC; do
  PRE=${GSE:0:6}nnn          # GSE282038 -> GSE282nnn
  DIR=$BASE/gse/$GSE
  mkdir -p "$DIR"
  [ -f "$DIR/.done" ] && { echo "skip $GSE (done)"; continue; }
  URLROOT="https://ftp.ncbi.nlm.nih.gov/geo/series/${PRE}/${GSE}/suppl/"
  echo "=== $GSE ==="
  LIST=$(curl -sf --retry 3 --max-time 60 "${URLROOT}/" | grep -oE 'href="[^"]+\.(tar|gz|tsv|h5|h5ad|mtx)"' | sed 's/href="//;s/"//' | grep -v '^/' | sort -u)
  if [ -z "$LIST" ]; then echo "$GSE: no files listed"; echo 1 > "$DIR/.done"; continue; fi
  for F in $LIST; do
    T="$DIR/$F"
    [ -f "$T" ] && continue
    SZ=$(curl -sI --retry 3 --max-time 30 "https://ftp.ncbi.nlm.nih.gov/geo/series/${PRE}/${GSE}/suppl/$F" | grep -i content-length | tail -1 | tr -dc 0-9)
    SZ=${SZ:-0}
    if [ "$SZ" -gt 150000000 ]; then
      N=8; CSZ=$((SZ / N))
      for i in $(seq 0 $((N-1))); do
        S=$((i * CSZ)); E=$((S + CSZ - 1)); [ $i -eq $((N-1)) ] && E=$((SZ - 1))
        nohup curl -sf -r $S-$E -o "$T.part$i" --retry 3 "${URLROOT}/$F" >/dev/null 2>&1 &
      done
      wait
      cat "$T".part* > "$T" && rm -f "$T".part*
    fi
    if [ ! -f "$T" ] || [ "$(stat -c%s "$T")" -ne "$SZ" ]; then
      curl -sf --retry 3 -o "$T" "${URLROOT}/$F" || echo "  FAIL $F"
    fi
    echo "  $F  $(stat -c%s "$T" 2>/dev/null || echo 0) bytes"
  done
  echo 1 > "$DIR/.done"
  echo "$GSE done: $(du -sm $DIR | cut -f1) MB"
done
echo BATCH-DL-DONE
