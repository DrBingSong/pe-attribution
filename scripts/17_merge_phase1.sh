#!/usr/bin/env bash
# Phase 1 分片合并 + 字节校验（分片名: <file>.p00, <file>.p01, ...）
set -u
cd /f/pe_signature/results/phase1

merge() {  # $1=文件名 $2=期望字节数 $3=分片数
  local f=$1 want=$2 n=$3
  rm -f "$f"
  for i in $(seq -f "%02g" 0 $((n-1))); do
    cat "$f.p$i" >> "$f"
  done
  local got=$(stat -c%s "$f")
  if [ "$got" -eq "$want" ]; then
    echo "OK  $f  $got"
    for i in $(seq -f "%02g" 0 $((n-1))); do rm -f "$f.p$i"; done
  else
    echo "FAIL $f  got=$got want=$want"
  fi
}

merge atts_GSE240922.npz 1294582205 19
merge atts_GSE290578.npz 1415714195 21
merge atts_GSE298602.npz  740855035 11
merge atts_GSE329173.npz  427054994  7
merge atts_GSE298119.npz  434321548  7
merge atts_GSE282038.npz  118081496  2

echo "--- 剩余散落分片（应为空） ---"
ls *.p[0-9][0-9] 2>/dev/null || echo none
echo "--- 最终文件 ---"
ls -la *.npz
