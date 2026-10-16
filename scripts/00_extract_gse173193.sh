#!/usr/bin/env bash
# 解压 GSE173193_RAW.tar 并重组为 scanpy read_10x_mtx 需要的目录结构
# 输出: F:/pe_signature/raw/extracted/<GSMxxxxx_C1>/{barcodes,features,matrix}.{tsv,mtx}.gz
set -e
cd /f/pe_signature/raw
mkdir -p staging extracted
tar -xf GSE173193_RAW.tar -C staging
for f in staging/*_barcodes.tsv.gz; do
  base=$(basename "$f" _barcodes.tsv.gz)
  mkdir -p "extracted/$base"
  mv "staging/${base}_barcodes.tsv.gz" "extracted/$base/barcodes.tsv.gz"
  mv "staging/${base}_features.tsv.gz" "extracted/$base/features.tsv.gz"
  mv "staging/${base}_matrix.mtx.gz"   "extracted/$base/matrix.mtx.gz"
  echo "OK $base"
done
rmdir staging 2>/dev/null || true
echo "--- extracted ---"
ls extracted
