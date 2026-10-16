#!/usr/bin/env bash
# 本机 8 线程分段下载 GSE173193_RAW.tar（健壮版）
# 每片硬校验字节数，失败自动整片重试（最多 6 轮）；全部通过才合并。
set -u
cd /f/pe_signature/raw
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173193/suppl/GSE173193_RAW.tar"
SIZE=503316480
N=8
CHUNK=$((SIZE / N))   # 62914560，恰好整除
rm -f part_*.bin GSE173193_RAW.tar

for round in 1 2 3 4 5 6; do
  pids=()
  for i in $(seq 0 $((N-1))); do
    S=$((i * CHUNK)); E=$((S + CHUNK - 1)); [ $i -eq $((N-1)) ] && E=$((SIZE - 1))
    want=$CHUNK; [ $i -eq $((N-1)) ] && want=$((SIZE - S))
    have=0; [ -f part_$i.bin ] && have=$(stat -c%s part_$i.bin)
    if [ "$have" -ne "$want" ]; then
      ( curl -sf -r $S-$E -o part_$i.bin --retry 2 --max-time 2400 "$URL" ) & pids+=($!)
    fi
  done
  wait "${pids[@]}" 2>/dev/null
  bad=0
  for i in $(seq 0 $((N-1))); do
    want=$CHUNK; [ $i -eq $((N-1)) ] && want=$((SIZE - i*CHUNK))
    have=0; [ -f part_$i.bin ] && have=$(stat -c%s part_$i.bin)
    if [ "$have" -ne "$want" ]; then echo "round $round: part $i bad ($have/$want)"; bad=1; fi
  done
  if [ $bad -eq 0 ]; then
    cat part_0.bin part_1.bin part_2.bin part_3.bin part_4.bin part_5.bin part_6.bin part_7.bin > GSE173193_RAW.tar
    got=$(stat -c%s GSE173193_RAW.tar)
    if [ "$got" -eq "$SIZE" ]; then
      rm -f part_*.bin
      echo "ALL OK: $got bytes"; exit 0
    fi
    echo "merged size mismatch: $got"; exit 1
  fi
done
echo "FAILED after 6 rounds"; exit 1
