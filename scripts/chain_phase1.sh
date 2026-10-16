#!/usr/bin/env bash
# Phase 1 链式任务：等批量下载完成 → 批量预处理 → 全量归因
set -u
LOG=/root/pe_signature/logs
while ! grep -q "BATCH-DL-DONE" "$LOG/batch_dl.log" 2>/dev/null; do sleep 30; done
export MPLBACKEND=Agg
P=/root/pe_signature/env/bin/python
"$P" /root/pe_signature/scripts/12_preprocess_batch.py > "$LOG/prep12.log" 2>&1
"$P" /root/pe_signature/scripts/13_phase1_attribution.py > "$LOG/phase1.log" 2>&1
echo CHAIN1-DONE > "$LOG/chain1.done"
