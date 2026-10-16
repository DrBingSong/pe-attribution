#!/usr/bin/env bash
# 链式任务：等 pip 完成 → 实例跑预处理 → 跑 Phase 0 归因
set -u
LOG=/root/pe_signature/logs
while ! grep -q "PIP-ALL-DONE" "$LOG/pip2.log" 2>/dev/null; do sleep 20; done
export MPLBACKEND=Agg
P=/root/pe_signature/env/bin/python
cd /root/pe_signature
sed -e 's|F:/pe_signature/raw/extracted|/root/pe_signature/data/raw/extracted|' \
    -e 's|F:/pe_signature/processed|/root/pe_signature/data/processed|' \
    -e 's|D:/zcode/z0927/pe_signature/results|/root/pe_signature/results|' \
    scripts/01_preprocess_gse173193.py > scripts/01_preprocess_instance.py
"$P" scripts/01_preprocess_instance.py > "$LOG/preprocess.log" 2>&1
"$P" scripts/02_phase0_attribution.py > "$LOG/phase0.log" 2>&1
echo CHAIN-DONE > "$LOG/chain.done"
