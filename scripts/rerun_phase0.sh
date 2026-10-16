#!/usr/bin/env bash
# 重跑 Phase 0：先补装 leidenalg，再跑归因（日志在 logs/，完成标记 chain.done）
set -u
/root/pe_signature/env/bin/pip install -i https://mirrors.aliyun.com/pypi/simple/ leidenalg > /root/pe_signature/logs/leiden.log 2>&1
rm -f /root/pe_signature/logs/chain.done /root/pe_signature/logs/phase0.log
MPLBACKEND=Agg /root/pe_signature/env/bin/python /root/pe_signature/scripts/02_phase0_attribution.py > /root/pe_signature/logs/phase0.log 2>&1
echo RERUN-DONE > /root/pe_signature/logs/chain.done
