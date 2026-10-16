#!/usr/bin/env bash
# 用系统 miniconda 的 torch 重建 venv（--system-site-packages），装小依赖
set -u
pkill -9 -x pip 2>/dev/null
sleep 1
rm -rf /root/pe_signature/env
bash -lc 'python -m venv --system-site-packages /root/pe_signature/env'
P=/root/pe_signature/env/bin/pip
$P install -i https://mirrors.aliyun.com/pypi/simple/ scanpy sc-signature captum
/root/pe_signature/env/bin/python -c "import torch, scanpy; print('venv-torch', torch.__version__, torch.cuda.is_available())"
echo SMALL-PIP-OK
echo "PIP-ALL-DONE (system-torch-2.8.0+cu128)" >> /root/pe_signature/logs/pip2.log
