# -*- coding: utf-8 -*-
"""
35 bulk 验证（GSE25906，Illumina GPL6102，本地 family 包）：
GPL6102-tbl-1.txt 建 ILMN→symbol 映射；family.xml 取分组；
每样本 tbl-1.txt 信号列合并 → 基因×样本 → IMM28 模块 PE vs control。
输出 stdout + results/bulk_validation/GSE25906_module.csv
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

SRC = Path("F:/pe_signature/results/bulk/bulk/GSE25906/x")
OUT = Path(r"D:/zcode/z0927/pe_signature/results/bulk_validation")
OUT.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT / "GSE25906_module.csv"
IMM28 = ['ANPEP', 'C1QA', 'C1QC', 'CTSS', 'EREG', 'F13A1', 'FCER1G', 'FCN1', 'GLUL', 'IER3',
         'ITGAX', 'LAPTM5', 'LST1', 'MCL1', 'MRC1', 'MS4A4A', 'MS4A6A', 'OLR1', 'PLAUR', 'PTPRC',
         'S100A8', 'SLC11A1', 'SLC2A3', 'SPI1', 'STAB1', 'TNFRSF1B', 'TYROBP', 'VSIG4']


def build_map():
    """GPL6102-tbl-1.txt: 无表头，列0=ILMN 探针，列4=Symbol"""
    f = SRC / "GPL6102-tbl-1.txt"
    mapping = {}
    for ln in f.read_text(encoding="utf-8", errors="replace").splitlines():
        cols = ln.split("\t")
        if len(cols) > 5 and cols[0].startswith("ILMN_"):
            sym = cols[5].strip().strip('"')
            if sym and not sym.startswith(("HS.", "LOC", "XM_", "AK")):
                mapping[cols[0]] = sym
    return mapping


def parse_groups():
    """MINiML family.xml: Sample iid=GSM，Characteristics 内含 classification"""
    xml = (SRC / "GSE25906_family.xml").read_text(encoding="utf-8", errors="replace")
    groups = {}
    pat = re.compile(r'<Sample iid="(GSM\d+)">(.*?)</Sample>', re.S)
    chap = re.compile(r"<Characteristics[^>]*>\s*(.*?)\s*</Characteristics>", re.S)
    for m in pat.finditer(xml):
        groups[m.group(1)] = " | ".join(chap.findall(m.group(2)))
    return groups


def load_expr(mapping, gsms):
    cols = {}
    for gsm in gsms:
        f = SRC / ("%s-tbl-1.txt" % gsm)
        if not f.exists():
            continue
        d = {}
        for ln in f.read_text(encoding="utf-8", errors="replace").splitlines():
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            p, val = parts[0], parts[1]
            if p in mapping:
                try:
                    d[mapping[p]] = float(val)
                except ValueError:
                    pass
        cols[gsm] = d
    df = pd.DataFrame(cols)
    return df.groupby(level=0).mean().T  # 样本 × 基因


def main():
    mapping = build_map()
    print("ILMN→symbol 映射: %d 探针" % len(mapping), flush=True)
    groups = parse_groups()
    print("样本分组解析: %d GSM" % len(groups), flush=True)
    for g, c in list(groups.items())[:3]:
        print("  ", g, c[:90], flush=True)

    gsms = sorted(groups)
    df = load_expr(mapping, gsms)
    print("表达矩阵: %d 样本 × %d 基因" % df.shape, flush=True)

    pe = [g for g in df.index if re.search(r"pre-?eclam", groups.get(g, ""), re.I)]
    ct = [g for g in df.index if re.search(r"control|normal", groups.get(g, ""), re.I)]
    print("PE %d / CT %d" % (len(pe), len(ct)), flush=True)
    if not pe or not ct:
        print("分组不足，退出")
        return

    genes = [g for g in IMM28 if g in df.columns]
    print("IMM28 命中: %d/%d →" % (len(genes), len(IMM28)), genes, flush=True)
    sub = df[genes]
    z = sub.sub(sub.mean(axis=0), axis=1).div(sub.std(axis=0).replace(0, 1), axis=1)
    score = z.mean(axis=1)
    t, p = ttest_ind(score[pe], score[ct], equal_var=False)
    up = sum(sub.loc[pe, g].mean() > sub.loc[ct, g].mean() for g in genes)
    res = pd.DataFrame({"n_genes": [len(genes)], "n_pe": [len(pe)], "n_ct": [len(ct)],
                        "module_t": [round(float(t), 2)], "module_p": ["%.2e" % p],
                        "mean_pe": [round(float(score[pe].mean()), 3)],
                        "mean_ct": [round(float(score[ct].mean()), 3)],
                        "genes_up_in_pe": ["%d/%d" % (up, len(genes))]})
    res.to_csv(OUT_CSV, index=False)
    print(res.to_string(index=False), flush=True)
    print("BULK-GSE25906-DONE")


if __name__ == "__main__":
    main()
