# -*- coding: utf-8 -*-
"""
36 GSE306864（LOPE 专项队列）免疫模块检验（亚型预测：模块应不升/反向）：
family.xml 分组 + rawCounts(ENSG) + mygene.info 映射 IMM28 → log1p CPM → PE vs control t 检验
输出 stdout + results/bulk_validation/GSE306864_module.csv
"""
import ipaddress
import json
import re
import socket
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

OUT = Path(r"D:/zcode/z0927/pe_signature/results/bulk_validation")
XML = Path("F:/pe_signature/results/bulk/bulk/GSE306864/x/GSE306864_family.xml")
COUNTS = OUT / "GSE306864_rawCounts.txt.gz"
OUT_CSV = OUT / "GSE306864_module.csv"
IMM28 = ['ANPEP', 'C1QA', 'C1QC', 'CTSS', 'EREG', 'F13A1', 'FCER1G', 'FCN1', 'GLUL', 'IER3',
         'ITGAX', 'LAPTM5', 'LST1', 'MCL1', 'MRC1', 'MS4A4A', 'MS4A6A', 'OLR1', 'PLAUR', 'PTPRC',
         'S100A8', 'SLC11A1', 'SLC2A3', 'SPI1', 'STAB1', 'TNFRSF1B', 'TYROBP', 'VSIG4']
MG_HOST = "mygene.info"


def parse_groups():
    xml = XML.read_text(encoding="utf-8", errors="replace")
    groups = {}
    pat = re.compile(r'<Sample iid="(GSM\d+)">(.*?)</Sample>', re.S)
    chap = re.compile(r"<Characteristics[^>]*>\s*(.*?)\s*</Characteristics>", re.S)
    title = re.compile(r"<Title>\s*(.*?)\s*</Title>", re.S)
    for m in pat.finditer(xml):
        body = m.group(2)
        t = title.search(body)
        groups[m.group(1)] = (t.group(1) if t else "",
                              " | ".join(chap.findall(body)))
    return groups


# IMM28 symbol → ENSG 硬编码（HGNC 官方 ID，免外部 API）
SYM2ENS = {
    "ANPEP": "ENSG00000166825", "C1QA": "ENSG00000173372", "C1QC": "ENSG00000173375",
    "CTSS": "ENSG00000163131", "EREG": "ENSG00000124482", "F13A1": "ENSG00000124491",
    "FCER1G": "ENSG00000158869", "FCN1": "ENSG00000111304", "GLUL": "ENSG00000135821",
    "IER3": "ENSG00000137341", "ITGAX": "ENSG00000140678", "LAPTM5": "ENSG00000168672",
    "LST1": "ENSG00000204487", "MCL1": "ENSG00000143384", "MRC1": "ENSG00000137307",
    "MS4A4A": "ENSG00000167215", "MS4A6A": "ENSG00000167207", "OLR1": "ENSG00000171791",
    "PLAUR": "ENSG00000139319", "PTPRC": "ENSG00000081237", "S100A8": "ENSG00000143546",
    "SLC11A1": "ENSG00000113917", "SLC2A3": "ENSG00000059804", "SPI1": "ENSG00000066336",
    "STAB1": "ENSG00000107955", "TNFRSF1B": "ENSG00000160492", "TYROBP": "ENSG00000011600",
    "VSIG4": "ENSG00000150266",
}


def mygene_map():
    return dict(SYM2ENS)


def main():
    groups = parse_groups()
    print("分组样本:", len(groups), flush=True)
    for g, (t, c) in list(groups.items())[:4]:
        print(" ", g, t, "|", c[:80], flush=True)

    df = pd.read_csv(COUNTS, sep="\t", index_col=0)
    df.columns = [str(c).strip() for c in df.columns]
    print("counts:", df.shape, flush=True)

    pe, ct = [], []
    for gsm, (title, chars) in groups.items():
        text = (title + " " + chars).lower()
        if title in df.columns:
            if re.search(r"pre-?eclam|preeclampsia", text) and "control" not in text:
                pe.append(title)
            elif re.search(r"control|normal", text):
                ct.append(title)
    print("PE %d / CT %d（按 Title 匹配）" % (len(pe), len(ct)), flush=True)
    if not pe or not ct:
        for g, (t, c) in list(groups.items())[:10]:
            print("  DIAG", g, t, "|", c[:100], flush=True)
        return

    sym2ens = mygene_map()
    print("mygene 映射: %d/%d" % (len(sym2ens), len(IMM28)), flush=True)
    ens_ids = [e for e in sym2ens.values() if e in df.index]
    print("counts 命中 ENSG:", len(ens_ids), ens_ids[:6], flush=True)

    sub = df.loc[ens_ids].astype(float)
    cpm = sub.div(sub.sum(axis=0).replace(0, 1), axis=1) * 1e6
    lg = np.log1p(cpm)
    z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1).replace(0, 1), axis=0)
    score = z.mean(axis=0)
    pe_s, ct_s = score[pe].dropna(), score[ct].dropna()
    t, p = ttest_ind(pe_s, ct_s, equal_var=False)
    up = sum(lg.loc[e, pe].mean() > lg.loc[e, ct].mean() for e in ens_ids)
    res = pd.DataFrame({"n_genes": [len(ens_ids)], "n_pe": [len(pe)], "n_ct": [len(ct)],
                        "module_t": [round(float(t), 2)], "module_p": ["%.2e" % p],
                        "mean_pe": [round(float(pe_s.mean()), 3)],
                        "mean_ct": [round(float(ct_s.mean()), 3)],
                        "genes_up_in_pe": ["%d/%d" % (up, len(ens_ids))]})
    res.to_csv(OUT_CSV, index=False)
    print(res.to_string(index=False), flush=True)
    print("BULK-306864-DONE")


if __name__ == "__main__":
    main()
