# -*- coding: utf-8 -*-
"""
27 GO/KEGG 富集（Enrichr REST API）：四细胞类型的跨队列共识 top200 交集基因。
安全：仅允许 https + maayanlab.cloud，解析 IP 阻断私网/环回。
库：GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022
输出：results/enrichr/<type>_<lib>.csv（top15）+ stdout 摘要。
"""
import ipaddress
import json
import socket
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

SRC = Path(r"D:/zcode/z0927/pe_signature/results/headtohead")
OUT = Path(r"D:/zcode/z0927/pe_signature/results/enrichr")
OUT.mkdir(parents=True, exist_ok=True)

TYPES = ["Trophoblast", "Immune", "Endothelial", "Stromal"]
LIBS = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]
BASE = "https://maayanlab.cloud/Enrichr"
ALLOWED_HOST = "maayanlab.cloud"
HDR = {"User-Agent": "pe-signature-analysis/1.0"}


def check_url(url):
    p = urllib.parse.urlparse(url)
    assert p.scheme == "https" and p.hostname == ALLOWED_HOST, "host not allowed"
    ip = socket.getaddrinfo(p.hostname, 443)[0][4][0]
    addr = ipaddress.ip_address(ip)
    assert not (addr.is_private or addr.is_loopback or addr.is_link_local), "private ip"


def post(url, data, headers):
    check_url(url)
    import requests
    resp = requests.post(url, files=data, timeout=60)
    resp.raise_for_status()
    return resp.json()


def get(url):
    check_url(url)
    req = urllib.request.Request(url, headers=HDR)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def enrich(genes, lib):
    resp = post(BASE + "/addList",
                {"list": (None, "\n".join(genes)), "description": (None, "pe_consensus")}, HDR)
    user = resp["userListId"]
    for attempt in range(3):
        try:
            r = get("%s/enrich?userListId=%s&backgroundType=%s"
                    % (BASE, user, urllib.parse.quote(lib)))
            return r.get(lib, [])
        except Exception:
            time.sleep(4)
    return []


def main():
    summary = {}
    for t in TYPES:
        genes = pd.read_csv(SRC / ("%s_attr_consensus_top200_intersect.csv" % t)).gene.tolist()
        print("== %s（%d 基因）==" % (t, len(genes)), flush=True)
        summary[t] = {}
        for lib in LIBS:
            rows = enrich(genes, lib)
            if not rows:
                print("  [%s] 无结果/失败" % lib, flush=True)
                continue
            df = pd.DataFrame(rows)
            df = df.iloc[:, :8]
            df.columns = ["rank", "term", "pval", "zscore", "combined",
                          "overlap_genes", "adj_pval", "old_pval"]
            df = df[df.adj_pval < 0.05].head(15)
            df.to_csv(OUT / ("%s_%s.csv" % (t, lib)), index=False)
            top3 = [(r.term[:60], "%.1e" % r.adj_pval) for r in df.head(3).itertuples()]
            summary[t][lib] = top3
            print("  [%s] %d 项 adj.p<0.05 | top: %s" % (lib, len(df), top3), flush=True)
            time.sleep(1.5)
    print("ENRICH-JSON-BEGIN")
    print(json.dumps(summary, indent=1, ensure_ascii=False))
    print("ENRICH-JSON-END")


if __name__ == "__main__":
    main()
