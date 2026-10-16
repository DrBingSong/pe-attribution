# -*- coding: utf-8 -*-
"""收尾回传：F1 图件 + scF 敏感性 + bulk 验证集（tar.gz 分片 6 线程）。"""
import math
import os
import tarfile
import threading
import time
from pathlib import Path

import paramiko

HOST, PORT, USER = "connect.cqa1.seetacloud.com", 33510, "root"
CHUNK = 70000000
FIG_LOCAL = Path(r"D:/zcode/z0927/pe_signature/results/figures")
SENS_LOCAL = Path(r"D:/zcode/z0927/pe_signature/results/sensitivity")
BULK_LOCAL = Path("F:/pe_signature/results/bulk")
N_THREADS = 6

SMALL = [
    ("/root/pe_signature/results/figures/fig1_annotation_umap.png", FIG_LOCAL / "fig1_annotation_umap.png"),
    ("/root/pe_signature/results/sensitivity/scf_immune_sig.csv", SENS_LOCAL / "scf_immune_sig.csv"),
]
for g in ["GSE173193", "GSE282038", "GSE298119", "GSE298602", "GSE240922", "GSE290578", "GSE329173"]:
    SMALL.append(("/root/pe_signature/results/figures/umap_%s.png" % g, FIG_LOCAL / ("umap_%s.png" % g)))

REMOTE_BULK_TGZ = "/root/pe_signature/data/bulk_all.tar.gz"


def connect(pw):
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    kh = Path(__file__).with_name(".ssh_known_hosts")
    kh.touch(exist_ok=True)
    c.load_host_keys(str(kh))
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, PORT, USER, pw, timeout=30)
    return c


def run_fresh(pw, cmd, timeout=300):
    c = connect(pw)
    try:
        _, o, e = c.exec_command(cmd, timeout=timeout)
        o.channel.recv_exit_status()
        return o.read().decode(errors="replace"), e.read().decode(errors="replace")
    finally:
        c.close()


def main():
    pw = os.environ["AUTODL_PW"]
    FIG_LOCAL.mkdir(parents=True, exist_ok=True)
    SENS_LOCAL.mkdir(parents=True, exist_ok=True)
    BULK_LOCAL.mkdir(parents=True, exist_ok=True)

    c = connect(pw)
    s = c.open_sftp()
    for rp, lp in SMALL:
        try:
            s.get(rp, str(lp))
            print("ok:", lp.name, lp.stat().st_size)
        except Exception as e:
            print("skip:", rp, type(e).__name__)
    s.close()
    c.close()

    out, err = run_fresh(pw, "cd /root/pe_signature/data && tar czf bulk_all.tar.gz bulk && stat -c%s bulk_all.tar.gz")
    size = int(out.strip().splitlines()[-1]) if out.strip() else 0
    print("bulk tar.gz:", size, "bytes")
    if size == 0:
        print("tar 失败:", err[:200])
        return
    n_parts = math.ceil(size / CHUNK)
    sub = "/root/pe_signature/data/bulk_split"
    run_fresh(pw, "rm -rf %s && mkdir -p %s && split -d -b %d %s %s/p_" % (sub, sub, CHUNK, REMOTE_BULK_TGZ, sub))

    def fetch(i):
        cc = connect(pw)
        try:
            ss = cc.open_sftp()
            ss.get("%s/p_%02d" % (sub, i), str(BULK_LOCAL / ("bulk.p%02d" % i)))
        finally:
            cc.close()

    threads = []
    for i in range(n_parts):
        t = threading.Thread(target=fetch, args=(i,))
        t.start()
        threads.append(t)
        time.sleep(1)
    for t in threads:
        t.join()
    total = sum((BULK_LOCAL / ("bulk.p%02d" % i)).stat().st_size for i in range(n_parts))
    print("bulk 分片:", n_parts, "总字节:", total, "OK" if total == size else "MISMATCH")
    run_fresh(pw, "rm -rf %s %s" % (sub, REMOTE_BULK_TGZ))
    print("FINAL-DONE")


if __name__ == "__main__":
    main()
