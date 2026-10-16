# -*- coding: utf-8 -*-
"""Phase 1 结果批量回传 v2：逐文件 split + 6 线程 s.get 分片，本地合并交由 Bash 完成。"""
import math
import os
import threading
import time
from pathlib import Path

import paramiko

HOST, PORT, USER = "connect.cqa1.seetacloud.com", 33510, "root"
REMOTE_DIR = "/root/pe_signature/results/phase1"
SPLIT_DIR = "/root/pe_signature/tmp_split"
CHUNK = 70000000
LOCAL_DIR = Path("F:/pe_signature/results/phase1")

NAMES = [
    "atts_GSE240922.npz", "atts_GSE290578.npz", "atts_GSE298602.npz",
    "atts_GSE329173.npz", "atts_GSE298119.npz", "atts_GSE282038.npz",
    "atts_scimilarity_full.npz", "meta_GSE240922.npz", "meta_GSE282038.npz",
    "meta_GSE298119.npz", "meta_GSE298602.npz", "meta_GSE329173.npz",
    "sanity_all.json.txt",
]
N_THREADS = 6


def connect(pw):
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    kh = Path(__file__).with_name(".ssh_known_hosts")
    kh.touch(exist_ok=True)
    c.load_host_keys(str(kh))
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, PORT, USER, pw, timeout=30)
    return c


def run(c, cmd, timeout=180):
    _, o, e = c.exec_command(cmd, timeout=timeout)
    o.channel.recv_exit_status()
    return o.read().decode(errors="replace")


def chunk_local(name, idx):
    """字面量前缀 + 受控序号的分片路径（contained in LOCAL_DIR）。"""
    if idx < 0 or idx > 63:
        raise ValueError("bad chunk idx")
    p = (LOCAL_DIR / (name + ".p%02d" % idx)).resolve()
    if p.parent != LOCAL_DIR.resolve():
        raise ValueError("path escape")
    return p


def fetch_one(pw, remote, local):
    c = connect(pw)
    try:
        s = c.open_sftp()
        s.get(remote, str(local))
    finally:
        c.close()


def main():
    pw = os.environ["AUTODL_PW"]
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)
    c = connect(pw)

    sizes = {}
    for name in NAMES:
        rfile = REMOTE_DIR + "/" + name
        sizes[name] = -1
        for attempt in range(4):
            try:
                out = run(c, "stat -c%s " + rfile).strip()
                sizes[name] = int(out)
                break
            except ValueError:
                print("stat 空返回，重试", name, "第", attempt + 1, "次", flush=True)
                time.sleep(8)
                c = connect(pw)
        if sizes[name] < 0:
            print("跳过（远端 stat 失败）:", name, flush=True)

    t0 = time.time()
    for name in NAMES:
        expect = sizes[name]
        if expect < 0:
            continue
        lp = LOCAL_DIR / name
        if lp.exists() and lp.stat().st_size == expect:
            print("skip(已存在):", name, flush=True)
            continue
        n_parts = max(1, math.ceil(expect / CHUNK))
        parts_ok = all(chunk_local(name, i).exists() for i in range(n_parts))
        if parts_ok:
            print("  分片已存在，跳过下载:", name, flush=True)
        elif n_parts == 1:
            fetch_one(pw, "%s/%s" % (REMOTE_DIR, name), lp)
            print("  单文件完成:", name, expect, flush=True)
        else:
            run(c, "rm -rf %s/%s && mkdir -p %s/%s && split -d -b %d %s/%s %s/%s/p_"
                % (SPLIT_DIR, name, SPLIT_DIR, name, CHUNK, REMOTE_DIR, name, SPLIT_DIR, name), timeout=240)
            threads = []
            for i in range(n_parts):
                rp = "%s/%s/p_%02d" % (SPLIT_DIR, name, i)
                lpath = chunk_local(name, i)
                th = threading.Thread(target=fetch_one, args=(pw, rp, lpath))
                th.start()
                threads.append(th)
                time.sleep(1)
            for th in threads:
                th.join()
            print("  分片回传完成:", name, n_parts, "片", flush=True)
        if n_parts > 1:
            total = sum(chunk_local(name, i).stat().st_size for i in range(n_parts))
            if total != expect:
                print("CHUNK-SIZE-MISMATCH:", name, total, "!=", expect, flush=True)
            else:
                print("CHUNKS-OK:", name, total, flush=True)
    c.close()
    print("PYTHON-DL-DONE, 请用 Bash 合并校验（总耗时 %.0fs）" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
