# -*- coding: utf-8 -*-
"""
完整性修复 v5b：远端 split → 远端逐片 SHA256 → 分片拉取（逐片校验+重试）→ 报告。
合并与最终校验由 Bash 脚本（17_merge_phase1.sh 扩展版）完成。
已确认完好的文件自动跳过。
"""
import hashlib
import json
import math
import os
import threading
import time
from pathlib import Path

import paramiko

HOST, PORT, USER = "connect.cqa1.seetacloud.com", 33510, "root"
REMOTE_DIR = "/root/pe_signature/results/phase1"
SPLIT_DIR = "/root/autodl-tmp/pe_sig_data/split_tmp"  # 数据盘，避免系统盘写满
CHUNK = 70000000
LOCAL_DIR = Path("F:/pe_signature/results/phase1")
N_THREADS = 4

FILES = {
    "atts_GSE282038.npz": 118081496,
    "atts_GSE298119.npz": 434321548,
}  # 其余 4 个本地加载验证已通过
SKIP_OK = set()


def connect(pw):
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    kh = Path(__file__).with_name(".ssh_known_hosts")
    kh.touch(exist_ok=True)
    c.load_host_keys(str(kh))
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, PORT, USER, pw, timeout=30)
    return c


def run_fresh(pw, cmd, timeout=600):
    c = connect(pw)
    try:
        _, o, e = c.exec_command(cmd, timeout=timeout)
        o.channel.recv_exit_status()
        return o.read().decode(errors="replace"), e.read().decode(errors="replace")
    finally:
        c.close()


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def chunk_local(name, idx):
    p = (LOCAL_DIR / (name + ".sh%02d" % idx)).resolve()
    if LOCAL_DIR.resolve() not in p.parents:
        raise ValueError("path escape")
    return p


def parse_remote(pw, sub):
    txt, _ = run_fresh(pw, "cd %s/%s && ls p_* | while read f; do printf '%s %s\\n' \"$f\" \"$(sha256sum \"$f\" | cut -d' ' -f1)\"; done"
                       % (SPLIT_DIR, sub, "%s", "%s"))
    out = {}
    for line in txt.strip().splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].startswith("p_"):
            out[SPLIT_DIR + "/" + sub + "/" + parts[0]] = parts[1]
    return out


def fetch_verified(pw, remote, local, want_sha, tries=3):
    for t in range(tries):
        try:
            c = connect(pw)
            try:
                s = c.open_sftp()
                s.get(remote, str(local))
            finally:
                c.close()
            if local.exists() and sha256_of(local) == want_sha:
                return True
        except Exception as e:
            print("    retry%d %s: %s" % (t + 1, local.name, type(e).__name__), flush=True)
        time.sleep(4)
    return local.exists() and sha256_of(local) == want_sha


def repair(name, expect, pw):
    sub = name + ".split"
    n_parts = int(math.ceil(expect / CHUNK))
    print("== %s（%d 片）==" % (name, n_parts), flush=True)

    remote_hashes = {}
    hashfile = "%s/%s/hashes.txt" % (SPLIT_DIR, sub)
    for attempt in range(3):
        run_fresh(pw, "rm -rf %s/%s && mkdir -p %s/%s && split -d -b %d %s/%s %s/%s/p_"
                  % (SPLIT_DIR, sub, SPLIT_DIR, sub, CHUNK, REMOTE_DIR, name, SPLIT_DIR, sub))
        run_fresh(pw, "cd %s/%s && sha256sum p_* > %s" % (SPLIT_DIR, sub, hashfile))
        # 哈希清单经 SFTP 读取（exec stdout 在中继上不可靠）
        local_hashlist = LOCAL_DIR / ("hashes_" + name + ".txt")
        try:
            c = connect(pw)
            s = c.open_sftp()
            s.get(hashfile, str(local_hashlist))
            s.close()
            c.close()
            htxt = local_hashlist.read_text(encoding="utf-8", errors="replace")
        except Exception as e:
            print("  hashes.txt 下载失败:", type(e).__name__, flush=True)
            htxt = "" 
        for line in htxt.strip().splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1].startswith("p_"):
                remote_hashes[SPLIT_DIR + "/" + sub + "/" + parts[1]] = parts[0]
        if len(remote_hashes) == n_parts:
            break
        print("  split/sha 第 %d 次尝试不完整（%d/%d 片），重试..." % (attempt + 1, len(remote_hashes), n_parts), flush=True)
        time.sleep(8)
    if len(remote_hashes) != n_parts:
        print("  远端分片数异常 %d != %d" % (len(remote_hashes), n_parts), flush=True)
        return False

    jobs = []
    for i in range(n_parts):
        rp = "%s/%s/p_%02d" % (SPLIT_DIR, sub, i)
        lpath = chunk_local(name, i)
        jobs.append((i, rp, lpath, remote_hashes[rp]))

    ok = {j[0]: False for j in jobs}
    todo_queue = list(jobs)
    lock = threading.Lock()

    def worker():
        while True:
            with lock:
                if not todo_queue:
                    return
                idx, rp, lpath, sha = todo_queue.pop(0)
            if fetch_verified(pw, rp, lpath, sha):
                ok[idx] = True
                print("  片 %02d ✓" % idx, flush=True)

    ts = [threading.Thread(target=worker) for _ in range(N_THREADS)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()

    bad = [j[0] for j in jobs if not ok[j[0]]]
    if bad:
        print("  %s 坏片: %s" % (name, bad), flush=True)
        return False
    print("  %s 全部 %d 片 SHA256 校验通过" % (name, n_parts), flush=True)
    return True


def main():
    pw = os.environ["AUTODL_PW"]
    results = {}
    for name, expect in FILES.items():
        lp = LOCAL_DIR / name
        if name in SKIP_OK and lp.exists():
            print("跳过（已验证）:", name, flush=True)
            results[name] = "skip-ok"
            continue
        results[name] = "ok" if repair(name, expect, pw) else "FAILED"
    print("REPAIR-SUMMARY:", json.dumps(results), flush=True)


if __name__ == "__main__":
    main()
