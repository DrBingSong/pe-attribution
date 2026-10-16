# -*- coding: utf-8 -*-
"""v4 终极稳健版：每条命令独立连接（规避长连接状态污染），分片逐片校验+重试。"""
import hashlib
import math
import os
import threading
import time
from pathlib import Path

import paramiko

HOST, PORT, USER = "connect.cqa1.seetacloud.com", 33510, "root"
REMOTE_DIR = "/root/pe_signature/results/phase1"
SPLIT_DIR = "/root/pe_signature/tmp_split2"
CHUNK = 70000000
LOCAL_DIR = Path("F:/pe_signature/results/phase1")
FILES = {"atts_GSE240922.npz": 1294582205, "atts_GSE290578.npz": 1415714195}
N_THREADS = 4


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


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def chunk_local(name, idx):
    if idx < 0 or idx > 63:
        raise ValueError("bad idx")
    p = (LOCAL_DIR / (name + ".rp%02d" % idx)).resolve()
    if p.parent != LOCAL_DIR.resolve():
        raise ValueError("path escape")
    return p


def parse_remote_chunks(pw, split_sub):
    out = {}
    txt, _ = run_fresh(pw, "ls -la %s/%s 2>/dev/null" % (SPLIT_DIR, split_sub))
    for line in txt.strip().splitlines():
        parts = line.split()
        if len(parts) >= 9 and parts[0].startswith("-"):
            base = parts[8]
            if base.startswith("p_"):
                out[SPLIT_DIR + "/" + split_sub + "/" + base] = int(parts[4])
    return out


def fetch_with_retry(pw, remote, local, want_size, tries=3):
    for t in range(tries):
        try:
            c = connect(pw)
            try:
                s = c.open_sftp()
                s.get(remote, str(local))
            finally:
                c.close()
            if local.exists() and local.stat().st_size == want_size:
                return True
        except Exception as e:
            print("    retry%d %s: %s" % (t + 1, local.name, type(e).__name__), flush=True)
        time.sleep(5)
    return local.exists() and local.stat().st_size == want_size


def main():
    pw = os.environ["AUTODL_PW"]
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)

    for name, expect in FILES.items():
        print("== %s (%d bytes) ==" % (name, expect), flush=True)
        sub = name + ".split"
        n_parts = int(math.ceil(expect / CHUNK))
        rchunks = parse_remote_chunks(pw, sub)
        if not (len(rchunks) == n_parts and sum(rchunks.values()) == expect):
            err = run_fresh(pw, "rm -rf %s/%s && mkdir -p %s/%s && split -d -b %d %s/%s %s/%s/p_"
                            % (SPLIT_DIR, sub, SPLIT_DIR, sub, CHUNK, REMOTE_DIR, name, SPLIT_DIR, sub))
            if err.strip():
                print("  split stderr:", err.strip()[:200], flush=True)
            rchunks = parse_remote_chunks(pw, sub)
        if len(rchunks) != n_parts or sum(rchunks.values()) != expect:
            print("  远端分片异常（%d 片 / %d 字节），跳过" % (len(rchunks), sum(rchunks.values())), flush=True)
            continue
        print("  远端 %d 片校验通过，开始拉取 ..." % n_parts, flush=True)

        jobs = []
        for i in range(n_parts):
            rp = "%s/%s/p_%02d" % (SPLIT_DIR, sub, i)
            want = rchunks[SPLIT_DIR + "/" + sub + "/p_%02d" % i]
            jobs.append((i, rp, chunk_local(name, i), want))

        ok = {i: False for i in range(n_parts)}
        for round_i in range(2):
            todo = [j for j in jobs if not ok[j[0]]]
            if not todo:
                break
            if round_i > 0:
                print("  重试轮次 %d（%d 片）" % (round_i, len(todo)), flush=True)
            threads = []
            for idx, rp, lpath, want in todo:
                th = threading.Thread(target=fetch_with_retry, args=(pw, rp, lpath, want, 2))
                th.start()
                threads.append(th)
                while sum(t.is_alive() for t in threads) >= N_THREADS:
                    time.sleep(0.5)
            for th in threads:
                th.join()
            for idx, rp, lpath, want in todo:
                ok[idx] = lpath.exists() and lpath.stat().st_size == want

        bad = [i for i in range(n_parts) if not ok[i]]
        if bad:
            print("  %s 仍有坏分片: %s" % (name, bad), flush=True)
            continue

        # 本地合并 + SHA256 对账（远端 vs 本地）
        import subprocess
        lp = LOCAL_DIR / name
        cat = subprocess.run(["bash", "-c",
                              "cat '%s' > '%s'" % ("' '".join(str(chunk_local(name, i)) for i in range(n_parts)), lp)],
                             capture_output=True, text=True)
        got = lp.stat().st_size
        if got != expect:
            print("  %s 合并大小不符 %d != %d" % (name, got, expect), flush=True)
            continue
        _, rsha = run_fresh(pw, "sha256sum %s/%s" % (REMOTE_DIR, name))
        lsha = sha256_of(lp)
        rsha = rsha.split()[0] if rsha.split() else "?"
        print("  %s SHA256 %s (%s)" % (name, lsha[:16] + "...", "匹配 ✓" if lsha == rsha else "不匹配 ✗ 远端:" + rsha[:16]),
              flush=True)
    print("PULL-V4-DONE", flush=True)


if __name__ == "__main__":
    main()
