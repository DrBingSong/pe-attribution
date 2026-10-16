# -*- coding: utf-8 -*-
"""从 AutoDL 拉取 Phase 0 结果：小文件直传 + 大文件远端 split 后 6 线程回传。"""
import os
import threading
import time
from pathlib import Path

import paramiko

HOST, PORT, USER = "connect.cqa1.seetacloud.com", 33510, "root"
BIG_REMOTE = "/root/pe_signature/results/phase0/atts_scimilarity_full.npz"
BIG_EXPECT = 556255994
CHUNK = 70000000  # 70MB × 8
SPLIT_DIR = "/root/pe_signature/tmp_split"

SMALL_FILES = {
    "/root/pe_signature/results/phase0/atts_scfoundation.npz": Path("F:/pe_signature/results/phase0/atts_scfoundation.npz"),
    "/root/pe_signature/results/phase0/atts_scimilarity.npz": Path("F:/pe_signature/results/phase0/atts_scimilarity.npz"),
    "/root/pe_signature/results/phase0/atts_scimilarity_full_meta.npz": Path("F:/pe_signature/results/phase0/atts_scimilarity_full_meta.npz"),
}
CHUNK_LOCAL = tuple(Path("F:/pe_signature/dl_tmp/chunk_%02d.bin" % i) for i in range(8))


def connect(pw):
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    kh = Path(__file__).with_name(".ssh_known_hosts")
    kh.touch(exist_ok=True)
    c.load_host_keys(str(kh))
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, PORT, USER, pw, timeout=30)
    return c


def run(c, cmd, timeout=120):
    _, o, e = c.exec_command(cmd, timeout=timeout)
    o.channel.recv_exit_status()
    return o.read().decode(errors="replace")


def main():
    pw = os.environ["AUTODL_PW"]
    for p in CHUNK_LOCAL:
        p.parent.mkdir(parents=True, exist_ok=True)
    c = connect(pw)

    for rp, lp in SMALL_FILES.items():
        lp.parent.mkdir(parents=True, exist_ok=True)
        s = c.open_sftp()
        s.get(rp, str(lp))
        s.close()
        print("小文件完成:", lp.name, lp.stat().st_size, "字节", flush=True)

    run(c, "rm -rf %s && mkdir -p %s && split -b %d %s %s/chunk_"
        % (SPLIT_DIR, SPLIT_DIR, CHUNK, BIG_REMOTE, SPLIT_DIR), timeout=180)
    names = sorted(run(c, "ls %s" % SPLIT_DIR).split())
    print("远端分片:", names, flush=True)
    c.close()

    threads = []
    for i, name in enumerate(names):
        th = threading.Thread(target=fetch_one, args=(pw, SPLIT_DIR + "/" + name, CHUNK_LOCAL[i]))
        th.start()
        threads.append(th)
        time.sleep(1)
    for th in threads:
        th.join()
    print("全部分片回传完成", flush=True)


def fetch_one(pw, remote, local):
    c = connect(pw)
    try:
        s = c.open_sftp()
        s.get(remote, str(local))
        print("  ", local.name, local.stat().st_size, "字节", flush=True)
    finally:
        c.close()


if __name__ == "__main__":
    main()
