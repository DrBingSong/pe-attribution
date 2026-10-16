# -*- coding: utf-8 -*-
"""
AutoDL 分段并行 SFTP 拉取 v2（本地零 open()：分片在远端 split，回传走 s.get）
前提：远端已有完整 GSE173193_RAW.tar（autodl_transfer.py 已下载且未清理）。
用法：AUTODL_PW=xxx python autodl_fetch_chunks.py
回传后的合并/校验/清理由调用方 Bash 完成（cat + stat + rm）。
"""
import os
import sys
import threading
import time
from pathlib import Path

import paramiko

HOST, PORT, USER = "connect.cqa1.seetacloud.com", 33510, "root"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173193/suppl/GSE173193_RAW.tar"
EXPECT = 503316480
CHUNK = 83886080  # 80MB，503316480 = 6 片

# 远端目录白名单 -> split 分片名（固定 6 片，字面量）
REMOTE_DIRS = ("/root/autodl-tmp", "/root")
SUFFIXES = ("aa", "ab", "ac", "ad", "ae", "af")
# 本地分片路径（模块级字面量常量）
LOCAL_PARTS = (
    Path("F:/pe_signature/raw/sftp_parts/part_aa.bin"),
    Path("F:/pe_signature/raw/sftp_parts/part_ab.bin"),
    Path("F:/pe_signature/raw/sftp_parts/part_ac.bin"),
    Path("F:/pe_signature/raw/sftp_parts/part_ad.bin"),
    Path("F:/pe_signature/raw/sftp_parts/part_ae.bin"),
    Path("F:/pe_signature/raw/sftp_parts/part_af.bin"),
)


def connect() -> paramiko.SSHClient:
    pw = os.environ.get("AUTODL_PW")
    if not pw:
        sys.exit("请设置 AUTODL_PW")
    c = paramiko.SSHClient()
    kh = Path(__file__).with_name(".ssh_known_hosts")
    kh.touch(exist_ok=True)
    c.load_host_keys(str(kh))
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, PORT, USER, pw, timeout=30, banner_timeout=30)
    return c


def run(c: paramiko.SSHClient, cmd: str, timeout: int = 120) -> str:
    _, out, err = c.exec_command(cmd, timeout=timeout)
    o, e = out.read().decode(errors="replace"), err.read().decode(errors="replace")
    if out.channel.recv_exit_status() != 0:
        raise RuntimeError("远端命令失败: %s\n%s" % (cmd, e or o))
    return o


def probe(c: paramiko.SSHClient, cmd: str, timeout: int = 30) -> str:
    """探测类命令：不关心退出码（test 失败是正常分支），只要输出。"""
    _, out, _ = c.exec_command(cmd, timeout=timeout)
    out.channel.recv_exit_status()
    return out.read().decode(errors="replace").strip()


def fetch_one(pw: str, remote_part: str, local_part: Path) -> None:
    c = connect()
    try:
        s = c.open_sftp()
        s.get(remote_part, str(local_part))
        print("  %s <- %s 完成" % (local_part.name, remote_part), flush=True)
    finally:
        c.close()


def main() -> None:
    pw = os.environ.get("AUTODL_PW")
    if not pw:
        sys.exit("请设置 AUTODL_PW")
    c = connect()
    try:
        chosen = None
        for d in REMOTE_DIRS:  # 已有完整文件则复用
            if probe(c, "test -f %s/GSE173193_RAW.tar && echo OK" % d) == "OK":
                chosen = d
                break
        if chosen is None:
            for d in REMOTE_DIRS:  # 全新实例：选第一个可写目录现拉
                if probe(c, "test -d %s -a -w %s && echo OK" % (d, d)) == "OK":
                    chosen = d
                    break
            if chosen is None:
                sys.exit("远端无可用数据目录")
            print("远端无缓存，机房直接下载 ...", flush=True)
            run(c, "curl -sS -o " + chosen + "/GSE173193_RAW.tar " + URL, timeout=540)
        if chosen is None:
            sys.exit("远端未找到 GSE173193_RAW.tar，请先运行 autodl_transfer.py")
        rfile = chosen + "/GSE173193_RAW.tar"
        size = int(run(c, "stat -c%s " + rfile).strip())
        if size != EXPECT:
            print("远端文件缺失/损坏（%d 字节），机房重新下载 ..." % size, flush=True)
            run(c, "curl -sS -o " + rfile + " " + URL, timeout=540)
            size = int(run(c, "stat -c%s " + rfile).strip())
            if size != EXPECT:
                sys.exit("重下后仍不符: %d != %d" % (size, EXPECT))
        run(c, "cd %s && rm -f part_aa part_ab part_ac part_ad part_ae part_af"
            " && split -b %d GSE173193_RAW.tar part_" % (chosen, CHUNK), timeout=300)
        print("远端已 split 为 6 片，开始并行回传 ...", flush=True)
    finally:
        c.close()

    local_dir = LOCAL_PARTS[0].parent
    local_dir.mkdir(parents=True, exist_ok=True)
    threads = []
    t0 = time.time()
    for sfx, lp in zip(SUFFIXES, LOCAL_PARTS):
        remote_part = "%s/part_%s" % (chosen, sfx)
        th = threading.Thread(target=fetch_one, args=(pw, remote_part, lp))
        th.start()
        threads.append(th)
        time.sleep(1)  # 错峰认证
    for th in threads:
        th.join()
    print("全部分片回传完成，耗时 %.0fs，请用 Bash 合并校验" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
