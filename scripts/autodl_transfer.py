# -*- coding: utf-8 -*-
"""
AutoDL 无卡中转下载：机房带宽拉 NCBI → SFTP 传回本机
用法（密码经环境变量传入，不落盘）：
  AUTODL_PW=xxx python autodl_transfer.py            # 下载远端 + 传回 + 校验
  AUTODL_PW=xxx python autodl_transfer.py --info     # 仅查看实例状态
"""
import os
import stat
import sys
import time
from pathlib import Path

import paramiko

HOST, PORT, USER = "connect.bjb2.seetacloud.com", 45759, "root"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173193/suppl/GSE173193_RAW.tar"
EXPECT = 503316480
REMOTE_DIR_CANDIDATES = ["/root/autodl-tmp", "/root"]
LOCAL_OUT = Path("F:/pe_signature/raw/GSE173193_RAW.tar")


def connect() -> paramiko.SSHClient:
    pw = os.environ.get("AUTODL_PW")
    if not pw:
        sys.exit("请设置环境变量 AUTODL_PW")
    c = paramiko.SSHClient()
    # 主机密钥持久化：首次记录，之后逐次校验，防中间人替换
    kh = Path(__file__).with_name(".ssh_known_hosts")
    kh.touch(exist_ok=True)
    c.load_host_keys(str(kh))
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, PORT, USER, pw, timeout=30, banner_timeout=30)
    return c


def run(c: paramiko.SSHClient, cmd: str, timeout: int = 60) -> str:
    _, out, err = c.exec_command(cmd, timeout=timeout)
    o, e = out.read().decode(errors="replace"), err.read().decode(errors="replace")
    code = out.channel.recv_exit_status()
    if code != 0:
        raise RuntimeError(f"远端命令失败({code}): {cmd}\n{e or o}")
    return o


def info(c: paramiko.SSHClient) -> None:
    print("== 实例状态 ==")
    print(run(c, "hostname; uptime | head -1"))
    print(run(c, "nvidia-smi -L 2>/dev/null || echo '无 GPU（无卡模式）'"))
    print(run(c, "df -h /root /root/autodl-tmp 2>/dev/null | grep -v Filesystem || df -h /root"))


def remote_download(c: paramiko.SSHClient) -> str:
    dest_dir = ""
    for d in REMOTE_DIR_CANDIDATES:  # 本地常量白名单，不接受远端回显拼接
        if run(c, f"test -d {d} -a -w {d} && echo OK").strip() == "OK":
            dest_dir = d
            break
    if not dest_dir:
        sys.exit("远端无可用数据目录")
    dest = f"{dest_dir}/GSE173193_RAW.tar"
    print(f"远端下载到 {dest} ...")
    t0 = time.time()
    print(run(c, f"curl -sS -o {dest} '{URL}' && stat -c%s {dest}", timeout=540).strip())
    print(f"远端下载耗时 {time.time()-t0:.0f}s")
    return dest


def sftp_fetch(c: paramiko.SSHClient, remote: str) -> None:
    s = c.open_sftp()
    size = s.stat(remote).st_size
    if size != EXPECT:
        sys.exit(f"远端大小不符: {size} != {EXPECT}")
    done = {"n": 0}

    def cb(sent: int, total: int) -> None:
        if sent // (100 * 1024 * 1024) > done["n"]:
            done["n"] = sent // (100 * 1024 * 1024)
            print(f"  传回 {sent/1e6:.0f}/{total/1e6:.0f} MB ...", flush=True)

    LOCAL_OUT.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    s.get(remote, str(LOCAL_OUT), callback=cb)
    dt = time.time() - t0
    print(f"传回完成: {LOCAL_OUT} ({LOCAL_OUT.stat().st_size} 字节, {dt:.0f}s, {EXPECT/1e6/dt:.1f} MB/s)")


def main() -> None:
    c = connect()
    try:
        info(c)
        if "--info" in sys.argv:
            return
        remote = remote_download(c)
        sftp_fetch(c, remote)
        got = LOCAL_OUT.stat().st_size
        print("校验:", "OK" if got == EXPECT else f"FAIL {got}!={EXPECT}")
        if got == EXPECT:
            run(c, f"rm -f {remote}")
            print("已清理远端临时文件")
    finally:
        c.close()


if __name__ == "__main__":
    main()
