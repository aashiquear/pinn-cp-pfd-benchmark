#!/usr/bin/env python3
"""Run a command in its own process group, log its memory, stop it before memory runs out.

usage: scripts/memory_guard.py <min available GB> <max seconds> <command ...>

Every 5 s one line: seconds, memory still available, summed resident memory of the
command and its children. The whole group is killed when the available memory falls
below the floor or the time limit passes, and at the end if a child was left behind.
"Available" is MemAvailable of the machine, or, inside a container with a memory limit
(cgroup v2), the limit minus the current use. The output of the command goes to
/tmp/guard_out.txt. Written after a run on the GPU host that left it unresponsive.
"""
import os, signal, subprocess, sys, time

def available_gb():
    try:
        limit = open("/sys/fs/cgroup/memory.max").read().strip()
        if limit != "max":
            return (int(limit) - int(open("/sys/fs/cgroup/memory.current").read())) / 2**30
    except OSError:
        pass
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable"):
            return int(line.split()[1]) / 2**20
    return 0.0

def group_rss_gb(pgid):
    total = 0
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            if os.getpgid(int(pid)) != pgid:
                continue
            for line in open(f"/proc/{pid}/status"):
                if line.startswith("VmRSS"):
                    total += int(line.split()[1])
        except (OSError, ValueError):
            pass
    return total / 2**20

floor, limit = float(sys.argv[1]), float(sys.argv[2])
start_free = available_gb()
proc = subprocess.Popen(sys.argv[3:], start_new_session=True, stdout=open("/tmp/guard_out.txt", "w"), stderr=subprocess.STDOUT)
tic, low, peak, killed, last = time.time(), start_free, 0.0, "", 0.0
while proc.poll() is None:
    free, rss = available_gb(), group_rss_gb(proc.pid)
    low, peak = min(low, free), max(peak, rss)
    now = time.time() - tic
    if now - last >= 5.0:
        print("t %5.1f s  available %5.1f GB  group RSS %5.2f GB" % (now, free, rss), flush=True)
        last = now
    if free < floor or now > limit:
        killed = "low memory" if free < floor else "time limit"
        os.killpg(proc.pid, signal.SIGKILL)
        break
    time.sleep(0.25)
proc.wait()
time.sleep(2)
left = group_rss_gb(proc.pid)
print("exit", proc.returncode, killed, "| wall %.0f s | available at start %.1f GB, lowest %.1f GB, after %.1f GB | peak group RSS %.2f GB | group RSS left behind %.2f GB" % (time.time() - tic, start_free, low, available_gb(), peak, left))
if left > 0.01:
    os.killpg(proc.pid, signal.SIGKILL)
    print("left-over processes of the group killed")
