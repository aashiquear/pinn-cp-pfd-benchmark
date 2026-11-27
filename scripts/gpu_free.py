#!/usr/bin/env python3
"""List the GPUs of this machine with their free memory and name those that can take a part.

usage: scripts/gpu_free.py [<GB needed per card>] [<most cards to use>]

Reads `nvidia-smi` (or `rocm-smi` on AMD). Prints one line per card and, for the cards
with at least the given free memory (default 0: all), the override to pass to a run:

    global.rank_devices="cuda:1 cuda:4 cuda:6"

The cards with the most free memory come first; at most the given number (default 4) are
named. Memory needed per card for the wide column decks: see the header of each deck
(about 36 kB per element, divided by the number of cards). Leave spare room: the solver
reports the memory it allocated, the driver reserves somewhat more.
"""

from __future__ import annotations

import shutil
import subprocess
import sys


def cards() -> list[tuple[int, str, float, float]]:
    """(index, name, free GB, total GB) of every card."""
    if shutil.which("nvidia-smi"):
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=index,name,memory.free,memory.total", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, check=True,
        ).stdout
        rows = []
        for line in out.strip().splitlines():
            index, name, free, total = [part.strip() for part in line.split(",")]
            rows.append((int(index), name, float(free) / 1024.0, float(total) / 1024.0))
        return rows
    if shutil.which("rocm-smi"):
        import json

        out = subprocess.run(["rocm-smi", "--showmeminfo", "vram", "--json"], capture_output=True, text=True, check=True).stdout
        rows = []
        for key, entry in sorted(json.loads(out).items()):
            if not key.startswith("card"):
                continue
            total = float(entry.get("VRAM Total Memory (B)", 0)) / 2**30
            used = float(entry.get("VRAM Total Used Memory (B)", 0)) / 2**30
            rows.append((int(key[4:]), "AMD", total - used, total))
        return rows
    raise SystemExit("neither nvidia-smi nor rocm-smi was found.")


def main() -> None:
    need = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
    most = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    rows = cards()
    for index, name, free, total in rows:
        mark = "free enough" if free >= need else "busy"
        print(f"cuda:{index}  {name:28s} {free:7.1f} GB free of {total:6.1f} GB   {mark}")
    chosen = sorted((r for r in rows if r[2] >= need), key=lambda r: -r[2])[:most]
    if not chosen:
        raise SystemExit(f"no card has {need:g} GB free.")
    devices = " ".join(f"cuda:{r[0]}" for r in sorted(chosen))
    print(f'\n{len(chosen)} card(s) with at least {need:g} GB free:\n  global.rank_devices="{devices}"')


if __name__ == "__main__":
    main()
