#!/usr/bin/env python3
"""Fill a base pipeline deck and study decks from a finished Gc(l0) calibration.

Usage: pfd_fill_base.py CAL_ROOT L0 BASE.i [STUDY.i ...]

Reads ``<CAL_ROOT>/calibrated.json``, takes the entry with ``l0 == L0`` and
writes its ``Gc`` into ``Gc = `` of BASE.i's [Plate] block (and ``l0``);
in every STUDY.i the placeholders ``@GC_HALF@`` and ``@GC_DOUBLE@`` become
half and twice that Gc. Every edit is printed, so the record shows where
the numbers came from.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path


def main(argv):
    if len(argv) < 4:
        print(__doc__)
        return 2
    cal = json.loads((Path(argv[1]) / "calibrated.json").read_text())
    l0 = float(argv[2])
    entry = next(c for c in cal if float(c["l0"]) == l0 and c.get("Gc") is not None)
    gc = float(entry["Gc"])
    base = Path(argv[3])
    text = base.read_text()
    text, n1 = re.subn(r"^(\s*Gc\s*=\s*)[^#\n]*", rf"\g<1>{gc:.6g}        ", text, count=1, flags=re.M)
    text, n2 = re.subn(r"^(\s*l0\s*=\s*)[^#\n]*", rf"\g<1>{l0:g}            ", text, count=1, flags=re.M)
    base.write_text(text)
    print(f"{base}: Gc = {gc:.6g}, l0 = {l0:g} (from {argv[1]}, fracture strain {entry.get('fracture_strain')})")
    for s in argv[4:]:
        p = Path(s)
        t = p.read_text().replace("@GC_HALF@", f"{0.5 * gc:.6g}").replace("@GC_DOUBLE@", f"{2.0 * gc:.6g}")
        p.write_text(t)
        print(f"{p}: filled Gc x0.5 = {0.5 * gc:.6g}, x2 = {2.0 * gc:.6g}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
