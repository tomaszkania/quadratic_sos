#!/usr/bin/env python3
"""Optional external-CAS availability probe.

This script is not part of the mandatory TOMS smoke-test path.  It records whether
SageMath and PARI/GP executables are visible on the local machine.  The papers include
self-contained generic baselines because external CAS installations vary substantially
between reviewers' systems.
"""

from __future__ import annotations

import shutil
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    out = root / "data" / "optional_cas_baselines.tsv"
    out.parent.mkdir(exist_ok=True)
    rows = [("system", "executable", "status")]
    for name in ["sage", "gp"]:
        exe = shutil.which(name)
        rows.append((name, exe or "", "available" if exe else "not-found"))
    out.write_text("\n".join("\t".join(row) for row in rows) + "\n")
    print(out.read_text(), end="")


if __name__ == "__main__":
    main()
