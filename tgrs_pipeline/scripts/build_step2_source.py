#!/usr/bin/env python3
"""
Step 1 (run via `mcosima -t N CosmicProtons.Step1.source`) writes one isotope
production file per thread: IsotopesLarge.dat.p1.inc1.dat ... incN.dat. The
repo's checked-in CosmicProtons.Step2.source hardcodes exactly 5 of these
(A.IsotopeProductionFile lines) because that's how many threads the original
example happened to be run with.

This script regenerates Step 2's source file by globbing for whatever
IsotopesLarge.dat.p1.inc*.dat files Step 1 actually produced, so it's
correct regardless of how many threads/instances you run Step 1 with.

Usage:
    python3 build_step2_source.py <template=CosmicProtons.Step2.source> <run_dir>
Writes:
    <run_dir>/CosmicProtons.Step2.generated.source
"""

import re
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 3:
        print("Usage: build_step2_source.py <template.source> <run_dir>")
        return 1

    template_path = Path(sys.argv[1])
    run_dir = Path(sys.argv[2])

    isotope_files = sorted(
        run_dir.glob("IsotopesLarge.dat.p1.inc*.dat"),
        key=lambda p: int(re.search(r"\.inc(\d+)\.dat$", p.name).group(1)),
    )

    if not isotope_files:
        print(f"ERROR: no IsotopesLarge.dat.p1.inc*.dat files found in {run_dir}. "
              f"Did Step 1 (mcosima on CosmicProtons.Step1.source) run and finish first?")
        return 1

    print(f"Found {len(isotope_files)} isotope production files:")
    for f in isotope_files:
        print(f"  {f.name}")

    text = template_path.read_text()

    # Drop every existing "A.IsotopeProductionFile ..." line from the
    # template (whatever fixed count it happened to hardcode) and replace
    # with the actual set we found.
    lines = [ln for ln in text.splitlines()
             if not ln.strip().startswith("A.IsotopeProductionFile")]

    insert_at = next(
        (i for i, ln in enumerate(lines) if ln.strip().startswith("A.ActivationFile")),
        len(lines),
    )

    new_lines = [f"A.IsotopeProductionFile   {f.name}" for f in isotope_files]
    lines[insert_at:insert_at] = new_lines

    out_path = run_dir / "CosmicProtons.Step2.generated.source"
    out_path.write_text("\n".join(lines) + "\n")
    print(f"\nWrote {out_path} with {len(isotope_files)} IsotopeProductionFile line(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
