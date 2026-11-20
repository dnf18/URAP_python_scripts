#!/usr/bin/env python3
"""Merge the per-task outputs of one stage into the run folder.

    python3 merge_task_outputs.py <run_dir> <photons|step1|step3>

For every <run_dir>/<stage>_NNNN folder it finds the instance files
(<base>.p1.inc*.id*.sim.gz, plus IsotopesLarge.dat.p1.inc*.dat for step1), links them
into <run_dir> under one global numbering, and writes a new stub <base>.p1.sim.gz that
lists them, in the same format as the stub mcosima wrote (header + one 'IN' line per file).
It stops with an error if two instance files share a SimulationSeed (identical events).
"""
import glob, gzip, hashlib, os, re, sys, zlib
from collections import Counter

STAGES = {
    "photons": ("CosmicPhotons", False),
    "step1": ("CosmicProtonsPrompt", True),
    "step3": ("CosmicProtons.Decays", False),
}


def inc_index(path):
    m = re.search(r"\.inc(\d+)\.", os.path.basename(path))
    return int(m.group(1)) if m else 0


def fingerprint(path):
    """Return (seed, event_hash, ok). ok=False if the file cannot be read at all.
    event_hash covers the first events of the file, so two instances that started from the
    same random seed (identical events) get the same hash even if no seed line exists."""
    seed, h, n, lines = None, hashlib.sha1(), 0, 0
    in_events = False
    try:
        with gzip.open(path, "rt", errors="replace") as f:
            for line in f:
                lines += 1
                if not in_events:
                    if line.startswith("SimulationSeed") and seed is None:
                        seed = line.split()[1]
                    if line.startswith("SE"):
                        in_events = True
                    elif lines > 5000:
                        break
                    else:
                        continue
                h.update(line.encode("utf-8", "replace"))
                n += 1
                if n >= 2000:
                    break
    except (OSError, EOFError, zlib.error):
        if lines == 0:
            return None, None, False          # unreadable from the first byte
    return seed, (h.hexdigest() if n else None), True


def isotope_ok(path):
    """False for isotope files cut off mid-write (job killed at the time limit): they do not end
    with a newline or have an 'RP' line with fewer than 4 columns, which crashes cosima Step 2."""
    try:
        with open(path, "rb") as f:
            data = f.read()
    except OSError:
        return False
    if not data.endswith(b"\n"):
        return False
    for line in data.split(b"\n"):
        if line.startswith(b"RP") and len(line.split()) < 4:
            return False
    return True


def last_ti(path):
    """Last 'TI' (event time) in an instance file, tolerant of truncation; None if there is none."""
    t = None
    try:
        with gzip.open(path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("TI "):
                    try:
                        t = float(line.split()[1])
                    except (ValueError, IndexError):
                        pass
    except (OSError, EOFError, zlib.error):
        pass
    return t


def isotope_tt(path):
    """Observation time (header line 'TT <seconds>') of an isotope file, None if missing."""
    try:
        with open(path, "r", errors="replace") as f:
            for i, line in enumerate(f):
                if line.startswith("TT"):
                    return float(line.split()[1])
                if i > 50:
                    break
    except (OSError, ValueError, IndexError):
        pass
    return None


def link(src, dst):
    if os.path.lexists(dst):
        os.remove(dst)
    os.symlink(os.path.abspath(src), dst)


def main():
    if len(sys.argv) != 3 or sys.argv[2] not in STAGES:
        print(__doc__)
        return 1
    run_dir, stage = os.path.abspath(sys.argv[1]), sys.argv[2]
    base, iso = STAGES[stage]
    tasks = sorted(glob.glob(os.path.join(run_dir, f"{stage}_[0-9][0-9][0-9][0-9]")))
    if not tasks:
        print(f"ERROR: no {stage}_NNNN folders in {run_dir}")
        return 1

    sims, isos, header, trailer, template = [], [], None, [], None
    for t in tasks:
        sims += sorted((p for p in glob.glob(os.path.join(t, f"{base}.p1.inc*.id*.sim.gz"))
                        if os.path.getsize(p) > 0), key=inc_index)
        if iso:
            isos += sorted((p for p in glob.glob(os.path.join(t, "IsotopesLarge.dat.p1.inc*.dat"))
                            if os.path.getsize(p) > 0), key=inc_index)
        stub = os.path.join(t, f"{base}.p1.sim.gz")
        if header is None and os.path.isfile(stub):
            try:
                lines = gzip.open(stub, "rt", errors="replace").read().splitlines()
            except (OSError, EOFError):
                continue
            idx = [i for i, l in enumerate(lines) if l.split()[:1] == ["IN"]]
            if idx:
                header = lines[:idx[0]]
                template = lines[idx[0]].split()[0]
                trailer = [l for l in lines[idx[-1] + 1:] if l.strip()]
    print(f"{stage}: {len(tasks)} task folders, {len(sims)} instance files"
          + (f", {len(isos)} isotope files" if iso else ""))
    if not sims or header is None:
        print("ERROR: no usable output or no stub header found (did every task fail?)")
        return 1
    if iso and not isos:
        print("ERROR: no isotope production files found")
        return 1
    if iso:
        bad_iso = [p for p in isos if not isotope_ok(p)]
        if bad_iso:
            print(f"WARNING: skipping {len(bad_iso)} truncated isotope file(s), e.g. {bad_iso[:3]}")
        isos = [p for p in isos if p not in set(bad_iso)]
        # a file with a missing or near-zero TT would blow up Step 2 (counts / ~0 s)
        tts = {p: isotope_tt(p) for p in isos}
        vals = sorted(t for t in tts.values() if t)
        if vals:
            med = vals[len(vals) // 2]
            low = [p for p, t in tts.items() if not t or t < 0.1 * med]
            if low:
                print(f"WARNING: skipping {len(low)} isotope file(s) with missing/near-zero TT "
                      f"(median TT {med:g} s), e.g. {[(os.path.basename(os.path.dirname(p)), tts[p]) for p in low[:5]]}")
                lowset = set(low)
                isos = [p for p in isos if p not in lowset]
            print(f"isotope TT: median {med:g} s, total {sum(t for p, t in tts.items() if p in set(isos) and t):.0f} s")
        print(f"usable isotope files: {len(isos)}")
        if not isos:
            print("ERROR: no usable isotope files")
            return 1
    # remove links from an earlier merge so no stale entries remain (only symlinks are touched)
    for pat in (f"{base}.p1.inc*.id1.sim.gz", "IsotopesLarge.dat.p1.inc*.dat"):
        if pat.startswith("Isotopes") and not iso:
            continue
        for old in glob.glob(os.path.join(run_dir, pat)):
            if os.path.islink(old):
                os.remove(old)

    fps = [fingerprint(p) for p in sims]
    bad = [p for p, fp in zip(sims, fps) if not fp[2]]
    if bad:
        print(f"WARNING: skipping {len(bad)} unreadable instance file(s), e.g. {bad[:3]}")
    good = [(p, fp) for p, fp in zip(sims, fps) if fp[2]]
    sims = [p for p, _ in good]
    seeds = [fp[0] for _, fp in good if fp[0]]
    hashes = [fp[1] for _, fp in good if fp[1]]
    print(f"usable instance files: {len(sims)}; with seed line: {len(seeds)} "
          f"({len(set(seeds))} distinct); with events: {len(hashes)} ({len(set(hashes))} distinct)")
    if not seeds and not hashes:
        print("ERROR: found neither seeds nor events in any file, cannot check for duplicates. "
              "Send the output of: zcat <one inc file> | head -40")
        return 4
    ndup = (len(seeds) - len(set(seeds))) + (len(hashes) - len(set(hashes)))
    if ndup:
        print(f"ERROR: {ndup} duplicate seed/event fingerprints, so some instances are identical copies. "
              f"Not merging.")
        return 3
    if not sims:
        print("ERROR: no usable instance files")
        return 1

    if stage == "step3":   # decays: report the observation time actually present
        tis = [last_ti(p) for p in sims]
        done = [t for t in tis if t]
        total = sum(done)
        print(f"step3 observation time: {total:.1f} s from {len(done)} of {len(sims)} instance files "
              f"(target {os.environ.get('TARGET_TIME', 'unknown')} s)")
        if not done:
            print("ERROR: none of the decay files contains a complete observation (all TI are 0).")
            return 5
        sims = [p for p, t in zip(sims, tis) if t]   # keep only files with events

    names = []
    for k, p in enumerate(sims, 1):
        name = f"{base}.p1.inc{k}.id1.sim.gz"
        link(p, os.path.join(run_dir, name))
        names.append(name)
    for k, p in enumerate(isos, 1):
        link(p, os.path.join(run_dir, f"IsotopesLarge.dat.p1.inc{k}.dat"))

    out = os.path.join(run_dir, f"{base}.p1.sim.gz")
    if os.path.lexists(out):
        os.remove(out)
    with gzip.open(out, "wt") as f:
        f.write("\n".join(header) + "\n")
        f.write("\n".join(f"{template} {n}" for n in names) + "\n")
        if trailer:
            f.write("\n".join(trailer) + "\n")
    print(f"Wrote {out} listing {len(names)} files (first line: {template} {names[0]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
