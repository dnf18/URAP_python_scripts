#!/usr/bin/env python3
"""Statistical comparison of Step 1 isotope production between two runs.

Each IsotopesLarge.dat.p1.inc*.dat holds integer COUNTS of produced isotopes ('RP id exc count') and the
simulated time 'TT' of that file. This script sums counts and times over all files of a run, gives the rate
N/T with a Poisson error, and compares two runs isotope by isotope (summed over volumes). A bootstrap over
files gives a second, independent error estimate.

Usage: python3 isotope_stats.py <run_A_dir> <run_B_dir> [--out DIR] [--nboot 2000] [--min-counts 20]
"""
import argparse, csv, glob, os, re, sys
from collections import defaultdict
import numpy as np

SYM = ("n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
       "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
       "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U").split()
def name(iso, exc):
    z, a = divmod(iso, 1000)
    return "%s-%d%s" % (SYM[z] if z < len(SYM) else "Z%d" % z, a, "" if exc == 0 else " (%.1f keV)" % exc)

def read_run(d):
    files = sorted(glob.glob(os.path.join(d, "IsotopesLarge.dat.p1.inc*.dat")),
                   key=lambda p: int(re.search(r"inc(\d+)\.dat$", p).group(1)))
    if not files: sys.exit("No IsotopesLarge.dat.p1.inc*.dat in %s" % d)
    T = np.zeros(len(files)); C = defaultdict(lambda: np.zeros(len(files)))
    for i, f in enumerate(files):
        with open(f) as fh:
            for ln in fh:
                p = ln.split()
                if not p or p[0].startswith("#"): continue
                if p[0] == "TT" and len(p) > 1: T[i] += float(p[1])
                elif p[0] == "RP" and len(p) >= 4:
                    C[(int(p[1]), round(float(p[2]), 1))][i] += float(p[3])
    return files, T, C

ap = argparse.ArgumentParser()
ap.add_argument("a"); ap.add_argument("b"); ap.add_argument("--out", default=".")
ap.add_argument("--nboot", type=int, default=2000); ap.add_argument("--min-counts", type=int, default=20)
o = ap.parse_args()
fa, Ta, Ca = read_run(o.a); fb, Tb, Cb = read_run(o.b)
print("A: %d files, total TT %.4g s; B: %d files, total TT %.4g s" % (len(fa), Ta.sum(), len(fb), Tb.sum()))
rng = np.random.default_rng(12345)
ia = rng.integers(0, len(fa), (o.nboot, len(fa))); ib = rng.integers(0, len(fb), (o.nboot, len(fb)))
TAb = Ta[ia].sum(1); TBb = Tb[ib].sum(1)
rows = []
for k in sorted(set(Ca) | set(Cb)):
    ca = Ca.get(k, np.zeros(len(fa))); cb = Cb.get(k, np.zeros(len(fb)))
    Na, Nb = ca.sum(), cb.sum()
    ra, rb = Na / Ta.sum(), Nb / Tb.sum()
    ea, eb = np.sqrt(max(Na, 1)) / Ta.sum(), np.sqrt(max(Nb, 1)) / Tb.sum()
    sig = (rb - ra) / np.hypot(ea, eb)
    if Na + Nb >= o.min_counts and Na > 0 and Nb > 0:
        qa = ca[ia].sum(1) / TAb; qb = cb[ib].sum(1) / TBb
        with np.errstate(divide="ignore", invalid="ignore"):
            q = qb / qa
        lo, hi = np.nanpercentile(q[np.isfinite(q)], [2.5, 97.5])
    else: lo = hi = np.nan
    rows.append(dict(isotope=name(*k), id=k[0], exc_keV=k[1], N_A=Na, N_B=Nb, rate_A=ra, rate_B=rb,
                     ratio_BoverA=(rb / ra if ra > 0 else np.inf if rb > 0 else np.nan),
                     boot_lo=lo, boot_hi=hi, poisson_sigma=sig))
os.makedirs(o.out, exist_ok=True)
path = os.path.join(o.out, "isotope_stats.csv")
with open(path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("Wrote", path, "(%d rows)" % len(rows))

def show(title, sel):
    print("\n" + title)
    print("%-22s %9s %9s %9s %16s %8s" % ("isotope", "N_A", "N_B", "B/A", "boot 95% CI", "sigma"))
    for r in sel:
        ci = "%.2f-%.2f" % (r["boot_lo"], r["boot_hi"]) if np.isfinite(r["boot_lo"]) else "n/a"
        print("%-22s %9.0f %9.0f %9.3g %16s %8.1f" % (r["isotope"], r["N_A"], r["N_B"], r["ratio_BoverA"], ci, r["poisson_sigma"]))
watch = {33071, 33072, 33073, 33074, 31066, 31067, 31068, 32069, 32071, 32073, 32075, 32077, 13026, 13028, 11022, 11024, 4007, 9018}
show("Watch list (As, Ga, Ge parents, Al, Na, 7Be, 18F)", [r for r in rows if r["id"] in watch])
big = [r for r in rows if r["N_A"] + r["N_B"] >= o.min_counts]
big.sort(key=lambda r: -abs(r["poisson_sigma"]))
show("Largest significant differences (N_A+N_B >= %d)" % o.min_counts, big[:30])
nsig = sum(abs(r["poisson_sigma"]) > 3 for r in big)
print("\n%d of %d isotope entries with enough counts differ by more than 3 sigma (Poisson only)" % (nsig, len(big)))
