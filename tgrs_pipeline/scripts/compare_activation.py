#!/usr/bin/env python3
"""Compare two Cosima activation files (CosmicProtons.activation.dat).

Usage: python3 compare_activation.py <A_file> <B_file> [--out DIR] [--top 40]

Format: 'VN <volume>' starts a block; 'RP <1000*Z+A> <excitation keV> <value>' lines follow.
Values are summed over volumes (per isotope and excitation), then B/A is reported.
Writes activation_comparison.csv and prints the biggest isotopes, a watch list of
parents of the lines that differ, and totals per element.
"""
import argparse, csv, os, sys
from collections import defaultdict

SYM = ("n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
       "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
       "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U").split()

def name(iso, exc):
    z, a = divmod(iso, 1000)
    s = SYM[z] if z < len(SYM) else "Z%d" % z
    return "%s-%d%s" % (s, a, "" if exc == 0 else " (%.1f keV)" % exc)

def read(path):
    d = defaultdict(float); vols = set(); vol = None; bad = 0
    with open(path) as f:
        for ln in f:
            p = ln.split()
            if not p or p[0].startswith("#"):
                continue
            if p[0] == "VN" and len(p) > 1:
                vol = p[1]; vols.add(vol)
            elif p[0] == "RP" and len(p) >= 4:
                try:
                    d[(int(p[1]), round(float(p[2]), 1))] += float(p[3])
                except ValueError:
                    bad += 1
    return d, vols, bad

ap = argparse.ArgumentParser()
ap.add_argument("a"); ap.add_argument("b")
ap.add_argument("--out", default="."); ap.add_argument("--top", type=int, default=40)
ap.add_argument("--a-label", default="A(old)"); ap.add_argument("--b-label", default="B(dc6)")
o = ap.parse_args()

A, va, ba = read(o.a); B, vb, bb = read(o.b)
print("Volumes: A %d, B %d; unreadable RP lines: A %d, B %d" % (len(va), len(vb), ba, bb))
print("Isotope entries: A %d, B %d; only in A %d, only in B %d" %
      (len(A), len(B), len(set(A) - set(B)), len(set(B) - set(A))))
sa, sb = sum(A.values()), sum(B.values())
print("Total of all values: A %.4g, B %.4g, B/A %.3f" % (sa, sb, sb / sa if sa else float("nan")))

keys = sorted(set(A) | set(B), key=lambda k: -(A.get(k, 0) + B.get(k, 0)))
rows = []
for k in keys:
    a, b = A.get(k, 0.0), B.get(k, 0.0)
    rows.append(dict(isotope=name(*k), id=k[0], exc_keV=k[1], A=a, B=b,
                     ratio_BoverA=(b / a if a > 0 else float("inf") if b > 0 else float("nan"))))
os.makedirs(o.out, exist_ok=True)
path = os.path.join(o.out, "activation_comparison.csv")
with open(path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("Wrote", path, "(%d rows)" % len(rows))

def show(title, sel):
    print("\n" + title)
    print("%-22s %12s %12s %9s" % ("isotope", o.a_label, o.b_label, "B/A"))
    for r in sel:
        print("%-22s %12.4g %12.4g %9.3g" % (r["isotope"], r["A"], r["B"], r["ratio_BoverA"]))

show("Largest %d entries" % o.top, rows[:o.top])

# parents of the lines that differ: 71As (175.7 keV), 73As, 67Ga, 68Ga, 69Ge, 71Ge, 75Ge (139.7), 71Ge isomer (198)
watch = {33071, 33072, 33073, 33074, 31066, 31067, 31068, 32069, 32071, 32073, 32075, 32077, 13026, 13028, 11022, 11024}
show("Watch list (As, Ga, Ge parents; Al and Na)", [r for r in rows if r["id"] in watch])

# per element
ea, eb = defaultdict(float), defaultdict(float)
for (iso, _), v in A.items(): ea[iso // 1000] += v
for (iso, _), v in B.items(): eb[iso // 1000] += v
print("\nTotal per element (B/A)")
for z in sorted(set(ea) | set(eb)):
    if ea[z] + eb[z] > 0.002 * (sa + sb):
        print("%-4s %12.4g %12.4g %9.3g" % (SYM[z] if z < len(SYM) else z, ea[z], eb[z], eb[z] / ea[z] if ea[z] else float("inf")))
