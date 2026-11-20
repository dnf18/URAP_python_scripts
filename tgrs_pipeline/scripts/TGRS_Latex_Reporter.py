#!/usr/bin/env python3
"""
TGRS_Latex_Reporter.py

Validation report in the LaTeX layout (table of contents, numbered sections, page numbers, links,
paper-style Table 5/6 appendix), built entirely from OUR pipeline's data and OUR analysis code:

    * histograms from the TGRS program's output ROOT files (TGRS, Sim, Component_N)
    * TGRS_Line_Compare.py -> net_rate (line = window +-1.5 FWHM minus sideband continuum),
                              group_summary (group-pattern analysis)
    * the A-versus-B tolerance test (95% of bins within 10% or 3 sigma, bands within 5% or 3 sigma)

Every section starts with a short "what this shows / how to read it" box.

Usage (same inputs as TGRS_ROOT_Report.py):
    python3 TGRS_Latex_Reporter.py A/TGRS_output.root B/TGRS_output.root --tgrs inputs_dc6/TGRS.dat \
        --ref-label "MEGAlib 10.2 (Geant4 10.2.3)" --test-label "MEGAlib dc6 (Geant4 11.4.3)" \
        --lines tgrs_lines.csv --out analysis/latex_report

Result: <out>/TGRS_Validation_Report.pdf (and the .tex and figures, so it can also be compiled elsewhere,
e.g. in Overleaf). Needs numpy, scipy, matplotlib, uproot and a LaTeX installation (pdflatex with booktabs,
longtable, pdflscape; on Savio: module load texlive, if available).
It can also be called from TGRS_Reporter:  TGRS_Reporter(...).generate_latex_pdf(root_a, root_b, tgrs, ...)
"""

import argparse
import csv
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as sst

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def root_hist(path, name):
    """(contents, edges, errors) of one histogram of a ROOT file."""
    import uproot
    with uproot.open(path) as f:
        keys = [k.split(";")[0] for k in f.keys()]
        if name not in keys:
            raise KeyError(f"'{name}' not in {path}. Available: {keys}")
        h = f[name]
        return (np.asarray(h.values(), dtype=float), np.asarray(h.axis().edges(), dtype=float),
                np.asarray(h.errors(), dtype=float))


def remap_density(values, errors, src_edges, dst_edges):
    """Map a per-keV density histogram onto other bins (width-averaged); bins not covered become NaN."""
    n = len(dst_edges) - 1
    out = np.full(n, np.nan); oerr = np.full(n, np.nan)
    for i in range(n):
        lo, hi = dst_edges[i], dst_edges[i + 1]
        if lo < src_edges[0] - 1e-9 or hi > src_edges[-1] + 1e-9:
            continue
        ov = np.clip(np.minimum(src_edges[1:], hi) - np.maximum(src_edges[:-1], lo), 0, None)
        w = ov / (hi - lo)
        out[i] = np.sum(values * w); oerr[i] = np.sqrt(np.sum((errors * w) ** 2))
    return out, oerr


import TGRS_Line_Compare as LC                                               # noqa: E402

plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 6.5,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})

SPLIT = 245.0                       # data bins below, simulation bins above
BANDS = [(0.0, 100.0), (100.0, 1000.0), (1000.0, math.inf)]
COMP_BANDS = [(50, 100), (100, 300), (300, 1000), (1000, 3000), (3000, 8195)]
COL = {"A": "tab:blue", "B": "tab:red", "D": "black"}


# ------------------------------------------------------------------ text helpers

def esc(s):
    s = str(s).replace("–", "--").replace("−", "-")
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("$", r"\$"), ("#", r"\#"),
                 ("_", r"\_"), ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"), ("^", r"\^{}")):
        s = s.replace(a, b)
    for a, b in (("β", r"$\beta$"), ("τ", r"$\tau$"), ("∼", r"$\sim$"), ("±", r"$\pm$"), ("γ", r"$\gamma$"),
                 ("α", r"$\alpha$"), ("µ", r"$\mu$"), ("μ", r"$\mu$"), ("→", r"$\to$"), ("≥", r"$\ge$"), ("≤", r"$\le$"),
                 ("×", r"$\times$"), ("°", r"$^\circ$")):
        s = s.replace(a, b)
    return s


def nuc(s):
    """73mGe(67-13) -> $^{73m}$Ge(67--13)"""
    s = esc(s)
    return re.sub(r"(?<![\d.])(\d+m?)([A-Z][a-z]?)\(", lambda m: f"$^{{{m.group(1)}}}$" + m.group(2) + "(", s)


def sci(x, n=4):
    if x is None or not np.isfinite(x):
        return "--"
    if x == 0:
        return "0"
    if 1e-2 <= abs(x) < 1e5:
        return f"{x:.{n}g}" if abs(x) < 1e3 else f"{x:,.1f}".replace(",", r"\,")
    m, e = f"{x:.{n - 1}e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def pval(p):
    return "$<10^{-300}$" if p < 1e-300 else (sci(p, 3))


def pm(v, e):
    if e is None or not np.isfinite(e) or e <= 0:
        return sci(v, 3)
    d = max(0, 1 - int(math.floor(math.log10(e))))
    return f"{v:.{d}f}$\\pm${e:.{d}f}"


def paren(v, e):
    """paper style 0.74(6)"""
    if e is None or not np.isfinite(e) or e <= 0:
        return sci(v, 3)
    lead1 = str(f"{e:.6e}")[0] == "1"
    d = max(0, -int(math.floor(math.log10(e))) + (1 if lead1 else 0))
    return f"{v:.{d}f}({int(round(e * 10 ** d))})"


def sigtxt(x):
    return f"{x:+.1f}".replace("-", "$-$")


def guide(text):
    return "\\guide{" + text + "}\n"


def table(rows, head=None, spec=None, caption=None):
    cols = len(rows[0]) if rows else len(head)
    spec = spec or ("l" + "r" * (cols - 1))
    out = ["\\begin{center}\\small", f"\\begin{{tabular}}{{{spec}}}", "\\toprule"]
    if head:
        out += [" & ".join(head) + " \\\\", "\\midrule"]
    for r in rows:
        out.append(r if isinstance(r, str) else " & ".join(str(c) for c in r) + " \\\\")
    out += ["\\bottomrule", "\\end{tabular}", "\\end{center}"]
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------ statistics helpers

def integ(v, e, edges, lo, hi):
    ov = np.clip(np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo), 0, None)
    ok = np.isfinite(v)
    return float(np.sum(v[ok] * ov[ok])), float(np.sqrt(np.sum((e[ok] * ov[ok]) ** 2)))


def gof(x, ex, y, ey, edges):
    """Pull = (x - y)/sigma per bin, chi2, per band numbers, KS."""
    ok = np.isfinite(x) & np.isfinite(y)
    sig = np.sqrt(np.where(ok, ex ** 2 + ey ** 2, 0.0))
    used = ok & (sig > 0) & ~((x == 0) & (y == 0))
    pulls = np.zeros_like(x)
    pulls[used] = (x[used] - y[used]) / sig[used]
    chi2 = float(np.sum(pulls[used] ** 2)); ndf = int(used.sum())
    w = np.diff(edges); mid = 0.5 * (edges[:-1] + edges[1:])
    k = int(np.argmax(np.abs(pulls)))
    cx = np.cumsum(np.where(ok, x, 0) * w); cy = np.cumsum(np.where(ok, y, 0) * w)
    r = dict(pulls=pulls, used=used, chi2=chi2, ndf=ndf, p=float(sst.chi2.sf(chi2, ndf)) if ndf else float("nan"),
             n3=int(np.sum(np.abs(pulls) > 3)), worst=(edges[k], edges[k + 1]), worst_pull=float(pulls[k]),
             ks=float(np.max(np.abs(cx / cx[-1] - cy / cy[-1]))) if cx[-1] > 0 and cy[-1] > 0 else float("nan"),
             tot_x=float(np.sum(np.where(ok, x, 0) * w)), tot_y=float(np.sum(np.where(ok, y, 0) * w)), bands=[])
    for lo, hi in BANDS:
        s = used & (mid >= lo) & (mid < hi)
        n = int(s.sum()); c = float(np.sum(pulls[s] ** 2))
        r["bands"].append(dict(lo=lo, hi=hi, n=n, chi2=c, p=float(sst.chi2.sf(c, n)) if n else float("nan"),
                               mean=float(np.mean(pulls[s])) if n else float("nan")))
    return r


def band_label(lo, hi):
    return f"$<{hi:g}$ keV" if lo == 0 else (f"{lo:g} -- {hi:g} keV" if np.isfinite(hi) else f"$>{lo:g}$ keV")


def merge_flags(flag, gap=1):
    """Consecutive flagged bins (gaps of <= gap bins closed) -> list of (first, last+1)."""
    idx = np.where(flag)[0]
    if idx.size == 0:
        return []
    out = []; s = p = idx[0]
    for i in idx[1:]:
        if i - p - 1 <= gap:
            p = i
        else:
            out.append((s, p + 1)); s = p = i
    out.append((s, p + 1))
    return out


# ------------------------------------------------------------------ the reporter

class TGRS_LatexReporter:

    def __init__(self, root_a, root_b, tgrs_file, out_dir="latex_report", ref_label="Reference (A)",
                 test_label="Test (B)", short_a="A", short_b="B", lines_csv=None,
                 comp_labels=("Cosmic Photons", "Cosmic Protons Prompt", "Cosmic Protons Activation"),
                 activation_component=2, bin_tol=0.10, min_fraction=0.95, band_tol=0.05, sim_merge=1,
                 emin=50.0, fwhm0=0.9, fwhm1=0.0013, range_sigma=5.0, n_zoom=12, min_rate=0.0):
        self.ra, self.rb, self.tgrs = str(root_a), str(root_b), str(tgrs_file)
        self.out = Path(out_dir); self.fig = self.out / "figs"; self.fig.mkdir(parents=True, exist_ok=True)
        self.la, self.lb, self.sa, self.sb = ref_label, test_label, short_a, short_b
        self.lines_csv = lines_csv or str(HERE / "tgrs_lines.csv")
        self.comp_labels = list(comp_labels); self.act = activation_component
        self.bin_tol, self.min_fraction, self.band_tol = bin_tol, min_fraction, band_tol
        self.sim_merge, self.emin, self.fwhm0, self.fwhm1 = sim_merge, emin, fwhm0, fwhm1
        self.range_sigma, self.n_zoom, self.min_rate = range_sigma, n_zoom, min_rate
        self.cat = []
        ga, gb = re.search(r"Geant4 (\d+\.\d+)", ref_label), re.search(r"Geant4 (\d+\.\d+)", test_label)
        self.ga = f"Geant4 {ga.group(1)}" if ga else "A"
        self.gb = f"Geant4 {gb.group(1)}" if gb else "B"
        self._load()

    @staticmethod
    def pc(x, d=0):
        return f"{x * 100:.{d}f}\\%"

    # ---------------------------------------------------------------- data
    def _load(self):
        self.H = {}
        for tag, p in (("A", self.ra), ("B", self.rb)):
            self.H[tag] = root_hist(p, "Sim")
            self.H[tag + "c"] = []
            i = 0
            while True:
                try:
                    self.H[tag + "c"].append(root_hist(p, f"Component_{i}"))
                except KeyError:
                    break
                i += 1
        self.H["D"] = root_hist(self.ra, "TGRS")
        sa, sb = self.H["A"][1], self.H["B"][1]
        if len(sa) != len(sb) or not np.allclose(sa, sb):
            raise SystemExit("Sim histograms of A and B have different binning.")
        d_edges, s_edges = self.H["D"][1], sa
        lo, hi = max(d_edges[0], s_edges[0], self.emin), min(d_edges[-1], s_edges[-1])
        de = d_edges[(d_edges >= lo - 1e-9) & (d_edges <= SPLIT + 1e-9)]
        se = s_edges[(s_edges > de[-1] + 1e-9) & (s_edges <= hi + 1e-9)][::max(self.sim_merge, 1)]
        self.edges = np.concatenate([de, se])
        self.mid = 0.5 * (self.edges[:-1] + self.edges[1:]); self.w = np.diff(self.edges)
        g = {}
        for k in ("A", "B", "D"):
            v, e, er = self.H[k]
            g[k], g["e" + k] = remap_density(v, er, e, self.edges)
        self.g = g
        self.gc = {t: [remap_density(v, er, e, self.edges) for v, e, er in self.H[t + "c"]] for t in ("A", "B")}
        self.ncomp = len(self.H["Ac"])
        while len(self.comp_labels) < self.ncomp:
            self.comp_labels.append(f"Component {len(self.comp_labels)}")
        try:
            with open(self.lines_csv) as f:
                self.cat = [dict(E=float(r["energy_keV"]), id=r["id"], table=r["table"]) for r in csv.DictReader(f)]
        except OSError:
            self.cat = []

    def cat_lines(self, lo, hi, limit=3):
        sel = [c for c in self.cat if lo - 0.5 <= c["E"] <= hi + 0.5]
        sel.sort(key=lambda c: c["E"])
        txt = "; ".join(f"{c['E']:g} {nuc(LC.split_id(c['id'])[0])}" for c in sel[:limit])
        if len(sel) > limit:
            txt += f" (+{len(sel) - limit})"
        return txt or "--"

    # ---------------------------------------------------------------- figures
    def _save(self, fig, name):
        fig.savefig(self.fig / name, dpi=200, bbox_inches="tight"); plt.close(fig); return f"figs/{name}"

    def _xaxis(self, ax):
        ax.set_xscale("log"); ax.set_xlim(self.emin, self.edges[-1])

    def fig_spectrum(self, tag, name):
        g = self.g; fig, ax = plt.subplots(figsize=(6.4, 3.6))
        cols = ["gold", "yellowgreen", "tab:blue", "tab:purple", "tab:brown"]
        for i, (v, e) in enumerate(self.gc[tag]):
            ax.stairs(np.where(v > 0, v, np.nan), self.edges, color=cols[i % 5], lw=0.7, label=self.comp_labels[i])
        ax.stairs(np.where(g[tag] > 0, g[tag], np.nan), self.edges, color="tab:red", lw=0.8,
                  label=f"{self.la if tag == 'A' else self.lb} (combined)")
        ax.stairs(np.where(g["D"] > 0, g["D"], np.nan), self.edges, color="black", lw=0.9, label="Measurements")
        self._xaxis(ax); ax.set_yscale("log"); ax.set_ylim(1e-5, 3 * np.nanmax(g["D"]))
        ax.set_xlabel("Energy (keV)"); ax.set_ylabel("Rate (cts/keV/s)")
        ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.2), frameon=False)
        return self._save(fig, name)

    def fig_ratio(self, tag, name, merge=5):
        g = self.g; fig, ax = plt.subplots(figsize=(6.4, 2.4))
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where((g[tag] > 0) & (g["D"] > 0), g["D"] / g[tag], np.nan)
        idx = np.arange(len(r)) // merge; ng = idx.max() + 1
        fin = np.isfinite(r)
        cnt = np.bincount(idx, weights=fin, minlength=ng); tot = np.bincount(idx, weights=np.where(fin, r, 0), minlength=ng)
        m = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan); em = self.edges[::merge]
        if em[-1] != self.edges[-1]:
            em = np.append(em, self.edges[-1])
        ax.stairs(np.clip(m, 0.2, 5), em, color="tab:blue", lw=0.8)
        ax.axhline(1, color="gray", ls=":", lw=0.8, label="Perfect agreement")
        self._xaxis(ax); ax.set_yscale("log"); ax.set_ylim(0.2, 5)
        ax.set_yticks([0.2, 0.5, 1, 2, 5]); ax.set_yticklabels(["0.2", "0.5", "1", "2", "5"])
        ax.set_xlabel("Energy (keV)"); ax.set_ylabel("Measurements / simulation"); ax.legend(loc="upper right")
        return self._save(fig, name)

    def fig_combined(self, ranges, name):
        g = self.g
        fig, ax = plt.subplots(3, 1, figsize=(6.4, 7.4), sharex=True, gridspec_kw=dict(height_ratios=[3, 1, 1], hspace=0.06))
        for a in ax:
            for lo, hi, *_ in ranges[:self.n_zoom]:
                a.axvspan(lo, hi, color="orange", alpha=0.25, lw=0)
        ax[0].stairs(np.where(g["D"] > 0, g["D"], np.nan), self.edges, color="k", lw=0.8, label="Measurements")
        ax[0].stairs(np.where(g["A"] > 0, g["A"], np.nan), self.edges, color=COL["A"], lw=0.7, label=self.la)
        ax[0].stairs(np.where(g["B"] > 0, g["B"], np.nan), self.edges, color=COL["B"], lw=0.7, label=self.lb)
        ax[0].set_yscale("log"); ax[0].set_ylabel("Rate (cts/keV/s)")
        ax[0].set_ylim(1e-5, 3 * np.nanmax(g["D"]))
        ax[0].plot([], [], color="orange", alpha=0.4, lw=6, label="ranges zoomed in the section on differences")
        ax[0].legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.2), frameon=False)
        with np.errstate(divide="ignore", invalid="ignore"):
            for t in ("A", "B"):
                ax[1].stairs(np.where(g["D"] > 0, g[t] / g["D"], np.nan), self.edges, color=COL[t], lw=0.6)
            ax[2].stairs(np.where(g["B"] > 0, g["A"] / g["B"], np.nan), self.edges, color="purple", lw=0.6)
        ax[1].axhline(1, color="gray", lw=0.5); ax[2].axhline(1, color="gray", lw=0.5)
        ax[1].set_ylim(0.5, 1.5); ax[2].set_ylim(0.75, 1.25)
        ax[1].set_ylabel("Sim / data"); ax[2].set_ylabel("A / B")
        self._xaxis(ax[2]); ax[2].set_xlabel("Energy (keV)")
        return self._save(fig, name)

    def fig_residual(self, st, name):
        g = self.g; fig, ax = plt.subplots(2, 1, figsize=(6.4, 3.6), sharex=True, gridspec_kw=dict(hspace=0.06))
        d = g["A"] - g["B"]; sig = np.sqrt(g["eA"] ** 2 + g["eB"] ** 2)
        ax[0].fill_between(self.mid, -3 * sig, 3 * sig, step="mid", color="lightgray", label="$\\pm3\\sigma$ (statistics)")
        ax[0].stairs(d, self.edges, color="purple", lw=0.6, label="A $-$ B")
        ax[0].set_yscale("symlog", linthresh=1e-4); ax[0].set_ylabel("Rate difference (cts/keV/s)"); ax[0].legend(loc="upper right")
        ax[1].stairs(st["pulls"], self.edges, color="purple", lw=0.6)
        for y in (-3, 3):
            ax[1].axhline(y, color="red", lw=0.5)
        ax[1].set_ylabel("Pull (A $-$ B)/$\\sigma$"); self._xaxis(ax[1]); ax[1].set_xlabel("Energy (keV)")
        return self._save(fig, name)

    def fig_component(self, i, name):
        (va, ea, ed), (vb, eb, _) = self.H["Ac"][i][0:3], self.H["Bc"][i][0:3]
        va, ed1, _e = self.H["Ac"][i]; vb = self.H["Bc"][i][0]
        fig, ax = plt.subplots(2, 1, figsize=(6.4, 2.0), sharex=True, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05))
        ax[0].stairs(np.where(va > 0, va, np.nan), ed1, color=COL["A"], lw=0.6, label=self.la)
        ax[0].stairs(np.where(vb > 0, vb, np.nan), ed1, color=COL["B"], lw=0.6, label=self.lb)
        ax[0].set_yscale("log"); ax[0].set_ylabel("cts/keV/s"); ax[0].set_title(self.comp_labels[i], fontsize=8)
        ax[0].legend(loc="upper right")
        with np.errstate(divide="ignore", invalid="ignore"):
            ax[1].stairs(np.where(vb > 0, va / vb, np.nan), ed1, color="purple", lw=0.6)
        ax[1].axhline(1, color="gray", lw=0.5); ax[1].set_ylim(0.5, 1.5); ax[1].set_ylabel("A / B")
        ax[1].set_xscale("log"); ax[1].set_xlim(self.emin, ed1[-1]); ax[1].set_xlabel("Energy (keV)")
        return self._save(fig, name)

    def fig_zoom(self, rng_list, name, title_fn):
        """Panels of energy ranges: data, A, B with 1-sigma bands, dashed = component changing most."""
        n = len(rng_list); nr = math.ceil(n / 2)
        fig = plt.figure(figsize=(6.6, 2.7 * nr))
        gs = fig.add_gridspec(nr * 3, 2, height_ratios=[3, 1, 1.1] * nr, hspace=0.08, wspace=0.28)
        g = self.g
        for k, r in enumerate(rng_list):
            lo, hi = r[0], r[1]; wd = hi - lo; x0, x1 = lo - 1.2 * wd, hi + 1.2 * wd
            a0 = fig.add_subplot(gs[3 * (k // 2), k % 2]); a1 = fig.add_subplot(gs[3 * (k // 2) + 1, k % 2], sharex=a0)
            sel = (self.mid > x0) & (self.mid < x1)
            for key, c, lab in (("D", "k", "data"), ("A", COL["A"], "A"), ("B", COL["B"], "B")):
                a0.stairs(g[key][sel], self.edges[np.append(sel, False) | np.append(False, sel)][:sel.sum() + 1] if False else self.edges[np.where(sel)[0][0]:np.where(sel)[0][-1] + 2],
                          color=c, lw=0.8, label=lab)
                if key != "D":
                    a0.fill_between(self.mid[sel], (g[key] - g["e" + key])[sel], (g[key] + g["e" + key])[sel], step="mid",
                                    color=c, alpha=0.2, lw=0)
            # component that changes most
            best, bd = None, 0
            for i in range(self.ncomp):
                d = abs(integ(self.gc["A"][i][0], self.gc["A"][i][1], self.edges, lo, hi)[0]
                        - integ(self.gc["B"][i][0], self.gc["B"][i][1], self.edges, lo, hi)[0])
                if d > bd:
                    best, bd = i, d
            if best is not None:
                for t, c in (("A", COL["A"]), ("B", COL["B"])):
                    a0.plot(self.mid[sel], self.gc[t][best][0][sel], color=c, ls="--", lw=0.6, alpha=0.8)
            a0.axvspan(lo, hi, color="orange", alpha=0.25, lw=0)
            a0.set_title(title_fn(r, best), fontsize=6); a0.tick_params(labelbottom=False)
            if k == 0:
                a0.legend(fontsize=5, loc="upper right")
            with np.errstate(divide="ignore", invalid="ignore"):
                a1.stairs((g["A"] / g["B"])[sel], self.edges[np.where(sel)[0][0]:np.where(sel)[0][-1] + 2], color="purple", lw=0.7)
            a1.axhline(1, color="gray", lw=0.5); a1.axvspan(lo, hi, color="orange", alpha=0.25, lw=0)
            a1.set_ylabel("A/B", fontsize=6); a1.set_xlim(x0, x1); a1.set_ylim(0.6, 1.4)
            if k >= n - 2:
                a1.set_xlabel("Energy (keV)", fontsize=6)
        return self._save(fig, name)

    def fig_lines(self, rows, name1, name2):
        E = np.array([r["energy"] for r in rows])
        fig, ax = plt.subplots(2, 1, figsize=(6.4, 5.0), sharex=True, gridspec_kw=dict(hspace=0.08))
        for a, t in zip(ax, ("A", "B")):
            for tb, mk in (("5", "o"), ("6", "s")):
                s = np.array([r["table"] == tb for r in rows])
                v = np.array([r["ratio_" + t + "_data"] for r in rows]); ul = np.array([r["ul_" + t] for r in rows])
                for flag, fc in ((False, COL[t]), (True, "none")):
                    m = s & (ul == flag) & np.isfinite(v) & (v > 0)
                    a.scatter(E[m], v[m], s=7, marker=mk, facecolors=fc, edgecolors=COL[t], lw=0.6,
                              label=f"Table {tb}{', upper limit' if flag else ''}")
            a.axhline(1, color="k", ls=":", lw=0.7); a.set_yscale("log"); a.set_xscale("log")
            a.set_ylabel(f"Sim / obs\n{self.la if t == 'A' else self.lb}", fontsize=7); a.legend(fontsize=5, ncol=2)
        ax[1].set_xlabel("Line energy (keV)")
        p1 = self._save(fig, name1)
        fig, ax = plt.subplots(2, 1, figsize=(6.4, 5.0), sharex=True, gridspec_kw=dict(hspace=0.08))
        det = [r for r in rows if r["det_A"] and r["det_B"] and r["A"] > 0 and r["B"] > 0]
        if det:
            ax[0].errorbar([r["energy"] for r in det], [r["ratio_AB"] for r in det], [r["eratio_AB"] for r in det],
                           fmt="o", ms=2.5, color="purple", lw=0.6)
        ax[0].axhline(1, color="k", ls=":", lw=0.7); ax[0].set_yscale("log"); ax[0].set_xscale("log")
        ax[0].set_ylabel("A / B (detected in both)")
        ax[1].scatter(E, [r["sig_AB"] for r in rows], s=6, color="purple")
        for y in (-3, 3):
            ax[1].axhline(y, color="red", lw=0.5)
        ax[1].set_ylabel("(A $-$ B)/$\\sigma$"); ax[1].set_xscale("log"); ax[1].set_xlabel("Line energy (keV)")
        return p1, self._save(fig, name2)

    def fig_forest(self, summ, name):
        s = [x for x in summ if x["n"] >= 3]
        fig, ax = plt.subplots(figsize=(6.4, max(2.0, 0.22 * len(s) + 0.8)))
        for i, x in enumerate(s):
            ax.errorbar(x["median_B_over_A"], i, xerr=[[x["median_B_over_A"] - x["ci_low"]], [x["ci_high"] - x["median_B_over_A"]]],
                        fmt="o", ms=3, color="tab:red" if x["differs"] else "gray", capsize=2, lw=0.8)
        ax.axvline(1, color="k", lw=0.7); ax.set_yticks(range(len(s)))
        ax.set_yticklabels([f"{x['group']} (n={x['n']})" for x in s], fontsize=6); ax.invert_yaxis()
        ax.set_xscale("log"); ax.set_xlabel("Median B/A with 95% bootstrap interval (red: interval excludes 1)")
        return self._save(fig, name)

    # ---------------------------------------------------------------- analysis
    def ranges_ab(self):
        g = self.g; st = self.st_ab
        with np.errstate(invalid="ignore"):
            flag = (np.abs(st["pulls"]) > 3) & (np.abs(g["A"] - g["B"]) > 0.03 * np.maximum(g["A"], g["B"]))
        out = []
        for i0, i1 in merge_flags(np.nan_to_num(flag, nan=0).astype(bool)):
            lo, hi = self.edges[i0], self.edges[i1]
            a, ea = integ(g["A"], g["eA"], self.edges, lo, hi); b, eb = integ(g["B"], g["eB"], self.edges, lo, hi)
            d, ed = integ(g["D"], g["eD"], self.edges, lo, hi)
            s = (a - b) / math.hypot(ea, eb) if (ea or eb) else 0.0
            if abs(s) > self.range_sigma:
                out.append((lo, hi, d, a, ea, b, eb, s))
        out.sort(key=lambda r: -abs(r[7]))
        return out

    def ranges_data(self):
        g = self.g; fl = {}
        for t in ("A", "B"):
            with np.errstate(invalid="ignore", divide="ignore"):
                sg = np.sqrt(g["e" + t] ** 2 + g["eD"] ** 2)
                fl[t] = np.nan_to_num((np.abs(g[t] - g["D"]) > 0.2 * g["D"]) & (np.abs(g[t] - g["D"]) > 3 * sg), nan=0).astype(bool)
        out = []
        for i0, i1 in merge_flags(fl["A"] | fl["B"]):
            lo, hi = self.edges[i0], self.edges[i1]
            d, _ = integ(g["D"], g["eD"], self.edges, lo, hi); a, _ = integ(g["A"], g["eA"], self.edges, lo, hi)
            b, _ = integ(g["B"], g["eB"], self.edges, lo, hi)
            if d <= 0:
                continue
            ba, bb = abs(a / d - 1) > 0.15, abs(b / d - 1) > 0.15
            if ba or bb:
                out.append((lo, hi, d, a, b, "both" if ba and bb else ("A" if ba else "B")))
        out.sort(key=lambda r: -max(abs(r[3] - r[2]), abs(r[4] - r[2])))
        return out

    def tolerance(self):
        g = self.g; a, b, ea, eb = g["A"], g["B"], g["eA"], g["eB"]
        sig = np.sqrt(ea ** 2 + eb ** 2); ok = np.isfinite(a) & np.isfinite(b)
        used = ok & (sig > 0) & ~((a == 0) & (b == 0))
        good = np.abs(a - b) <= self.bin_tol * np.maximum(np.abs(a), np.abs(b)) + 3 * sig
        frac = float(np.mean(good[used]))
        rows = []; all_ok = True
        for lo, hi in BANDS:
            s = ok & (self.mid >= lo) & (self.mid < hi)
            ia, sa = integ(np.where(s, a, 0), np.where(s, ea, 0), self.edges, 0, 1e9)
            ib, sb = integ(np.where(s, b, 0), np.where(s, eb, 0), self.edges, 0, 1e9)
            gd = abs(ia - ib) <= self.band_tol * max(abs(ia), abs(ib)) + 3 * math.hypot(sa, sb)
            all_ok &= gd; rows.append((band_label(lo, hi), ia / ib if ib else float("nan"), gd))
        self.tol = dict(frac=frac, bands=rows, passed=bool(frac >= self.min_fraction and all_ok))
        return self.tol

    def line_rows(self):
        H = {"A": self.H["A"], "B": self.H["B"], "D": self.H["D"]}
        if 0 <= self.act < self.ncomp:
            H["Ac"], H["Bc"] = self.H["Ac"][self.act], self.H["Bc"][self.act]
        rows = []
        for r in self.cat:
            e0 = r["E"]; fw = self.fwhm0 + self.fwhm1 * e0
            res = {k: LC.net_rate(v, er, ed, e0, fw) for k, (v, ed, er) in H.items()}
            if any(x is None for x in res.values()):
                continue
            (nA, eA, _), (nB, eB, _), (nD, eD, _) = res["A"], res["B"], res["D"]
            paper = next((float(c["rate_per_s"]) for c in self._raw if float(c["energy_keV"]) == e0 and c["id"] == r["id"]), float("nan"))
            d = nA - nB; sd = math.hypot(eA, eB)
            rAB, erAB = LC.ratio(nA, eA, nB, eB); rAD, _ = LC.ratio(nA, eA, nD, eD); rBD, _ = LC.ratio(nB, eB, nD, eD)
            row = dict(table=r["table"], energy=e0, paper=paper, id=r["id"], A=nA, eA=eA, B=nB, eB=eB, data=nD, edata=eD,
                       diff=d, sig_AB=d / sd if sd > 0 else 0.0, ratio_AB=rAB, eratio_AB=erAB,
                       ratio_A_data=rAD, ratio_B_data=rBD,
                       det_A=nA > 2 * eA, det_B=nB > 2 * eB, det_D=nD > 2 * eD)
            row["ul_A"], row["ul_B"], row["ul_D"] = not row["det_A"], not row["det_B"], not row["det_D"]
            row["ulv_A"], row["ulv_B"], row["ulv_D"] = (max(nA, 0) + 2 * eA, max(nB, 0) + 2 * eB, max(nD, 0) + 2 * eD)
            if "Ac" in res:
                row.update(cA=res["Ac"][0], ecA=res["Ac"][1], cB=res["Bc"][0], ecB=res["Bc"][1])
                row["det_cA"], row["det_cB"] = row["cA"] > 2 * row["ecA"], row["cB"] > 2 * row["ecB"]
                row["comp_sig"] = (row["cA"] - row["cB"]) / max(math.hypot(row["ecA"], row["ecB"]), 1e-30)
            rows.append(row)
        return rows

    # ---------------------------------------------------------------- formatting of line values
    @staticmethod
    def rate_pm(r, k):
        v, e = (r["data"], r["edata"]) if k == "D" else (r[k], r["e" + k])
        return pm(v, e) if r["det_" + k] else f"$<${r['ulv_' + k]:.2g}"

    @staticmethod
    def rate_paren(r, k):
        v, e = (r["data"], r["edata"]) if k == "D" else (r[k], r["e" + k])
        return paren(v, e) if r["det_" + k] else f"$<${r['ulv_' + k]:.2g}"

    def sim_over_obs(self, r, k, paren_style=True):
        """Sim/Obs with upper limits, as in the paper."""
        if r["det_" + k] and r["det_D"] and r["data"] > 0:
            v = r[k] / r["data"]; e = abs(v) * math.hypot(r["e" + k] / r[k] if r[k] else 0, r["edata"] / r["data"])
            return paren(v, e) if paren_style else pm(v, e)
        if (not r["det_" + k]) and r["det_D"] and r["data"] > 0:
            return f"$<${r['ulv_' + k] / r['data']:.2g}"
        return "$\\cdots$"

    def a_over_b(self, r, paren_style=True):
        if r["det_A"] and r["det_B"] and r["B"] > 0:
            return paren(r["ratio_AB"], r["eratio_AB"]) if paren_style else pm(r["ratio_AB"], r["eratio_AB"])
        if not r["det_A"] and r["det_B"]:
            return f"$<${r['ulv_A'] / r['B']:.2g}"
        if r["det_A"] and not r["det_B"] and r["ulv_B"] > 0:
            return f"$>${r['A'] / r['ulv_B']:.2g}"
        return "$\\cdots$"

    def differs_ab(self, r):
        return abs(r["sig_AB"]) > 3 and abs(r["A"] - r["B"]) > 0.1 * max(abs(r["A"]), abs(r["B"]), 1e-30)

    def differs_data(self, r, k):
        if not r["det_D"] or r["data"] <= 0:
            return False
        if not r["det_" + k]:
            return r["ulv_" + k] / r["data"] < 0.5
        v = r[k] / r["data"]
        sg = (r[k] - r["data"]) / math.hypot(r["e" + k], r["edata"])
        return abs(v - 1) > 0.3 and abs(sg) > 3

    # ---------------------------------------------------------------- build
    def build(self, compile_pdf=True, title="TGRS MEGAlib Validation Report"):
        with open(self.lines_csv) as f:
            self._raw = list(csv.DictReader(f)) if Path(self.lines_csv).exists() else []
        g, ed = self.g, self.edges
        okb = np.isfinite(g["A"]) & np.isfinite(g["B"]) & np.isfinite(g["D"])
        self.st_ab = gof(g["A"], g["eA"], g["B"], g["eB"], ed)
        self.st_am = gof(g["A"], g["eA"], g["D"], g["eD"], ed); self.st_bm = gof(g["B"], g["eB"], g["D"], g["eD"], ed)
        tol = self.tolerance()
        R8, R9 = self.ranges_ab(), self.ranges_data()
        rows = self.line_rows()
        T = []
        ap = T.append
        ap(self._preamble(title))

        # 1 Configuration
        ap("\\section{Configuration Summary}\n")
        ap(guide("Which two builds are compared and with which settings. A is the reference build, B the test build. "
                 "Only the A-versus-B tolerance test decides PASS or FAIL; every comparison with the measurements is information only."))
        cfg = [("Reference (A)", esc(self.la)), ("Test (B)", esc(self.lb)), ("Reference file", "\\texttt{" + esc(self.ra) + "}"),
               ("Test file", "\\texttt{" + esc(self.rb) + "}"), ("TGRS data file", "\\texttt{" + esc(self.tgrs) + "}"),
               ("Measurement scaling", "$\\times35$ cm$^2$ (as in TGRS.cxx)"),
               ("Energy range (keV)", f"{ed[0]:.2f} -- {ed[-1]:.2f}"),
               ("Comparison bins", f"{len(self.mid)} (data bins below {SPLIT:g} keV, simulation bins above)"),
               ("Simulation components", esc(", ".join(self.comp_labels[:self.ncomp]))),
               ("Line comparison", f"net rate in $\\pm1.5$ FWHM minus the sideband continuum (FWHM = {self.fwhm0} + {self.fwhm1}$\\,E$ keV); "
                                   "upper limit = net + 2$\\sigma$"),
               ("Overall PASS/FAIL decided by",
                f"A versus B tolerance test: at least {self.pc(self.min_fraction)} of the bins within {self.pc(self.bin_tol)} (or 3$\\sigma$) and each band "
                f"($<100$, 100--1000, $>1000$ keV) within {self.pc(self.band_tol)} (or 3$\\sigma$); thresholds await approval")]
        ap("\\begin{center}\\small\\begin{tabularx}{\\linewidth}{lX}\\toprule Parameter & Value\\\\\\midrule\n"
           + "\\\\\n".join(f"{a} & {b}" for a, b in cfg) + "\\\\\\bottomrule\\end{tabularx}\\end{center}\n")

        # 2-3 sim vs measurement ; 4-5 ratios
        for n, (tag, st, lab) in enumerate((("A", self.st_am, self.la), ("B", self.st_bm, self.lb))):
            ap(f"\\section{{Sim {tag} vs Measurements (information only)}}\n")
            ap(guide("Spectrum of the simulation (red) with its components and the measurements (black), then the goodness of fit. "
                     "Pull = (simulation $-$ data)/uncertainty per bin. $\\chi^2/\\nu$ near 1 would mean agreement within statistics; with this many simulated counts "
                     "it is large even for small differences, so it measures the size of the differences, not a verdict."))
            ap(f"\\begin{{center}}\\includegraphics[width=0.92\\linewidth]{{{self.fig_spectrum(tag, f'spec_{tag}.png')}}}\\end{{center}}\n")
            ap(self._gof_table(st, esc(lab), "Measurements"))
        for tag, st in (("A", self.st_am), ("B", self.st_bm)):
            lab = self.la if tag == "A" else self.lb
            ap(f"\\section{{Ratio {tag} vs Measurements (information only)}}\n")
            ap(guide("Measurements divided by the simulation. 1 (dotted) is perfect agreement; above 1 the simulation is lower than the data. "
                     "Narrow spikes are strong lines (the measured line shape has a low-energy tail, the simulated one is Gaussian)."))
            ap(f"\\begin{{center}}\\includegraphics[width=0.92\\linewidth]{{{self.fig_ratio(tag, f'ratio_{tag}.png')}}}\\end{{center}}\n")
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where((g[tag] > 0) & (g["D"] > 0), g["D"] / g[tag], np.nan)
            rr = r[np.isfinite(r)]
            ap(table([["Bins with a ratio", len(rr)], ["Median ratio", sci(float(np.median(rr)), 4)],
                      ["Mean ratio", sci(float(np.mean(rr)), 4)], ["Fraction of bins within 20\\% of 1", sci(float(np.mean(np.abs(rr - 1) < 0.2)), 3)]],
                     ["Quantity", "Value"], "lr"))

        # 6 A vs B
        ap("\\section{Sim A vs Sim B (combined)}\n")
        ap(guide("Both builds and the measurements on one energy axis. Top: spectra. Middle: each build divided by the data (1 = agreement). "
                 "Bottom: A/B, the quantity the tolerance test judges. Orange bands mark the ranges zoomed in the section on where A and B differ. "
                 "On the next pages: the residuals (difference and pull per bin), the metrics, and the tolerance test that decides PASS or FAIL."))
        ap(f"\\begin{{center}}\\includegraphics[width=0.9\\linewidth]{{{self.fig_combined(R8, 'combined.png')}}}\\end{{center}}\n")
        ap("\\clearpage\n")
        ap(f"\\begin{{center}}\\includegraphics[width=0.82\\linewidth]{{{self.fig_residual(self.st_ab, 'residual.png')}}}\\end{{center}}\n")
        ap(self._gof_table(self.st_ab, esc(self.la), esc(self.lb)))
        ap("\\subsection*{Tolerance test (decides PASS/FAIL)}\n")
        ap(guide(f"Two criteria: (1) at least {self.pc(self.min_fraction)} of the bins within {self.pc(self.bin_tol)} (or 3$\\sigma$) of each other; "
                 f"(2) in each band the integrated A/B within {self.pc(self.band_tol)} (or 3$\\sigma$). Each line below shows the measured value and whether it passes."))
        rows_t = [["Bins within tolerance", sci(tol["frac"], 4), f"required $\\ge$ {self.min_fraction}", "yes" if tol["frac"] >= self.min_fraction else "no"]]
        for lab, ratio, gd in tol["bands"]:
            rows_t.append([f"Integrated rate {lab}", f"A/B = {sci(ratio, 4)}", f"within {self.pc(self.band_tol)}", "yes" if gd else "no"])
        ap(table(rows_t, None, "llll"))
        ap(f"\\noindent Status: \\textbf{{\\textcolor{{{'green!50!black' if tol['passed'] else 'red'}}}{'PASS' if tol['passed'] else 'FAIL'}}}\n")

        # 7 Components
        ap("\\section{Components, build against build (information only)}\n")
        ap(guide("The simulation has components (photons, prompt protons, activation). Each is plotted for A and B with their ratio A/B, and the integrated rates per energy band "
                 "are tabulated. A ratio different from 1 in one component shows which part of the simulation causes a difference in the total."))
        for i in range(self.ncomp):
            ap(f"\\begin{{center}}\\includegraphics[width=0.8\\linewidth]{{{self.fig_component(i, f'comp_{i}.png')}}}\\end{{center}}\n")
        hdr = ["Energy band", "Data (s$^{-1}$)", "A (s$^{-1}$)", "B (s$^{-1}$)", "A/Data", "B/Data", "A/B"]
        rws = []
        for lo, hi in COMP_BANDS:
            d, a, b = (integ(g[k], g["e" + k], ed, lo, hi)[0] for k in ("D", "A", "B"))
            rws.append([f"{lo} -- {hi} keV", sci(d, 4), sci(a, 4), sci(b, 4), sci(a / d, 4), sci(b / d, 4), sci(a / b, 4)])
        ap("\\subsection*{Integrated rates, total}\n" + table(rws, hdr, "lrrrrrr"))
        rws = []
        for i in range(self.ncomp):
            for lo, hi in COMP_BANDS:
                a = integ(self.gc["A"][i][0], self.gc["A"][i][1], ed, lo, hi)[0]; b = integ(self.gc["B"][i][0], self.gc["B"][i][1], ed, lo, hi)[0]
                rws.append([esc(self.comp_labels[i]) if lo == COMP_BANDS[0][0] else "", f"{lo} -- {hi} keV", sci(a, 4), sci(b, 4), sci(a / b, 4) if b else "--"])
        ap("\\subsection*{Integrated rates per component}\n" + table(rws, ["Component", "Energy band", "A (s$^{-1}$)", "B (s$^{-1}$)", "A/B"], "llrrr"))

        # 8 Where A and B differ
        ap("\\section{Where A and B differ (information only)}\n")
        ap(guide(f"Energy ranges where the builds disagree bin by bin ($|$pull$|>3$ and more than 3\\% apart), merged into ranges and kept when the integrated difference exceeds "
                 f"{self.range_sigma:g}$\\sigma$. Columns: rates in counts/s (continuum included), A/B, significance $\\sigma$ of A$-$B (negative: B higher), ratios to the data and the catalogue "
                 f"lines in the range. {len(R8)} ranges, largest first; the {min(self.n_zoom, len(R8))} largest are drawn after the table."))
        rws = [[f"{r[0]:.1f}--{r[1]:.1f}", sci(r[2], 3), pm(r[3], r[4]), pm(r[5], r[6]), sci(r[3] / r[5], 3), f"{r[7]:+.0f}".replace("-", "$-$"),
                sci(r[3] / r[2], 3) if r[2] else "--", sci(r[5] / r[2], 3) if r[2] else "--", self.cat_lines(r[0], r[1])] for r in R8]
        ap(self._long(rws, ["Range (keV)", "Data", "A", "B", "A/B", "$\\sigma$", "A/Data", "B/Data", "Catalogue lines"], "lrrrrrrrp{3.3cm}"))
        top = R8[:self.n_zoom]
        for chunk, nm in ((top[:6], "zoom_ab_1.png"), (top[6:12], "zoom_ab_2.png")):
            if chunk:
                ap("\\clearpage\n\\begin{center}\\includegraphics[width=0.98\\linewidth]{" + self.fig_zoom(
                    chunk, nm, lambda r, b: f"{r[0]:.0f}--{r[1]:.0f} keV, A/B={r[3] / r[5]:.2f} ({r[7]:+.0f}$\\sigma$)".replace("$\\sigma$", " sigma").replace("-", "-")) + "}\\end{center}\n")
        if top:
            ap(guide("Zoom panels: data black, A blue, B red with 1$\\sigma$ bands; the dashed curves are the component (photons, prompt or activation) that changes most in the range, "
                     "which points to the origin of the difference; the lower panel is A/B."))

        # 9 Where each build differs from the data
        ap("\\section{Where each build differs from the data (information only)}\n")
        ap(guide("Here each build is compared with the measurements. Ranges are listed where the simulation is more than 20\\% away from the data bin by bin (beyond 3$\\sigma$) and the "
                 "integrated rate is more than 15\\% off. Blocks: \\textbf{both builds} (a shared difference, not a regression), \\textbf{only " + self.ga + "}, \\textbf{only " + self.gb + "} (the builds behave differently). "
                 "Ranges at strong lines can partly reflect the line shape of the data."))
        for cat, name in (("both", "Both builds differ from the data"), ("A", f"Only {self.ga} differs"), ("B", f"Only {self.gb} differs")):
            sel = [r for r in R9 if r[5] == cat]
            rws = [[f"{r[0]:.1f}--{r[1]:.1f}", sci(r[2], 3), sci(r[3], 3), sci(r[4], 3), sci(r[3] / r[2], 3), sci(r[4] / r[2], 3), self.cat_lines(r[0], r[1])] for r in sel]
            ap(f"\\subsection*{{{name} ({len(sel)})}}\n" + (self._long(rws, ["Range (keV)", "Data", "A", "B", "A/Data", "B/Data", "Catalogue lines"], "lrrrrrp{4.2cm}") if rws else "None.\n"))
        top9 = R9[:self.n_zoom]
        for chunk, nm in ((top9[:6], "zoom_d_1.png"), (top9[6:12], "zoom_d_2.png")):
            if chunk:
                ap("\\clearpage\n\\begin{center}\\includegraphics[width=0.98\\linewidth]{" + self.fig_zoom(
                    chunk, nm, lambda r, b: f"{r[0]:.0f}--{r[1]:.0f} keV, A/D={r[3] / r[2]:.2f}, B/D={r[4] / r[2]:.2f}") + "}\\end{center}\n")

        # 10 Lines
        ap("\\section{Lines, Tables 5 and 6 (information only)}\n")
        ap(self._lines_section(rows))

        # 11 Overall
        ap("\\section{Overall Validation Status}\n")
        ap(guide("PASS or FAIL comes only from the A-versus-B tolerance test of the section on A versus B. The comparisons with the data, the chi-squared values and the line tables "
                 "locate and describe differences; they do not change the status."))
        ap(f"\\noindent\\textbf{{\\textcolor{{{'green!50!black' if tol['passed'] else 'red'}}}{'PASS' if tol['passed'] else 'FAIL'}}}\\\\[2pt]\n"
           f"Decided by the A versus B tolerance test: at least {self.pc(self.min_fraction)} of the bins within {self.pc(self.bin_tol)} (or 3$\\sigma$) and each band within {self.pc(self.band_tol)} (or 3$\\sigma$). "
           f"Measured: {self.pc(tol['frac'],1)} of the bins; band ratios A/B: " + ", ".join(f"{sci(r, 4)}" for _, r, _ in tol["bands"]) + ".\n")

        # Appendix A (all lines, paper layout) and B (group patterns)
        ap(self._appendix_lines(rows))
        ap(self._appendix_groups(rows))
        ap("\\end{document}\n")
        tex = self.out / "TGRS_Validation_Report.tex"
        tex.write_text("\n".join(T), encoding="utf-8")
        print(f"Wrote {tex}")
        if compile_pdf:
            return self.compile(tex)
        return tex

    # ---------------------------------------------------------------- pieces
    def _preamble(self, title):
        return (r"""\documentclass[11pt]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage{lmodern}
\usepackage[letterpaper,margin=0.9in]{geometry}
\usepackage{booktabs,longtable,graphicx,pdflscape,amsmath,xcolor,array,tabularx,hyperref,fancyvrb}
\hypersetup{colorlinks=true,linkcolor=blue!60!black,urlcolor=blue!60!black,pdftitle={""" + title + r"""}}
\usepackage{titlesec}\newcommand{\sectionbreak}{\clearpage}\titleformat{\section}{\Large\bfseries}{\thesection}{1em}{}
\newcommand{\guide}[1]{\begin{center}\fcolorbox{blue!25!gray}{blue!4}{\parbox{0.95\linewidth}{\small\textbf{How to read this.} #1}}\end{center}}
\setlength{\tabcolsep}{4pt}
\begin{document}
\begin{center}{\LARGE\bfseries """ + title + r"""}\\[4pt]
A = """ + esc(self.la) + r""" \quad B = """ + esc(self.lb) + r"""\\[2pt]\today\end{center}
\renewcommand{\sectionbreak}{}\tableofcontents\clearpage
\renewcommand{\sectionbreak}{\clearpage}
""").replace("<A>", "A").replace("<B>", "B")

    def _gof_table(self, st, nx, ny):
        rows = [["\\textbf{Spectra compared}", ""], ["Bins compared", st["ndf"]], [f"Total rate, {nx} (cts/s)", sci(st["tot_x"], 4)],
                [f"Total rate, {ny} (cts/s)", sci(st["tot_y"], 4)], ["\\textbf{Goodness of fit}", ""], ["$\\chi^2$", sci(st["chi2"], 4)],
                ["Degrees of freedom $\\nu$", st["ndf"]], ["$\\chi^2/\\nu$", sci(st["chi2"] / st["ndf"], 4)], ["$p$-value", pval(st["p"])],
                ["Bins beyond $3\\sigma$", st["n3"]], ["Worst bin (keV)", f"{st['worst'][0]:.1f} -- {st['worst'][1]:.1f}"],
                ["Pull of the worst bin", sci(st["worst_pull"], 4).replace("-", "$-$")], ["Kolmogorov--Smirnov distance (cumulative spectra)", sci(st["ks"], 4)]]
        t = table(rows, ["Quantity", "Value"], "lr")
        br = [[band_label(b["lo"], b["hi"]), b["n"], sci(b["chi2"], 4), b["n"], sci(b["chi2"] / b["n"], 4) if b["n"] else "--", pval(b["p"]),
               sci(b["mean"], 3).replace("-", "$-$")] for b in st["bands"]]
        return t + table(br, ["Energy band", "Bins", "$\\chi^2$", "$\\nu$", "$\\chi^2/\\nu$", "$p$-value", "Mean pull"], "lrrrrrr")

    @staticmethod
    def _long(rows, head, spec):
        out = [f"{{\\scriptsize\\setlength{{\\tabcolsep}}{{3pt}}\\begin{{longtable}}{{{spec}}}\\toprule", " & ".join(head) + "\\\\\\midrule\\endfirsthead",
               "\\toprule " + " & ".join(head) + "\\\\\\midrule\\endhead", "\\bottomrule\\endfoot"]
        out += [" & ".join(str(c) for c in r) + "\\\\" for r in rows]
        out.append("\\end{longtable}}\n")
        return "\n".join(out)

    def _lines_section(self, rows):
        t5 = [r for r in rows if r["table"] == "5"]; t6 = [r for r in rows if r["table"] == "6"]
        d_ab = [r for r in rows if r["det_A"] and r["det_B"]]
        flagged = sorted([r for r in rows if self.differs_ab(r)], key=lambda r: -abs(r["sig_AB"]))
        med = lambda L, k: (float(np.median([r[k] / r["data"] for r in L if r["det_" + k] and r["det_D"] and r["data"] > 0] or [np.nan])))
        fw1000 = self.fwhm0 + self.fwhm1 * 1000
        s = [guide("Every line of Tables 5 and 6 of Weidenspointner et al. (2005) is fitted in the data and in both simulations with the same method: window $\\pm1.5$ FWHM "
                   "minus the continuum from two sidebands. Rates in counts/s; $<$ is an upper limit (net + 2$\\sigma$, line not detected at 2$\\sigma$); ratios Sim/Data as in the paper "
                   "(1 = the simulation reproduces the observed rate, below 1 it is too low); A/B compares the builds. Appendix A lists all lines.")]
        sm = [["Lines compared", len(rows)], ["Table 5 / Table 6 lines", f"{len(t5)} / {len(t6)}"], ["Detector FWHM at 1 MeV (keV)", f"{fw1000:.3g}"],
              ["Lines detected in the data ($>2\\sigma$)", sum(r["det_D"] for r in rows)],
              [f"Lines detected in A", sum(r["det_A"] for r in rows)], [f"Lines detected in B", sum(r["det_B"] for r in rows)],
              ["Median Sim/Data, A", sci(med(rows, "A"), 3)], ["Median Sim/Data, B", sci(med(rows, "B"), 3)],
              [f"Median A/B ({len(d_ab)} lines detected in both)", sci(float(np.median([r["A"] / r["B"] for r in d_ab])) if d_ab else float("nan"), 3)],
              ["Lines with $|$A$-$B$|>3\\sigma$ and $>10\\%$", len(flagged)]]
        s.append(table(sm, ["Quantity", "Value"], "lr"))
        p1, p2 = self.fig_lines(rows, "lines_simobs.png", "lines_ab.png")
        s.append(f"\\begin{{center}}\\includegraphics[width=0.88\\linewidth]{{{p1}}}\\end{{center}}\n\\begin{{center}}\\includegraphics[width=0.88\\linewidth]{{{p2}}}\\end{{center}}\n")
        s.append("\\subsection*{Lines where A and B differ}\n")
        has_c = "cA" in rows[0] if rows else False
        s.append(guide("Lines whose rates differ by more than 3$\\sigma$ and more than 10\\%, largest difference first." +
                       (" The last two columns give the same line fitted in the activation component alone: when they carry the difference, it comes from isotope production or decay, not from photons or prompt protons." if has_c else "")))
        rws = []
        for r in flagged:
            row = [f"{r['energy']:.1f}", nuc(LC.split_id(r["id"])[0]), self.rate_pm(r, "A"), self.rate_pm(r, "B"), self.a_over_b(r, False), sigtxt(r["sig_AB"])]
            if has_c:
                row += [pm(r["cA"], r["ecA"]) if r["det_cA"] else f"$<${max(r['cA'], 0) + 2 * r['ecA']:.2g}",
                        pm(r["cB"], r["ecB"]) if r["det_cB"] else f"$<${max(r['cB'], 0) + 2 * r['ecB']:.2g}"]
            rws.append(row)
        s.append(self._long(rws, ["E (keV)", "Transition", "A (s$^{-1}$)", "B (s$^{-1}$)", "A/B", "(A$-$B)/$\\sigma$"] + (["A (act.)", "B (act.)"] if has_c else []),
                            "lp{2.6cm}rrrr" + ("rr" if has_c else "")))
        s.append("\\subsection*{Lines where each build differs from the data}\n")
        s.append(guide(f"A build differs from the data on a line if the line is detected in the data and the simulation shows it with a rate more than 30\\% off (and beyond 3$\\sigma$), "
                       "or does not show it at a level below half of the observed rate (upper limit). Blocks: both builds, only {self.ga}, only {self.gb}."))
        for name, fa, fb in (("Both builds differ from the data", True, True), (f"Only {self.ga} differs from the data", True, False), (f"Only {self.gb} differs from the data", False, True)):
            sel = [r for r in rows if self.differs_data(r, "A") == fa and self.differs_data(r, "B") == fb]
            sel.sort(key=lambda r: -r["data"])
            rws = [[f"{r['energy']:.1f}", nuc(LC.split_id(r["id"])[0]), self.rate_pm(r, "D"), self.rate_pm(r, "A"), self.rate_pm(r, "B"),
                    self.sim_over_obs(r, "A", False), self.sim_over_obs(r, "B", False)] for r in sel]
            s.append(f"\\paragraph{{{name} ({len(sel)}).}}\n" + (self._long(rws, ["E (keV)", "Transition", "Data (s$^{-1}$)", "A (s$^{-1}$)", "B (s$^{-1}$)", "A/Data", "B/Data"], "lp{2.6cm}rrrrr") if rws else "None.\n"))
        return "\n".join(s)

    def _appendix_lines(self, rows):
        s = ["\\appendix\\section{All lines (layout of Tables 5 and 6)}\n",
             guide("All lines of Tables 5 and 6 of the paper in the layout of the paper: line energy, published rate, rate fitted in the data, simulated rates of A and B, ratios Sim/Obs and A/B, transition and lab energy. "
                   "$<$ marks an upper limit; $\\cdots$ means the ratio is not defined; $\\dagger$ marks lines where A and B differ by more than 3$\\sigma$ and 10\\%."),
             "\\newgeometry{margin=0.6in}\\pdfpagewidth=11in\\pdfpageheight=8.5in\\paperwidth=11in\\paperheight=8.5in \\setlength{\\textwidth}{9.8in}\\setlength{\\textheight}{7.1in}\\setlength{\\oddsidemargin}{-0.4in}\\setlength{\\evensidemargin}{-0.4in}\\setlength{\\topmargin}{-0.4in}\\setlength{\\headheight}{0pt}\\setlength{\\headsep}{0pt}\\setlength{\\footskip}{0.35in}\\setlength{\\linewidth}{9.8in}\\setlength{\\columnwidth}{9.8in}\\hsize=9.8in\n"]
        for tb, name in (("5", "Identified lines; blends of identified lines"), ("6", "Unidentified or suspected lines")):
            sel = sorted([r for r in rows if r["table"] == tb], key=lambda r: r["energy"])
            s.append(f"\\begin{{center}}\\textbf{{TABLE {tb}}}\\\\ {name}\\end{{center}}\n{{\\scriptsize\\setlength{{\\tabcolsep}}{{3pt}}\\begin{{longtable}}{{rrrrrrrrp{{4.0cm}}rp{{3.3cm}}}}\\toprule\n"
                     "Energy & Paper & Fit (data) & A & B & A/Obs & B/Obs & A/B & Transition & Lab & Comments\\\\ (keV) & (s$^{-1}$) & (s$^{-1}$) & (s$^{-1}$) & (s$^{-1}$) & & & & (levels) & (keV) & \\\\\\midrule\\endfirsthead\n"
                     "\\toprule Energy & Paper & Fit (data) & A & B & A/Obs & B/Obs & A/B & Transition & Lab & Comments\\\\\\midrule\\endhead\\bottomrule\\endfoot\n")
            for r in sel:
                tr, lab, com = LC.split_id(r["id"])
                dag = "$^\\dagger$" if self.differs_ab(r) else ""
                s.append(" & ".join([f"{r['energy']:.1f}", sci(r["paper"], 3), self.rate_paren(r, "D"), self.rate_paren(r, "A"), self.rate_paren(r, "B"),
                                     self.sim_over_obs(r, "A"), self.sim_over_obs(r, "B"), self.a_over_b(r) + dag, nuc(tr), esc(lab), esc(com)]) + "\\\\")
            s.append("\\end{longtable}}\n\\clearpage\n")
        s.append("\\clearpage\\pdfpagewidth=8.5in\\pdfpageheight=11in\\paperwidth=8.5in\\paperheight=11in \\restoregeometry\\setlength{\\linewidth}{\\textwidth}\\setlength{\\columnwidth}{\\textwidth}\\hsize=\\textwidth\n")
        return "\n".join(s)

    def _appendix_groups(self, rows):
        has_c = bool(rows) and "cA" in rows[0]
        # group_summary reads r["id"], r["energy"], r["A"], r["B"], r["eA"], r["eB"], optional compA/compB
        for r in rows:
            r["compA"], r["compB"] = r.get("cA"), r.get("cB")
        summ, use = LC.group_summary(rows, has_c)
        s = ["\\section{Where the line differences come from (group patterns)}\n",
             guide("The line tables list each line separately. Here the lines are grouped by properties known in advance (the origin flag of the paper, the energy, whether the line comes mostly from the "
                   "activation component, the element of the transition and the decay parent) to see whether the differences between A and B follow a pattern. For every group the median of B/A is given with a "
                   "95\\% bootstrap interval. Only lines detected in both builds (net rate more than 3$\\sigma$ above zero) are used. A group is marked \\textbf{differs} if it has at least 3 lines and the interval "
                   "excludes 1; groups with fewer lines are not interpreted. The analysis is descriptive: it shows where differences cluster, it does not prove a cause."),
             f"\\begin{{center}}\\includegraphics[width=0.92\\linewidth]{{{self.fig_forest(summ, 'forest.png')}}}\\end{{center}}\n"]
        rws = [[nuc(x["group"]), x["n"], sci(x["median_B_over_A"], 3), (f"{x['ci_low']:.2f} -- {x['ci_high']:.2f}" if x["n"] >= 3 else "--"),
                f"{x['n_B_higher']} of {x['n']}", "\\textbf{differs}" if x["differs"] else ("n $<$ 3" if x["n"] < 3 else "no"), esc(x["lines"])] for x in summ]
        s.append(self._long(rws, ["Group", "Lines", "Median B/A", "95\\% interval", "B higher in", "Result", "Line energies (keV)"], "p{3.9cm}rrrrlp{3.2cm}"))
        return "\n".join(s)

    def compile(self, tex):
        exe = shutil.which("pdflatex")
        if not exe:
            print("pdflatex not found. The .tex file and figs/ folder are ready: compile them elsewhere "
                  "(Overleaf, or `module load texlive` then `pdflatex`).")
            return tex
        for _ in range(2):
            p = subprocess.run([exe, "-interaction=nonstopmode", "-halt-on-error", tex.name], cwd=tex.parent,
                               capture_output=True, text=True)
            if p.returncode:
                print(p.stdout[-3000:]); raise SystemExit("pdflatex failed, see above (log: " + str(tex.with_suffix(".log")) + ")")
        pdf = tex.with_suffix(".pdf"); print(f"Wrote {pdf}")
        return pdf


def main():
    ap = argparse.ArgumentParser(description="LaTeX-layout validation report from TGRS_output.root of two builds")
    ap.add_argument("root_a"); ap.add_argument("root_b")
    ap.add_argument("--tgrs", required=True, help="TGRS.dat (kept for the record; data are read from the TGRS histogram)")
    ap.add_argument("--ref-label", default="Reference (A)"); ap.add_argument("--test-label", default="Test (B)")
    ap.add_argument("--lines", default=str(HERE / "tgrs_lines.csv"))
    ap.add_argument("--out", default="latex_report")
    ap.add_argument("--comp-labels", nargs="*", default=["Cosmic Photons", "Cosmic Protons Prompt", "Cosmic Protons Activation"])
    ap.add_argument("--activation-component", type=int, default=2)
    ap.add_argument("--sim-merge", type=int, default=1, help="merge this many simulation bins above 245 keV (default 1)")
    ap.add_argument("--bin-tol", type=float, default=0.10); ap.add_argument("--min-fraction", type=float, default=0.95)
    ap.add_argument("--band-tol", type=float, default=0.05)
    ap.add_argument("--fwhm0", type=float, default=0.9); ap.add_argument("--fwhm1", type=float, default=0.0013)
    ap.add_argument("--no-pdf", action="store_true", help="only write the .tex file and the figures")
    a = ap.parse_args()
    rep = TGRS_LatexReporter(a.root_a, a.root_b, a.tgrs, a.out, a.ref_label, a.test_label, lines_csv=a.lines,
                             comp_labels=a.comp_labels, activation_component=a.activation_component, sim_merge=a.sim_merge,
                             bin_tol=a.bin_tol, min_fraction=a.min_fraction, band_tol=a.band_tol, fwhm0=a.fwhm0, fwhm1=a.fwhm1)
    rep.build(compile_pdf=not a.no_pdf)


if __name__ == "__main__":
    main()
