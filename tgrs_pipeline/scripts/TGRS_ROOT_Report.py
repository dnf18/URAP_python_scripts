#!/usr/bin/env python3
"""
TGRS_ROOT_Report.py

Builds the validation PDF for two MEGAlib builds from the ROOT files written by the
MEGAlib TGRS program (patched to save its histograms), instead of parsing the *.sim.gz
files in Python. The TGRS program applies the detector energy resolution and noise and
fills fine log bins, which is the figure Andreas asked for.

Pipeline per build X (A = reference, B = test), run in that build's own environment:
    TGRS -t TGRS.dat -g TGRS.geo.setup \
         -b "Cosmic Photons"            CosmicPhotons.p1.sim.gz \
         -b "Cosmic Protons Prompt"     CosmicProtonsPrompt.p1.sim.gz \
         -b "Cosmic Protons Activation" CosmicProtons.Decays.p1.sim.gz
    -> TGRS_output.root   (histograms TGRS, Sim, Component_0 ... Component_N)

This script then
  * runs Lindsey's Comparator class (ratio mode) for A vs measurements and B vs measurements
    -> her spectrum plots and ratio plots (sections 1-4 of the report),
  * compares A and B directly (our statistics: Poisson pulls, chi2, KS, per-band chi2),
    makes the A-B residual plot and an overlay of both builds (section 5, decides PASS/FAIL),
  * runs her Comparator in chi2 mode for A vs B as a cross-check of the numbers,
  * writes everything in the layout TGRS_Reporter.py expects and builds the PDF.

Example:
    python3 TGRS_ROOT_Report.py a/TGRS_output.root b/TGRS_output.root --tgrs inputs_dc6/TGRS.dat \
        --ref-label "MEGAlib 10.2 (Geant4 10.2.3)" --test-label "MEGAlib dc6 (Geant4 11.4.3)" \
        --comparator-dir /path/to/folder/with/Comparator.py \
        --out analysis/results_root --reports analysis/reports_root

Needs: numpy, scipy, matplotlib, uproot, reportlab (module load anaconda3; pip install --user uproot reportlab).
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from TGRS_Comparator import (load_tgrs, rebin_measurement, compare_rates, compare_ratio,   # noqa: E402
                             plot_residuals, SIGMA_THRESHOLD)


# ------------------------------------------------------------------ reading ROOT

def read_hist(path, name):
    """(contents, edges, errors) of one histogram of a ROOT file."""
    import uproot
    with uproot.open(path) as f:
        keys = [k.split(";")[0] for k in f.keys()]
        if name not in keys:
            raise KeyError(f"'{name}' not in {path}. Available: {keys}")
        h = f[name]
        return (np.asarray(h.values(), dtype=float),
                np.asarray(h.axis().edges(), dtype=float),
                np.asarray(h.errors(), dtype=float))


def check_poisson_errors(label, values, errors):
    """The stored errors must be Poisson on the counts, scaled like the contents.
    If TGRS.cxx did not call Sumw2() before filling, ROOT reports sqrt(content) of the
    already scaled rates, which is wrong by a large factor and makes chi2 meaningless."""
    nz = values > 0
    if not np.any(nz):
        print(f"[check] {label}: histogram is empty")
        return False
    r = errors[nz] / np.sqrt(values[nz])
    med = float(np.median(r))
    ok = med < 0.5
    print(f"[check] {label}: median error / sqrt(rate) = {med:.3g} -> "
          + ("looks like Poisson on counts" if ok else
             "WARNING: errors look like sqrt(rate). Call Sumw2() on the histograms in TGRS.cxx"))
    return ok


# ------------------------------------------------------------------ binning helpers

def remap_density(values, errors, src_edges, dst_edges):
    """Map a per-keV density histogram onto other bins (width-averaged). Bins not fully
    covered by the source become NaN."""
    n = len(dst_edges) - 1
    out = np.full(n, np.nan)
    oerr = np.full(n, np.nan)
    for i in range(n):
        lo, hi = dst_edges[i], dst_edges[i + 1]
        if lo < src_edges[0] - 1e-9 or hi > src_edges[-1] + 1e-9:
            continue
        ov = np.clip(np.minimum(src_edges[1:], hi) - np.maximum(src_edges[:-1], lo), 0, None)
        w = ov / (hi - lo)
        out[i] = np.sum(values * w)
        oerr[i] = np.sqrt(np.sum((errors * w) ** 2))
    return out, oerr


# ------------------------------------------------------------------ plotting

def plot_overlay_ab(path, edges_a, a, ea, edges_b, b, eb, d_edges, d, de, label_a, label_b, emin):
    """Measurements and both simulated spectra at full resolution, paper style."""
    fig, ax = plt.subplots(figsize=(12, 7))
    for edges, v, e, color, label in ((edges_a, a, ea, "tab:blue", f"{label_a} (combined)"),
                                      (edges_b, b, eb, "tab:red", f"{label_b} (combined)")):
        ax.stairs(np.where(v > 0, v, np.nan), edges, color=color, linewidth=0.8, label=label)
        mid = 0.5 * (edges[:-1] + edges[1:])
        ax.fill_between(mid, np.maximum(v - e, 1e-12), v + e, step="mid", color=color,
                        alpha=0.2, linewidth=0)
    ax.stairs(np.where(d > 0, d, np.nan), d_edges, color="black", linewidth=1.0, label="Measurements (x35 cm2)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(emin, d_edges[-1])
    sel = (d_edges[:-1] >= emin) & (d > 0)
    if np.any(sel):
        ax.set_ylim(0.2 * float(np.min(d[sel])), 3.0 * float(np.max(d[sel])))
    ax.set_xlabel("Energy in keV", loc="right")
    ax.set_ylabel("Hits in cts/keV/s")
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.legend(loc="upper right", frameon=True, edgecolor="black", fancybox=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"Saved {path}")


# ------------------------------------------------------------------ main logic

def import_lindsey(comparator_dir):
    if comparator_dir:
        sys.path.insert(0, str(Path(comparator_dir).resolve()))
    try:
        from Comparator import Comparator
    except ImportError as e:
        raise SystemExit("Could not import Lindsey's Comparator.py ({}). Pass --comparator-dir "
                         "with the folder that contains it.".format(e))
    return Comparator


def make_paper_comparator(Base, emin, merge):
    """Lindsey's Comparator with the same plots, but framed like the reference figure:
    x starts at emin (default 50 keV, the simulation has nothing below the detector threshold),
    sensible y limits, and the ratio plot on a log axis with off-scale values clipped (a single
    near-empty bin at the upper end made the linear ratio plot useless).
    Statistics are not touched."""

    class PaperComparator(Base):

        def plot_overlay(self, edges, ref_n, test_orig, components):
            fig, ax = plt.subplots(figsize=(10, 7))
            for i, (label, vals, e) in enumerate(components):
                ax.stairs(np.where(vals > 0, vals, np.nan), e, color=self.component_color(label, i),
                          linewidth=0.8, label=label)
            test_n, test_e = test_orig
            ax.stairs(np.where(test_n > 0, test_n, np.nan), test_e, color="#e41a1c",
                      linewidth=0.9, label=self.test_label)
            ax.stairs(np.where(ref_n > 0, ref_n, np.nan), edges, color="black",
                      linewidth=1.2, label=self.ref_label)
            handles, labels = ax.get_legend_handles_labels()
            order = [len(labels) - 1, len(labels) - 2] + list(range(len(components)))
            ax.legend([handles[i] for i in order], [labels[i] for i in order],
                      loc="upper right", frameon=True, edgecolor="black", fancybox=False)
            ax.set_xlabel("Energy in keV", loc="right")
            ax.set_ylabel("Hits in cts/keV/s" if self.normalization == "none" and self.density
                          else f"Hits ({self.normalization})")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(emin, max(edges[-1], test_e[-1]))
            sel = np.isfinite(ref_n) & (edges[:-1] >= emin) & (ref_n > 0)
            top = float(np.max(ref_n[sel])) if np.any(sel) else 1.0
            ax.set_ylim(1e-5, top * 3.0)
            ax.tick_params(which="both", direction="in", top=True, right=True)
            fig.tight_layout()
            output_path = self.hist_dir / "energy_overlay.png"
            fig.savefig(output_path, dpi=200)
            plt.close(fig)
            print(f"[Comparator] Spectrum plot saved -> {output_path}")
            return output_path

        def plot_ratio(self, edges, panel):
            if self.mode != "ratio":
                return Base.plot_ratio(self, edges, panel)
            # merge neighbouring bins for display only (statistics use the unmerged ratio)
            idx = np.arange(len(panel)) // max(merge, 1)
            ng = int(idx.max()) + 1
            fin = np.isfinite(panel)
            cnt = np.bincount(idx, weights=fin.astype(float), minlength=ng)
            tot = np.bincount(idx, weights=np.where(fin, panel, 0.0), minlength=ng)
            merged = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
            e_m = edges[::max(merge, 1)]
            if e_m[-1] != edges[-1]:
                e_m = np.append(e_m, edges[-1])
            lo, hi = 0.2, 5.0
            fig, ax = plt.subplots(figsize=(10, 3.5))
            ax.stairs(np.clip(merged, lo, hi), e_m, color="tab:blue")
            ax.axhspan(1 - self.tolerance, 1 + self.tolerance, color="gray", alpha=0.15,
                       label=f"+-{self.tolerance:.0%}")
            ax.axhline(1.0, color="black", linewidth=0.6)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(emin, edges[-1])
            ax.set_ylim(lo, hi)
            ax.set_yticks([0.2, 0.5, 1, 2, 5])
            ax.yaxis.set_minor_formatter(NullFormatter())
            ax.set_yticklabels(["0.2", "0.5", "1", "2", "5"])
            ax.set_ylabel(f"{self.test_label} / {self.ref_label}", fontsize=8)
            ax.set_xlabel("Energy in keV", loc="right")
            ax.legend(loc="upper right")
            ax.tick_params(which="both", direction="in", top=True, right=True)
            fig.tight_layout()
            output_path = self.hist_dir / "energy_ratio.png"
            fig.savefig(output_path, dpi=200)
            plt.close(fig)
            print(f"[Comparator] Ratio plot saved -> {output_path}")
            return output_path

    return PaperComparator


def clean(stats):
    """Drop plot paths that would only clutter the PDF tables."""
    return {k: v for k, v in stats.items() if k not in ("overlay_plot", "ratio_plot")}


TOL_BANDS = [(0.0, 100.0), (100.0, 1000.0), (1000.0, float("inf"))]


def tolerance_test(a, ea, b, eb, edges, bin_tol, min_fraction, band_tol, nsig=3.0):
    """Build-to-build test that stays meaningful at high statistics.

    With ~1e5 s of simulation every bin has a tiny statistical error, so a chi2 p-value fails for any
    real difference of a few percent. This test instead asks whether the two spectra agree within a
    stated physical tolerance, allowing for the statistical error:
      * per bin (merged fine bins): |A-B| <= bin_tol*max(A,B) + nsig*sigma  -> fraction of such bins
        must be >= min_fraction
      * per energy band (<100, 100-1000, >1000 keV): the integrated rate of A and B agree within
        band_tol (plus nsig*sigma of the integrals)
    Individual lines are judged separately (line comparison pages, information only).
    """
    sigma = np.sqrt(ea ** 2 + eb ** 2)
    used = (sigma > 0) & ~((a == 0) & (b == 0))
    ok = np.abs(a - b) <= bin_tol * np.maximum(np.abs(a), np.abs(b)) + nsig * sigma
    frac = float(np.mean(ok[used])) if np.any(used) else float("nan")
    rows = {
        "Tolerance test: per-bin tolerance": f"{bin_tol:.0%} (+ {nsig:g} sigma)",
        "Tolerance test: fraction of bins within tolerance": round(frac, 4),
        "Tolerance test: required fraction": min_fraction,
    }
    width = np.diff(edges)
    mid = 0.5 * (edges[:-1] + edges[1:])
    bands_ok = True
    for lo, hi in TOL_BANDS:
        sel = (mid >= lo) & (mid < hi)
        label = f"< {hi:g} keV" if lo == 0 else (f"{lo:g} - {hi:g} keV" if np.isfinite(hi) else f"> {lo:g} keV")
        if not np.any(sel):
            continue
        ia, ib = float(np.sum(a[sel] * width[sel])), float(np.sum(b[sel] * width[sel]))
        sia = float(np.sqrt(np.sum((ea[sel] * width[sel]) ** 2)))
        sib = float(np.sqrt(np.sum((eb[sel] * width[sel]) ** 2)))
        good = abs(ia - ib) <= band_tol * max(abs(ia), abs(ib)) + nsig * float(np.hypot(sia, sib))
        bands_ok &= good
        rows[f"Tolerance test, band {label}: integrated rate A/B"] = round(ia / ib, 4) if ib else float("nan")
        rows[f"Tolerance test, band {label}: within {band_tol:.0%}"] = bool(good)
    passed = bool(frac >= min_fraction and bands_ok)
    rows["Tolerance test: pass"] = passed
    return rows, passed


def main():
    ap = argparse.ArgumentParser(description="PDF validation report from TGRS_output.root of two builds")
    ap.add_argument("root_a", help="TGRS_output.root of the reference build (A)")
    ap.add_argument("root_b", help="TGRS_output.root of the test build (B)")
    ap.add_argument("--tgrs", required=True, help="TGRS.dat (measurements, used for the statistics)")
    ap.add_argument("--ref-label", default="Reference (A)")
    ap.add_argument("--test-label", default="Test (B)")
    ap.add_argument("--tgrs-hist", default="TGRS", help="measurement histogram name in the ROOT files")
    ap.add_argument("--sim-hist", default="Sim", help="combined simulation histogram name")
    ap.add_argument("--comparator-dir", default=None, help="folder containing Lindsey's Comparator.py")
    ap.add_argument("--rebin", type=int, default=10,
                    help="merge this many TGRS.dat bins for the simulation-vs-measurement statistics (default 10)")
    ap.add_argument("--ab-rebin", type=int, default=5,
                    help="merge this many fine simulation bins for the A-vs-B statistics and residual plot "
                         "(default 5; fine bins hold few counts at high energy)")
    ap.add_argument("--plot-emin", type=float, default=50.0)
    ap.add_argument("--bin-tol", type=float, default=0.10,
                    help="tolerance test: a merged bin agrees if |A-B| <= bin-tol*max(A,B) + 3 sigma (default 0.10)")
    ap.add_argument("--min-fraction", type=float, default=0.95,
                    help="tolerance test: required fraction of agreeing bins (default 0.95)")
    ap.add_argument("--band-tol", type=float, default=0.05,
                    help="tolerance test: integrated rate in each energy band must agree within this fraction "
                         "(plus 3 sigma) (default 0.05)")
    ap.add_argument("--chi2-decides", action="store_true",
                    help="decide PASS/FAIL by the old chi2 p-value instead of the tolerance test")
    ap.add_argument("--out", default="results_root")
    ap.add_argument("--reports", default="./reports_root")
    ap.add_argument("--no-report", action="store_true")
    args = ap.parse_args()

    Comparator = make_paper_comparator(import_lindsey(args.comparator_dir), args.plot_emin, args.rebin)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # --- read histograms
    a, ea_edges, ea = None, None, None
    sim_a, edges_a, err_a = read_hist(args.root_a, args.sim_hist)
    sim_b, edges_b, err_b = read_hist(args.root_b, args.sim_hist)
    if len(edges_a) != len(edges_b) or not np.allclose(edges_a, edges_b):
        raise SystemExit("The Sim histograms of A and B have different binning; "
                         "run TGRS with the same settings for both builds.")
    ok_a = check_poisson_errors(f"{args.ref_label} Sim", sim_a, err_a)
    ok_b = check_poisson_errors(f"{args.test_label} Sim", sim_b, err_b)

    d_edges, d_rate, d_err = load_tgrs(args.tgrs)

    # --- sections 1-4: Lindsey's comparator, simulation vs measurements, one run per build
    results = {}
    for tag, root, label in (("A", args.root_a, args.ref_label), ("B", args.root_b, args.test_label)):
        sub = out / f"lindsey_{tag}"
        cmp_ = Comparator(root, root, sub / "comparison_results.json", sub,
                          ref_hist=args.tgrs_hist, test_hist=args.sim_hist, mode="ratio",
                          normalization="none", ref_label="Measurements", test_label=label,
                          density=True)
        results[tag] = cmp_.compare()
        shutil.copy(sub / "energy_overlay.png", out / ("plot1_sim_a.png" if tag == "A" else "plot2_sim_b.png"))
        shutil.copy(sub / "energy_ratio.png", out / ("plot3_ratio_a.png" if tag == "A" else "plot4_ratio_b.png"))

    # --- our statistics of simulation vs measurements (on the TGRS.dat bins, covered range only)
    stats_vs_meas = {}
    for tag, sim, err, label in (("A", sim_a, err_a, args.ref_label), ("B", sim_b, err_b, args.test_label)):
        m, em = remap_density(sim, err, edges_a, d_edges)
        good = np.where(np.isfinite(m))[0]
        lo, hi = good[0], good[-1] + 1
        e_cov = d_edges[lo:hi + 1]
        e_r, m_r, em_r = rebin_measurement(e_cov, m[lo:hi], em[lo:hi], args.rebin)
        _, d_r, dd_r = rebin_measurement(e_cov, d_rate[lo:hi], d_err[lo:hi], args.rebin)
        s_rates, _ = compare_rates(m_r, em_r, d_r, dd_r, e_r, label, "Measurements")
        stats_vs_meas[tag] = s_rates
    plots = {"A": out / "plot1_sim_a.png", "B": out / "plot2_sim_b.png"}

    def save(name, stats, plot):
        stats = clean(stats)
        stats["overlay_plot"] = str(plot)
        with open(out / name, "w") as f:
            json.dump(stats, f, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))
        print(f"Saved {name}")

    save("reference_stats.json", stats_vs_meas["A"], plots["A"])
    save("test_stats.json", stats_vs_meas["B"], plots["B"])
    save("reference_ratio_stats.json", results["A"], out / "plot3_ratio_a.png")
    save("test_ratio_stats.json", results["B"], out / "plot4_ratio_b.png")

    # --- section 5: A versus B (decides PASS/FAIL), residuals, overlay
    e_m, a_m, ea_m = rebin_measurement(edges_a, sim_a, err_a, args.ab_rebin)
    _, b_m, eb_m = rebin_measurement(edges_a, sim_b, err_b, args.ab_rebin)
    s5, pulls_ab = compare_rates(a_m, ea_m, b_m, eb_m, e_m, args.ref_label, args.test_label)

    cross = Comparator(args.root_a, args.root_b, out / "lindsey_AB" / "comparison_results.json",
                       out / "lindsey_AB", ref_hist=args.sim_hist, test_hist=args.sim_hist, mode="chi2",
                       normalization="none", ref_label=args.ref_label, test_label=args.test_label,
                       density=True).compare()
    s5["Cross-check, Comparator.py chi2 mode: Chi2/NDF (fine bins)"] = cross.get("Chi2/NDF")
    s5["Cross-check, Comparator.py chi2 mode: p-value (fine bins)"] = cross.get("Chi2 p-value")
    s5["Statistics bins"] = f"{len(a_m)} (fine bins merged x{args.ab_rebin})"
    # PASS/FAIL: tolerance test by default; the chi2 p-value is kept as information
    tol_rows, tol_pass = tolerance_test(a_m, ea_m, b_m, eb_m, e_m, args.bin_tol, args.min_fraction, args.band_tol)
    chi2_pass = bool(s5.get("pass", False))
    s5.pop("pass", None)
    s5["Chi2 test (information only): pass"] = chi2_pass
    s5.update(tol_rows)
    s5["pass"] = chi2_pass if args.chi2_decides else tol_pass

    p5 = out / "plot5_all.png"
    plot_overlay_ab(p5, edges_a, sim_a, err_a, edges_b, sim_b, err_b, d_edges, d_rate, d_err,
                    args.ref_label, args.test_label, args.plot_emin)
    save("combined_stats.json", s5, p5)
    plot_residuals(out / "plot6_residual_ab.png", e_m, a_m, ea_m, b_m, eb_m, pulls_ab,
                   args.ref_label, args.test_label)

    config = {
        "Reference (A)": args.ref_label,
        "Test (B)": args.test_label,
        "Source of the spectra": "ROOT files written by the MEGAlib TGRS program (with Lindsey's output patch)",
        "Reference file": str(args.root_a),
        "Test file": str(args.root_b),
        "TGRS data file": str(args.tgrs),
        "Energy resolution applied": "Yes - detector resolution and noise as in TGRS.cxx",
        "Plot resolution": f"{len(sim_a)} simulation bins, {len(d_rate)} measurement bins",
        "Statistics bins vs measurements": f"TGRS.dat bins merged x{args.rebin}",
        "Statistics bins A vs B": f"{len(a_m)} (fine bins merged x{args.ab_rebin})",
        "Errors checked as Poisson": f"A: {ok_a}, B: {ok_b}",
        "Overall PASS/FAIL decided by": ("A versus B, chi2 p-value" if args.chi2_decides else
            f"A versus B tolerance test: >= {args.min_fraction:.0%} of bins within {args.bin_tol:.0%} (or 3 sigma) and each band (<100, 100-1000, >1000 keV) within {args.band_tol:.0%} (or 3 sigma)"),
    }
    with open(out / "run_config.json", "w") as f:
        json.dump(config, f, indent=2)

    if not args.no_report:
        from TGRS_Reporter import TGRS_Reporter
        TGRS_Reporter(config_json=str(out / "run_config.json"),
                      output_dir=args.reports).generate_pdf(stats_dir=out)


if __name__ == "__main__":
    main()
