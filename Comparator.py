"""
Comparator class (generalized)

Two kinds of comparison:

  1. Data vs. model  (mode="ratio")  -- the TGRS benchmark
     Reference = measured TGRS spectrum  ("TGRS" histogram, from TGRS.dat)
     Test      = summed simulation        ("Sim" histogram = photons + prompt protons + activation)
     Both live in the same TGRS_output.root and are already in cts/keV/s, so use normalization "none".
     The two histograms have different binning (TGRS.dat bins vs. 1500 fine log bins); the test
     spectrum is rebinned onto the reference bins automatically.
     Metrics: per-bin ratio Sim/Data, median ratio, typical deviation, fraction of bins within
     +-tolerance, worst energy range.

  2. Simulation vs. simulation  (mode="chi2")  -- e.g. Geant4 10.2 vs 11.4
     Per-bin pulls, chi2/ndf, chi2 p-value, bins beyond N sigma, KS distance.

Inputs can be JSON ({"bins": [...], "edges": [...]}, optional "errors", "time") or ROOT files.
For ROOT files, ref_hist / test_hist choose the histogram (both may come from the same file).

Normalization (applied before comparing):
  "none"   as stored (TGRS histograms are already cts/keV/s)
  "area"   unit integral -> shape-only comparison
  "time"   divide by observation time (needs ref_time / test_time)
  "width"  divide by bin width
  combinable with "+", e.g. "time+width"

Command line examples:
    # TGRS benchmark: measured data vs summed simulation, both from one file
    python3 Comparator.py TGRS_output.root --ref-hist TGRS --test-hist Sim \
        --mode ratio --norm none --ref-label "TGRS data" --test-label "Simulation (sum)"

    # Geant4 10.2 vs 11.4: same histogram from two files
    python3 Comparator.py g4_10.2.root g4_11.4.root --ref-hist Sim --test-hist Sim \
        --mode ratio --norm none --ref-label "Geant4 10.2" --test-label "Geant4 11.4"
"""

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")  # no display needed (works on Savio)
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_dist


class Comparator:

    def __init__(self, ref_json, test_json, output_json, histogram_output_dir,
                 sigma_threshold=3.0, hist_name="Sim", normalization="area",
                 ref_time=None, test_time=None, p_threshold=0.05,
                 ref_label="Reference", test_label="Test",
                 log_x=True, log_y=True,
                 ref_hist=None, test_hist=None, mode="chi2",
                 tolerance=0.25, pass_fraction=0.8, density=True,
                 component_prefix="Component_", dpi=400):
        # ref_json / test_json keep their old names for compatibility,
        # but they may point to .json or .root files (and may be the same file).
        self.ref_path = Path(ref_json)
        self.test_path = Path(test_json)
        self.output_json = Path(output_json)
        self.hist_dir = Path(histogram_output_dir)
        self.hist_dir.mkdir(parents=True, exist_ok=True)

        self.ref_hist = ref_hist or hist_name
        self.test_hist = test_hist or hist_name
        self.mode = mode
        self.sigma_threshold = sigma_threshold
        self.p_threshold = p_threshold
        self.tolerance = tolerance          # ratio mode: allowed fractional deviation per bin
        self.pass_fraction = pass_fraction  # ratio mode: fraction of bins that must be within tolerance
        self.density = density              # True if contents are per-keV densities (TGRS output)
        self.normalization = normalization.lower()
        self.ref_time = ref_time
        self.test_time = test_time
        self.ref_label = ref_label
        self.test_label = test_label
        self.log_x = log_x
        self.log_y = log_y
        self.component_prefix = component_prefix
        self.dpi = dpi  # PNG resolution; a vector PDF copy is saved too

    # ------------------------------------------------------------------ loading

    def load_components(self, path):
        """Individual simulation components (photons, prompt protons, activation, ...)
        saved by TGRS.cxx as Component_0, Component_1, ... with the -b name as title.
        Returns a list of (label, contents, edges); empty if none are stored."""
        path = Path(path)
        if path.suffix != ".root":
            return []
        import uproot
        comps = []
        with uproot.open(path) as f:
            names = sorted({k.split(";")[0] for k in f.keys() if k.startswith(self.component_prefix)})
            for name in names:
                h = f[name]
                title = h.member("fTitle") or name
                comps.append((title, np.asarray(h.values(), dtype=float),
                              np.asarray(h.axis().edges(), dtype=float)))
        return comps

    def load_hist(self, path, name):
        """Returns (contents, edges, errors, time). errors/time may be None."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Histogram file not found: {path}")

        if path.suffix == ".root":
            import uproot
            with uproot.open(path) as f:
                keys = [k.split(";")[0] for k in f.keys()]
                if name not in keys:
                    raise KeyError(f"'{name}' not in {path}. Available: {keys}")
                h = f[name]
                contents = np.asarray(h.values(), dtype=float)
                edges = np.asarray(h.axis().edges(), dtype=float)
                errors = np.asarray(h.errors(), dtype=float)
            return contents, edges, errors, None

        with open(path) as f:
            data = json.load(f)
        if "bins" not in data or "edges" not in data:
            raise KeyError(f"{path} needs keys 'bins' and 'edges'. Found: {list(data.keys())}")
        errors = np.asarray(data["errors"], dtype=float) if "errors" in data else None
        return (np.asarray(data["bins"], dtype=float),
                np.asarray(data["edges"], dtype=float),
                errors,
                data.get("time"))

    # ------------------------------------------------------------ normalization

    def normalize(self, contents, errors, edges, obs_time):
        if errors is None:
            errors = np.sqrt(np.abs(contents))  # Poisson fallback on raw counts
        contents = contents.copy()
        errors = errors.copy()

        for mode in [m.strip() for m in self.normalization.split("+") if m.strip()]:
            if mode == "none":
                continue
            elif mode == "width":
                widths = np.diff(edges)
                contents /= widths
                errors /= widths
            elif mode == "time":
                if not obs_time:
                    raise ValueError("Normalization 'time' needs ref_time and test_time (seconds)")
                contents /= obs_time
                errors /= obs_time
            elif mode == "area":
                total = np.sum(contents * (np.diff(edges) if self.density else 1.0))
                if total == 0:
                    raise ValueError("Cannot area-normalize a histogram with zero total content")
                contents /= total
                errors /= total
            else:
                raise ValueError(f"Unknown normalization mode '{mode}'")
        return contents, errors

    # ---------------------------------------------------------------- rebinning

    def rebin(self, values, errors, src_edges, dst_edges):
        """Maps a histogram onto new bin edges. Bins not fully covered become NaN.
        Densities (per keV) are width-averaged; counts are split by overlap fraction."""
        n = len(dst_edges) - 1
        out = np.full(n, np.nan)
        out_err = np.full(n, np.nan)
        src_w = np.diff(src_edges)
        for i in range(n):
            lo, hi = dst_edges[i], dst_edges[i + 1]
            if lo < src_edges[0] - 1e-9 or hi > src_edges[-1] + 1e-9:
                continue
            overlap = np.clip(np.minimum(src_edges[1:], hi) - np.maximum(src_edges[:-1], lo), 0, None)
            w = overlap / (hi - lo) if self.density else overlap / src_w
            out[i] = np.sum(values * w)
            out_err[i] = np.sqrt(np.sum((errors * w) ** 2))
        return out, out_err

    # --------------------------------------------------------------- statistics

    def moments(self, contents, edges):
        weights = contents * (np.diff(edges) if self.density else 1.0)
        ok = np.isfinite(weights) & (weights > 0)
        if not np.any(ok):
            return 0.0, 0.0
        centers = 0.5 * (edges[:-1] + edges[1:])[ok]
        mean = np.average(centers, weights=weights[ok])
        sigma = np.sqrt(np.average((centers - mean) ** 2, weights=weights[ok]))
        return float(mean), float(sigma)

    def compare(self):
        ref_c, ref_e, ref_err, ref_t = self.load_hist(self.ref_path, self.ref_hist)
        test_c, test_e, test_err, test_t = self.load_hist(self.test_path, self.test_hist)

        ref_n, ref_s = self.normalize(ref_c, ref_err, ref_e, self.ref_time or ref_t)
        test_n, test_s = self.normalize(test_c, test_err, test_e, self.test_time or test_t)
        test_orig = (test_n.copy(), test_e.copy())  # kept for plotting at full resolution

        rebinned = len(ref_e) != len(test_e) or not np.allclose(ref_e, test_e)
        if rebinned:
            test_n, test_s = self.rebin(test_n, test_s, test_e, ref_e)
        edges = ref_e

        ref_mu, ref_sig = self.moments(ref_n, edges)
        test_mu, test_sig = self.moments(test_n, edges)

        results = {
            "Mode": self.mode,
            "Reference histogram": f"{self.ref_path.name}:{self.ref_hist}",
            "Test histogram": f"{self.test_path.name}:{self.test_hist}",
            "Normalization": self.normalization,
            "Test rebinned onto reference bins": bool(rebinned),
            "Reference Mean (keV)": ref_mu,
            "Test Mean (keV)": test_mu,
            "Reference Standard Deviation (keV)": ref_sig,
            "Test Standard Deviation (keV)": test_sig,
        }

        if self.mode == "ratio":
            used = np.isfinite(test_n) & (ref_n > 0) & (test_n > 0)
            ratio = np.full_like(ref_n, np.nan)
            ratio[used] = test_n[used] / ref_n[used]
            logr = np.abs(np.log(ratio[used]))
            within = np.abs(ratio[used] - 1.0) <= self.tolerance
            frac_within = float(np.mean(within)) if np.any(used) else 0.0
            worst = int(np.nanargmax(np.abs(np.log(np.where(used, ratio, 1.0)))))
            panel = ratio
            results.update({
                "Bins compared": int(np.sum(used)),
                "Median ratio (test/ref)": float(np.median(ratio[used])),
                "Typical deviation (exp(median |ln ratio|) - 1)": float(np.exp(np.median(logr)) - 1),
                f"Fraction of bins within +-{self.tolerance:.0%}": frac_within,
                "Worst bin energy range (keV)": [float(edges[worst]), float(edges[worst + 1])],
                "Worst bin ratio": float(ratio[worst]),
                "pass": bool(frac_within >= self.pass_fraction),
            })
        else:
            sigma_tot = np.sqrt(ref_s ** 2 + test_s ** 2)
            used = np.isfinite(test_n) & (sigma_tot > 0)
            pulls = np.full_like(ref_n, np.nan)
            pulls[used] = (ref_n[used] - test_n[used]) / sigma_tot[used]
            chi2_val = float(np.sum(pulls[used] ** 2))
            ndf = int(np.sum(used)) - (1 if "area" in self.normalization else 0)
            chi2_p = float(chi2_dist.sf(chi2_val, ndf)) if ndf > 0 else float("nan")
            worst = int(np.nanargmax(np.abs(pulls)))
            r = np.where(used, ref_n, 0.0)
            t = np.where(used, test_n, 0.0)
            ks_d = float(np.max(np.abs(np.cumsum(r) / np.sum(r) - np.cumsum(t) / np.sum(t))))
            panel = pulls
            results.update({
                "Chi2": chi2_val,
                "NDF": ndf,
                "Chi2/NDF": chi2_val / ndf if ndf > 0 else float("nan"),
                "Chi2 p-value": chi2_p,
                f"Bins beyond {self.sigma_threshold:g} sigma": int(np.sum(np.abs(pulls[used]) > self.sigma_threshold)),
                "Worst bin energy range (keV)": [float(edges[worst]), float(edges[worst + 1])],
                "Worst bin pull": float(pulls[worst]),
                "KS distance": ks_d,
                "pass": bool(chi2_p > self.p_threshold),
            })

        components = self.load_components(self.test_path) if self.normalization == "none" else []
        results["Simulation components plotted"] = [c[0] for c in components]
        results["overlay_plot"] = str(self.plot_overlay(edges, ref_n, test_orig, components))
        results["ratio_plot"] = str(self.plot_ratio(edges, panel))

        with open(self.output_json, "w") as f:
            json.dump(results, f, indent=4)
        print(f"[Comparator] Results saved -> {self.output_json}")
        return results

    # ----------------------------------------------------------------- plotting

    @staticmethod
    def component_color(label, i):
        """Colors as in Zoglauer et al. (Cosima paper): prompt green, activation blue, photons yellow."""
        l = label.lower()
        if "activ" in l or "decay" in l or "delayed" in l:
            return "#0000ff"
        if "prompt" in l or "proton" in l:
            return "#33cc33"
        if "photon" in l:
            return "#e6c800"
        return ["#9467bd", "#8c564b", "#17becf", "#e377c2"][i % 4]

    def plot_overlay(self, edges, ref_n, test_orig, components):
        """Paper-style spectrum: measurements, combined simulation and each component."""
        fig, ax = plt.subplots(figsize=(10, 7))

        for i, (label, vals, e) in enumerate(components):
            ax.stairs(np.where(vals > 0, vals, np.nan), e, color=self.component_color(label, i),
                      linewidth=0.8, label=label)
        test_n, test_e = test_orig
        ax.stairs(np.where(test_n > 0, test_n, np.nan), test_e, color="#e41a1c",
                  linewidth=0.9, label=self.test_label)
        ax.stairs(np.where(ref_n > 0, ref_n, np.nan), edges, color="black",
                  linewidth=1.2, label=self.ref_label)

        # Legend in the paper's order: measurements, combined simulation, components
        handles, labels = ax.get_legend_handles_labels()
        order = [len(labels) - 1, len(labels) - 2] + list(range(len(components)))
        ax.legend([handles[i] for i in order], [labels[i] for i in order],
                  loc="upper right", frameon=True, edgecolor="black", fancybox=False)

        ax.set_xlabel("Energy in keV", loc="right")
        ax.set_ylabel("Hits in cts/keV/s" if self.normalization == "none" and self.density
                      else f"Hits ({self.normalization})")
        if self.log_x and edges[0] > 0:
            ax.set_xscale("log")
        if self.log_y:
            ax.set_yscale("log")
        ax.set_xlim(min(edges[0], test_e[0]), max(edges[-1], test_e[-1]))
        ax.tick_params(which="both", direction="in", top=True, right=True)

        fig.tight_layout()
        output_path = self.hist_dir / "energy_overlay.png"
        fig.savefig(output_path, dpi=self.dpi)
        fig.savefig(output_path.with_suffix(".pdf"))  # vector copy, sharp at any zoom
        plt.close(fig)
        print(f"[Comparator] Spectrum plot saved -> {output_path}")
        return output_path

    def plot_ratio(self, edges, panel):
        """Quantitative check: Sim/Data ratio (ratio mode) or pulls (chi2 mode)."""
        fig, ax = plt.subplots(figsize=(10, 3.5))
        finite = np.where(np.isfinite(panel))[0]  # skip uncovered bins at the ends
        if len(finite):
            a, b = finite[0], finite[-1] + 1
            ax.stairs(panel[a:b], edges[a:b + 1], color="tab:blue")
        if self.mode == "ratio":
            ax.axhspan(1 - self.tolerance, 1 + self.tolerance, color="gray", alpha=0.15,
                       label=f"+-{self.tolerance:.0%}")
            ax.axhline(1.0, color="black", linewidth=0.6)
            ax.set_ylabel(f"{self.test_label} / {self.ref_label}", fontsize=9)
            ax.legend(loc="upper right")
        else:
            for y in (self.sigma_threshold, -self.sigma_threshold):
                ax.axhline(y, color="gray", linestyle="--", linewidth=0.8)
            ax.axhline(0, color="black", linewidth=0.5)
            ax.set_ylabel("Pull (sigma)")
        ax.set_xlabel("Energy in keV", loc="right")
        if self.log_x and edges[0] > 0:
            ax.set_xscale("log")
        fig.tight_layout()
        output_path = self.hist_dir / "energy_ratio.png"
        fig.savefig(output_path, dpi=self.dpi)
        fig.savefig(output_path.with_suffix(".pdf"))  # vector copy, sharp at any zoom
        plt.close(fig)
        print(f"[Comparator] Ratio plot saved -> {output_path}")
        return output_path


def main():
    parser = argparse.ArgumentParser(description="Compare two MEGAlib / TGRS energy spectra")
    parser.add_argument("ref", help="Reference .json or .root file")
    parser.add_argument("test", nargs="?", help="Test file (default: same file as ref)")
    parser.add_argument("--ref-hist", default="TGRS", help="Reference histogram name in ROOT file")
    parser.add_argument("--test-hist", default="Sim", help="Test histogram name in ROOT file")
    parser.add_argument("--mode", choices=["ratio", "chi2"], default="ratio")
    parser.add_argument("--norm", default="none",
                        help="none | area | time | width, combinable with + (e.g. time+width)")
    parser.add_argument("--tolerance", type=float, default=0.25,
                        help="ratio mode: allowed fractional deviation per bin (0.25 = 25%%)")
    parser.add_argument("--pass-fraction", type=float, default=0.8,
                        help="ratio mode: fraction of bins that must be within tolerance to pass")
    parser.add_argument("--counts", action="store_true",
                        help="contents are raw counts, not per-keV densities")
    parser.add_argument("--ref-time", type=float)
    parser.add_argument("--test-time", type=float)
    parser.add_argument("--sigma", type=float, default=3.0)
    parser.add_argument("--ref-label", default="Measurements")
    parser.add_argument("--test-label", default="Combined Simulations")
    parser.add_argument("--dpi", type=int, default=400, help="PNG resolution (a PDF is saved too)")
    parser.add_argument("--out", default="comparison", help="Output directory")
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    comp = Comparator(args.ref, args.test or args.ref, out / "comparison_results.json", out,
                      sigma_threshold=args.sigma, normalization=args.norm,
                      ref_time=args.ref_time, test_time=args.test_time,
                      ref_label=args.ref_label, test_label=args.test_label,
                      ref_hist=args.ref_hist, test_hist=args.test_hist, mode=args.mode,
                      tolerance=args.tolerance, pass_fraction=args.pass_fraction, dpi=args.dpi,
                      density=not args.counts)
    for k, v in comp.compare().items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()