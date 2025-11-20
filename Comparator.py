import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import ks_2samp, anderson_ksamp


class Comparator:
    def __init__(self, dir_path: str):
        # Base directory where test.py lives (algo root or nearby)
        self.base = Path(dir_path).resolve()

        # Try to normalize: if we are already inside run_ref or run_test,
        # go one level up so base points to the algo directory.
        if self.base.name in ["run_ref", "run_test"]:
            self.base = self.base.parent

        self.comp_dir = self.base / "comparison"
        self.comp_dir.mkdir(exist_ok=True)

        print("[Comparator] Base directory:", self.base)
        print("[Comparator] Comparison directory:", self.comp_dir)

    # ---------------------------------------------------
    #  Locate spectrum.C for a given run (run_ref / run_test)
    # ---------------------------------------------------
    def find_spectrum(self, run_name: str) -> Path | None:
        """
        Try a few reasonable locations for spectrum.C:

        1) <base>/<run_name>/results/spectrum.C          (normal case)
        2) <base>/results/spectrum.C                     (if dir_path IS run_ref)
        3) <base>/<run_name>/<run_name>/results/spectrum.C  (nested case)
        """
        candidates = []

        # normal layout: algo/run_ref/results/spectrum.C
        candidates.append(self.base / run_name / "results" / "spectrum.C")

        # if someone passed base=/home/.../algo/run_ref
        if self.base.name == run_name:
            candidates.append(self.base / "results" / "spectrum.C")

        # weird nested layout: algo/run_ref/run_ref/results/spectrum.C
        candidates.append(self.base / run_name / run_name / "results" / "spectrum.C")

        for c in candidates:
            if c.exists():
                print(f"[Comparator] Found {run_name} spectrum at {c}")
                return c

        print(f"[Comparator] Could not find spectrum for {run_name}. Tried:")
        for c in candidates:
            print("   ", c)
        return None

    # ---------------------------------------------------
    #  Extract SetBinContent(bin, value) from spectrum.C
    # ---------------------------------------------------
    def extract_bin_contents(self, spectrum_c_path: Path):
        vals = []
        with open(spectrum_c_path, "r") as f:
            for line in f:
                line = line.strip()
                if "SetBinContent(" not in line:
                    continue
                try:
                    inside = line.split("SetBinContent(")[1].split(")")[0]
                    bin_idx, val = inside.split(",")
                    vals.append(float(val))
                except Exception:
                    continue
        return vals

    # ---------------------------------------------------
    #  Make overlaid histogram plot
    # ---------------------------------------------------
    def plot_hist(self, ref_vals, test_vals):
        plt.figure(figsize=(8, 5))

        bins = np.linspace(0, 2000, 100)  # 100 bins like ROOT

        # Reference spectrum – filled green, thin blue outline
        plt.hist(ref_vals, bins=bins, histtype='stepfilled', alpha=0.4,
                label='Reference', color='green')
        plt.hist(ref_vals, bins=bins, histtype='step',
                linewidth=1.0, color='blue')

        # Test spectrum – filled, different fill but matching style
        plt.hist(test_vals, bins=bins, histtype='stepfilled', alpha=0.2,
                label='Test', color='orange')
        plt.hist(test_vals, bins=bins, histtype='step',
                linewidth=1.0, color='darkorange', linestyle='--')

        plt.xlabel("Energy [keV]", fontsize=12)
        plt.ylabel("Counts / bin", fontsize=12)
        plt.title("Energy Spectrum Comparison", fontsize=14)
        plt.legend(fontsize=10)
        plt.grid(alpha=0.2)
        plt.tight_layout()

        out = self.comp_dir / "energy_comparison.png"
        plt.savefig(out, dpi=300)      # high resolution
        plt.close()
        return str(out)

    # ---------------------------------------------------
    #  Main comparison entry point
    # ---------------------------------------------------
    def compare_energy_hist(self):
        print("[Comparator] Starting energy spectrum comparison...")

        ref_spec = self.find_spectrum("run_ref")
        test_spec = self.find_spectrum("run_test")

        if ref_spec is None or test_spec is None:
            result = {
                "pass": False,
                "details": "Missing spectrum.C for run_ref or run_test."
            }
            json_path = self.comp_dir / "comparison_results.json"
            with open(json_path, "w") as f:
                json.dump({"results": result}, f, indent=4)
            print("[Comparator] Aborting: spectrum.C missing.")
            return str(json_path)

        ref_vals = self.extract_bin_contents(ref_spec)
        test_vals = self.extract_bin_contents(test_spec)

        if len(ref_vals) == 0 or len(test_vals) == 0:
            result = {
                "pass": False,
                "details": "Empty histogram data in one or both spectra."
            }
            json_path = self.comp_dir / "comparison_results.json"
            with open(json_path, "w") as f:
                json.dump({"results": result}, f, indent=4)
            print("[Comparator] Aborting: empty histogram arrays.")
            return str(json_path)

        # stats
        ks_stat, ks_p = ks_2samp(ref_vals, test_vals)
        ad = anderson_ksamp([ref_vals, test_vals])

        mean_diff = float(abs(np.mean(ref_vals) - np.mean(test_vals)))
        max_diff = float(abs(np.max(ref_vals) - np.max(test_vals)))

        passed = (ks_p > 0.05) and (ad.significance_level > 0.05)

        img_path = self.plot_hist(ref_vals, test_vals)

        results = {
            "pass": bool(passed),
            "ks_statistic": float(ks_stat),
            "ks_p_value": float(ks_p),
            "anderson_statistic": float(ad.statistic),
            "anderson_significance": float(ad.significance_level),
            "mean_difference": float(mean_diff),
            "max_difference": float(max_diff),
            "details": "Spectra consistent." if bool(passed) else "Significant difference detected."
        }

        json_path = self.comp_dir / "comparison_results.json"
        with open(json_path, "w") as f:
            json.dump(
                {"results": results, "histograms": {"Energy Spectrum": img_path}},
                f,
                indent=4
            )

        print(f"[Comparator] Done. JSON at {json_path}")
        return str(json_path)
