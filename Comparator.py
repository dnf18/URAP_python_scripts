import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import ks_2samp


class Comparator:
    def __init__(self, dir_path: str):
        # Base directory where runs live
        self.base = Path(dir_path).resolve()
        if self.base.name in ["run_ref", "run_test"]:
            self.base = self.base.parent

        self.comp_dir = self.base / "comparison"
        self.comp_dir.mkdir(exist_ok=True)

        print("[Comparator] Base directory:", self.base)
        print("[Comparator] Comparison directory:", self.comp_dir)

    # ---------------------------------------------------
    # Locate spectrum.C for a given run
    # ---------------------------------------------------
    def find_spectrum(self, run_name: str) -> Path | None:
        candidates = [
            self.base / run_name / "results" / "spectrum.C",
            self.base / run_name / run_name / "results" / "spectrum.C"
        ]
        if self.base.name == run_name:
            candidates.append(self.base / "results" / "spectrum.C")

        for c in candidates:
            if c.exists():
                print(f"[Comparator] Found {run_name} spectrum at {c}")
                return c

        print(f"[Comparator] Could not find spectrum for {run_name}. Tried:")
        for c in candidates:
            print("   ", c)
        return None

    # ---------------------------------------------------
    # Extract energies from spectrum.C
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
    # Make overlaid histogram plot (step + fill)
    # ---------------------------------------------------
    def plot_hist(self, ref_vals, test_vals, bins=40):
        plt.figure(figsize=(7, 4))

        # Compute histograms with the same bin edges
        hist_ref, bin_edges = np.histogram(ref_vals, bins=bins)
        hist_test, _ = np.histogram(test_vals, bins=bin_edges)

        # Step plot with fill_between for overlay
        plt.step(bin_edges[:-1], hist_ref, where='mid', label="Reference", color='blue')
        plt.fill_between(bin_edges[:-1], hist_ref, step='mid', alpha=0.3, color='blue')

        plt.step(bin_edges[:-1], hist_test, where='mid', label="Test", color='orange')
        plt.fill_between(bin_edges[:-1], hist_test, step='mid', alpha=0.3, color='orange')

        plt.xlabel("Counts per bin")
        plt.ylabel("Bin frequency")
        plt.title("Energy Spectrum Comparison")
        plt.legend()

        out = self.comp_dir / "energy_comparison.png"
        plt.tight_layout()
        plt.savefig(out)
        plt.close()
        return str(out)

    # ---------------------------------------------------
    # Main comparison entry point
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

        # Statistics
        ks_stat, ks_p = ks_2samp(ref_vals, test_vals)

        mean_diff = float(abs(np.mean(ref_vals) - np.mean(test_vals)))
        sigma_diff = float(abs(np.std(ref_vals) - np.std(test_vals)))
        max_diff = float(abs(np.max(ref_vals) - np.max(test_vals)))

        # Define a "pass" threshold: KS p-value > 0.05 and sigma difference < 10%
        passed = (ks_p > 0.05) and (sigma_diff / np.std(ref_vals) < 0.1)

        img_path = self.plot_hist(ref_vals, test_vals)

        results = {
            "pass": bool(passed),
            "ks_statistic": float(ks_stat),
            "ks_p_value": float(ks_p),
            "mean_difference": float(mean_diff),
            "sigma_difference": float(sigma_diff),
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
