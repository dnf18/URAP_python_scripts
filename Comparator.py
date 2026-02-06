import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from pathlib import Path

"""""
This is the Comparator Class

Using the energy JSONs, we found the Means, Standard Devations(STD), STD difference, Relative STD Difference, and Kolmogorov-Smirnov statistic and pvalue

We also created an overlay histogram of the spectrums

All relevant file paths/data are then saved on a dictionary and sent back to Supervisor
"""""

class Comparator:

    def __init__(self, ref_json, test_json, output_json,
                 histogram_output_dir, sigma_threshold=3.0):

        self.ref_json = Path(ref_json)
        self.test_json = Path(test_json)
        self.output_json = Path(output_json)
        self.hist_dir = Path(histogram_output_dir)
        self.hist_dir.mkdir(parents=True, exist_ok=True)

        self.sigma_threshold = sigma_threshold

    #Load histograms

    def load_hist(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Histogram JSON not found: {path}")

        with open(path, "r") as f:
            data = json.load(f)

        #ensure correct keys exist
        if "bins" not in data or "edges" not in data:
            raise KeyError(
                f"Histogram file {path} does not contain required keys "
                f"'bins' and 'edges'. Keys found: {list(data.keys())}"
            )

        return np.array(data["bins"]), np.array(data["edges"])


    #compute mu sigma from histogram
    
    def compute_moments(self, counts, bins):
        if np.sum(counts) == 0:
            print("[Comparator] WARNING: Histogram contains zero total counts")
            return 0.0, 0.0

        centers = 0.5 * (bins[:-1] + bins[1:])
        mean = np.average(centers, weights=counts)
        variance = np.average((centers - mean) ** 2, weights=counts)
        sigma = np.sqrt(variance)
        return mean, sigma

    #compare and return dict
   
    def compare(self):

        ref_counts, ref_bins = self.load_hist(self.ref_json)
        test_counts, test_bins = self.load_hist(self.test_json)

        ref_mu, ref_sigma = self.compute_moments(ref_counts, ref_bins)
        test_mu, test_sigma = self.compute_moments(test_counts, test_bins)

        sigma_diff = abs(ref_sigma - test_sigma)
        sigma_sig = sigma_diff / ref_sigma if ref_sigma > 0 else 999

        ks_stat, ks_p = ks_2samp(ref_counts, test_counts)

        passed = bool((sigma_sig < self.sigma_threshold) and (ks_p > 0.05))

        overlay_path = self.hist_dir / "energy_overlay.png"

        results = {
            "Reference Mean": float(ref_mu),
            "Reference Standard Deviation": float(ref_sigma),
            "Test Mean": float(test_mu),
            "Test Standard Deviation": float(test_sigma),
            "Standard Deviation Difference": float(sigma_diff),
            "Relative Standard Deviation Difference": float(sigma_sig),
            "Kolmogorov-Smirnov statistic": float(ks_stat),
            "Kolmogorov-Smirnov pvalue": float(ks_p),
            "pass": bool(passed),
            "overlay_plot": str(overlay_path)
        }

        with open(self.output_json, "w") as f:
            json.dump(results, f, indent=4)

        print(f"[Comparator] Results saved → {self.output_json}")
        return results

    #overlay histogram plot
 
    def plot_overlay(self):
        ref_counts, ref_bins = self.load_hist(self.ref_json)
        test_counts, test_bins = self.load_hist(self.test_json)

        plt.figure(figsize=(10, 6))
        #bins from the JSON file which a data scrapped from spectrum.C file
        plt.hist(ref_bins[:-1], bins=ref_bins, weights=ref_counts,
                 alpha=0.5, label="Reference")
        plt.hist(test_bins[:-1], bins=test_bins, weights=test_counts,
                 alpha=0.5, label="Test")

        plt.xlabel("Energy (keV)")
        plt.ylabel("Counts")
        plt.title("Energy Distribution Comparison")
        plt.legend()

        output_path = self.hist_dir / "energy_overlay.png"
        plt.savefig(output_path, dpi=200)
        plt.close()

        print(f"[Comparator] Overlay histogram saved → {output_path}")
        return str(output_path)
