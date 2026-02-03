import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from pathlib import Path


class Comparator:

    def __init__(self, ref_json, test_json, output_json,
                 histogram_output_dir, sigma_threshold=3.0):

        self.ref_json = Path(ref_json)
        self.test_json = Path(test_json)
        self.output_json = Path(output_json)
        self.hist_dir = Path(histogram_output_dir)
        self.hist_dir.mkdir(parents=True, exist_ok=True)

        self.sigma_threshold = sigma_threshold

    # ----------------------------------------
    # Load histograms
    # ----------------------------------------
    def load_hist(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Histogram JSON not found: {path}")

        with open(path, "r") as f:
            data = json.load(f)

        # FIXED: ensure correct keys exist
        if "bins" not in data or "edges" not in data:
            raise KeyError(
                f"Histogram file {path} does not contain required keys "
                f"Keys found: {list(data.keys())}"
            )
        counts = np.asarray(data["bins"], dtype = float)
        edges = np.asarray(data["edges"], dtype = float)
        
        return counts, edges

    # ----------------------------------------
    # Compute mu sigma from histogram
    # ----------------------------------------
    def compute_moments(self, counts, bins):
    	total = np.sum(counts)
    	if total == 0:
    		print("[Comparator] WARNING: Histogram contains zero total counts")
    		return 0.0,0.0
    	centers = 0.5 * (edges[:-1] + edges[1:])
    	mean = np.average(centers, weights=counts)
    	variance = np.average((centers - mean) ** 2, weights=counts)
    	sigma = np.sqrt(variance)
    	return mean, sigma
        
    def expand_samples(self, counts, edges):
    	if np.sum(counts) == 0:
    		return np.array([])
    		
    	centers = 0.5 * (edges[:-1] + edges[1:])
    	return np.repeat(centers, counts.astype(int))

    # ----------------------------------------
    # Compare and return dict
    # ----------------------------------------
    def compare(self):

        ref_counts, ref_edges = self.load_hist(self.ref_json)
        test_counts, test_edges = self.load_hist(self.test_json)

        ref_mu, ref_sigma = self.compute_moments(ref_counts, ref_edges)
        test_mu, test_sigma = self.compute_moments(test_counts, test_edges)

        sigma_diff = abs(ref_sigma - test_sigma)
        sigma_sig = sigma_diff / ref_sigma if ref_sigma > 0 else 999
        
        ref_samples = self.expand_samples(ref_counts, ref_edges)
        test_samples = self.expand_samples(test_counts, test_edges)
        
        if len(ref_samples) == 0 or len(test_samples) == 0:
        	ks_stat, ks_p = 0.0, 0.0
        	passed = False
        	print("[Comparator] WARNING: Empty samples -> automatic FAIL")
        else:
        	ks_stat, ks_p = ks_2samp(ref_samples, test_samples)
        	passed = (sigma_sig < self.sigma_threshold) and (ks_p > 0.05)

        
        results = {
            "reference_sigma": ref_sigma,
            "test_sigma": test_sigma,
            "sigma_diff": sigma_diff,
            "sigma_significance": sigma_sig,
            "sigma_threshold": self.sigma_threshold,
            "ks_statistic": float(ks_stat),
            "ks_pvalue": float(ks_p),
            "pass": bool(passed)
        }
        
        #save JSON
        with open(self.output_json, "w") as f:
            json.dump(results, f, indent=4)

        print(f"[Comparator] Results saved → {self.output_json}")
        return results

    # ----------------------------------------
    # Overlay histogram plot
    # ----------------------------------------
    def plot_overlay(self):
        ref_counts, ref_edges = self.load_hist(self.ref_json)
        test_counts, test_edges = self.load_hist(self.test_json)

        plt.figure(figsize=(10, 6))

        plt.hist(ref_edges[:-1], bins=ref_edges, weights=ref_counts,
                 alpha=0.5, label="Reference")
        plt.hist(test_edges[:-1], bins=test_edges, weights=test_counts,
                 alpha=0.5, label="Test")

        plt.xlabel("Energy (keV)")
        plt.ylabel("Counts")
        plt.title("Energy Distribution Comparison")
        plt.legend()

        output_path = self.hist_dir / "energy_overlay.png"
        plt.tight_layout()
        plt.savefig(output_path, dpi=200)
        plt.close()

        print(f"[Comparator] Overlay histogram saved → {output_path}")
        return str(output_path)
