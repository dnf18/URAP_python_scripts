import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import ks_2samp, anderson_ksamp
import gzip
from Reporter_TGRS import TGRS_Reporter



def extract_energy_list(filepath):
    measurements = Path(filepath)
    if not measurements.exists():
        print(f"[DataFlow] .dat not found: {measurements}")
        return
    data = []
    with open(measurements) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    data.append({
                        "channel":     int(parts[0]),
                        "emin":        float(parts[1]),
                        "emax":        float(parts[2]),
                        "count_rate":  float(parts[3]),
                        "uncertainty": float(parts[4])
                    })
                except ValueError:
                    continue
    return data


def parse_sim_to_json(filepath, output_json):
    events = []
    current_event = None

    with gzip.open(filepath, 'rt') as f:
        for line in f:
            line = line.strip()

            if line.startswith('SE'):
                current_event = {"id": None, "time": None, "energy": None, "triggered": False}

            elif line.startswith('ID') and current_event is not None:
                parts = line.split()
                current_event["id"] = int(parts[1])

            elif line.startswith('TI') and current_event is not None:
                current_event["time"] = float(line.split()[1])

            elif line.startswith('ED') and current_event is not None:
                current_event["energy"] = float(line.split()[1])

            elif line.startswith('HTsim') and current_event is not None:
                current_event["triggered"] = True
                events.append(current_event)
                current_event = None

    with open(output_json, 'w') as f:
        json.dump(events, f, indent=2)

    print(f"Saved {len(events)} triggered events to {output_json}")
    return events


def load_sim_rate(sim_json, emin, emax, bins, bin_width):
    with open(sim_json) as f:
        events = json.load(f)

    obs_time = events[-1]["time"]
    energies = np.array([e["energy"] for e in events if e["energy"] is not None], dtype=float)
    energies = energies[(energies >= emin[0]) & (energies <= emax[-1])]

    hist, _ = np.histogram(energies, bins=bins)
    rate = hist / (obs_time * bin_width)
    return rate


def compute_stats(rate_a, rate_b):
    ref_mu    = np.mean(rate_a)
    ref_sigma = np.std(rate_a)
    test_mu   = np.mean(rate_b)
    test_sigma= np.std(rate_b)

    sigma_diff = abs(ref_sigma - test_sigma)
    sigma_sig  = sigma_diff / ref_sigma if ref_sigma > 0 else np.nan

    ks_stat, ks_p = ks_2samp(rate_a, rate_b)

    try:
        ad_result = anderson_ksamp([rate_a, rate_b])
        ad_stat   = float(ad_result.statistic)
        ad_p      = float(ad_result.significance_level)
    except Exception:
        ad_stat, ad_p = np.nan, np.nan

    passed = bool(ks_p > 0.05 and ad_p > 0.05)

    return {
        "Reference Mean":                        float(ref_mu),
        "Reference Standard Deviation":          float(ref_sigma),
        "Test Mean":                             float(test_mu),
        "Test Standard Deviation":               float(test_sigma),
        "Standard Deviation Difference":         float(sigma_diff),
        "Relative Standard Deviation Difference":float(sigma_sig),
        "Kolmogorov-Smirnov statistic":          float(ks_stat),
        "Kolmogorov-Smirnov pvalue":             float(ks_p),
        "Anderson-Darling Statistic":            float(ad_stat),
        "Anderson-Darling p-value":              float(ad_p),
        "pass":                                  passed,
    }


def graph(tgrs_json, sim_a_json, sim_b_json, output_dir="."):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Load real TGRS data
    with open(tgrs_json) as f:
        data = json.load(f)

    emin        = np.array([d["emin"]        for d in data])
    emax        = np.array([d["emax"]        for d in data])
    count_rate  = np.array([d["count_rate"]  for d in data])
    uncertainty = np.array([d["uncertainty"] for d in data])
    emid        = (emin + emax) / 2
    bin_width   = emax - emin
    bins        = np.append(emin, emax[-1])

    # Load sim rates
    rate_a = load_sim_rate(sim_a_json, emin, emax, bins, bin_width)
    rate_b = load_sim_rate(sim_b_json, emin, emax, bins, bin_width)

    # Ratios
    ratio_a = np.where(rate_a > 0, count_rate / rate_a, np.nan)
    ratio_b = np.where(rate_b > 0, count_rate / rate_b, np.nan)

    # --------------------------------------------------
    # Figure 1 - Sim A vs Measurements
    # --------------------------------------------------
    plot1_path = str(output_dir / 'plot1_sim_a.png')
    plt.figure(figsize=(12, 6))
    plt.step(emid, count_rate, where='mid', color='black', linewidth=1.5, label='Measurements')
    plt.errorbar(emid, count_rate, yerr=uncertainty, fmt='none', color='black', alpha=0.3)
    plt.step(emid, rate_a, where='mid', color='blue', linewidth=1, alpha=0.7, label='Cosmic Photons A')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Hits (cts/keV/s)')
    plt.ylim(bottom=1e-7)
    plt.title('Cosmic Photons A vs Measurements')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot1_path, dpi=150)
    plt.close()
    print(f"Saved {plot1_path}")

    stats_a = compute_stats(rate_a, count_rate)
    stats_a["overlay_plot"] = plot1_path
    with open(output_dir / 'reference_stats.json', 'w') as f:
        json.dump(stats_a, f, indent=2)
    print("Saved reference_stats.json")

    # --------------------------------------------------
    # Figure 2 - Sim B vs Measurements
    # --------------------------------------------------
    plot2_path = str(output_dir / 'plot2_sim_b.png')
    plt.figure(figsize=(12, 6))
    plt.step(emid, count_rate, where='mid', color='black', linewidth=1.5, label='Measurements')
    plt.errorbar(emid, count_rate, yerr=uncertainty, fmt='none', color='black', alpha=0.3)
    plt.step(emid, rate_b, where='mid', color='red', linewidth=1, alpha=0.7, label='Cosmic Photons B')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Hits (cts/keV/s)')
    plt.ylim(bottom=1e-7)
    plt.title('Cosmic Photons B vs Measurements')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot2_path, dpi=150)
    plt.close()
    print(f"Saved {plot2_path}")

    stats_b = compute_stats(rate_b, count_rate)
    stats_b["overlay_plot"] = plot2_path
    with open(output_dir / 'test_stats.json', 'w') as f:
        json.dump(stats_b, f, indent=2)
    print("Saved test_stats.json")

    # --------------------------------------------------
    # Figure 3 - Ratio A vs Measurements
    # --------------------------------------------------
    plot3_path = str(output_dir / 'plot3_ratio_a.png')
    plt.figure(figsize=(12, 6))
    plt.step(emid, ratio_a, where='mid', color='blue', linewidth=1.5, label='Measurements / Sim A')
    plt.axhline(1.0, color='black', linestyle='--', linewidth=1, label='Perfect agreement')
    plt.xscale('log')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Ratio (Measurements / Sim)')
    plt.title('Ratio - Measurements vs Cosmic Photons A')
    plt.ylim(0, 3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot3_path, dpi=150)
    plt.close()
    print(f"Saved {plot3_path}")

    stats_ratio_a = compute_stats(rate_a, count_rate)
    stats_ratio_a["overlay_plot"] = plot3_path
    with open(output_dir / 'reference_ratio_stats.json', 'w') as f:
        json.dump(stats_ratio_a, f, indent=2)
    print("Saved reference_ratio_stats.json")

    # --------------------------------------------------
    # Figure 4 - Ratio B vs Measurements
    # --------------------------------------------------
    plot4_path = str(output_dir / 'plot4_ratio_b.png')
    plt.figure(figsize=(12, 6))
    plt.step(emid, ratio_b, where='mid', color='red', linewidth=1.5, label='Measurements / Sim B')
    plt.axhline(1.0, color='black', linestyle='--', linewidth=1, label='Perfect agreement')
    plt.xscale('log')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Ratio (Measurements / Sim)')
    plt.title('Ratio - Measurements vs Cosmic Photons B')
    plt.ylim(0, 3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot4_path, dpi=150)
    plt.close()
    print(f"Saved {plot4_path}")

    stats_ratio_b = compute_stats(rate_b, count_rate)
    stats_ratio_b["overlay_plot"] = plot4_path
    with open(output_dir / 'test_ratio_stats.json', 'w') as f:
        json.dump(stats_ratio_b, f, indent=2)
    print("Saved test_ratio_stats.json")

    # --------------------------------------------------
    # Figure 5 - All together
    # --------------------------------------------------
    plot5_path = str(output_dir / 'plot5_all.png')
    plt.figure(figsize=(12, 6))
    plt.step(emid, count_rate, where='mid', color='black', linewidth=1.5, label='Measurements')
    plt.errorbar(emid, count_rate, yerr=uncertainty, fmt='none', color='black', alpha=0.3)
    plt.step(emid, rate_a, where='mid', color='blue', linewidth=1, alpha=0.7, label='Cosmic Photons A')
    plt.step(emid, rate_b, where='mid', color='red', linewidth=1, alpha=0.7, label='Cosmic Photons B')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Hits (cts/keV/s)')
    plt.ylim(bottom=1e-7)
    plt.title('Cosmic Photons A & B vs Measurements')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot5_path, dpi=150)
    plt.close()
    print(f"Saved {plot5_path}")

    stats_all = compute_stats(rate_a, rate_b)
    stats_all["overlay_plot"] = plot5_path
    with open(output_dir / 'combined_stats.json', 'w') as f:
        json.dump(stats_all, f, indent=2)
    print("Saved combined_stats.json")

    plt.show()


# --- Main ---
data = extract_energy_list('TGRS.dat')
with open('TGRS.json', 'w') as f:
    json.dump(data, f, indent=2)
print(f"Saved {len(data)} bins to TGRS.json")

parse_sim_to_json('a/CosmicPhotons.inc1.id1.sim.gz', 'CosmicPhotons_a.json')
parse_sim_to_json('b/CosmicPhotons.inc1.id1.sim.gz', 'CosmicPhotons_b.json')

graph('TGRS.json', 'CosmicPhotons_a.json', 'CosmicPhotons_b.json', output_dir='results')
reporter = TGRS_Reporter(
    config_json="steering_config.json",
    output_dir="./reports"
)

pdf_path = reporter.generate_pdf(stats_dir="./results")

