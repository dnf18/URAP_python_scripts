"""
TGRS_Comparator

Compares two MEGAlib builds (A = reference, B = test) against each other and
against the measured TGRS spectrum (TGRS.dat).

What it does, following MEGAlib's own TGRS example (TGRS.cxx):
  * The measured spectrum is in ph/cm2/s/keV, so it is multiplied by the
    35 cm2 detector area (TGRS.cxx: TGRSHist->Scale(35)) to get cts/keV/s.
  * Each simulated event is filled at the SUM of its HTsim hit energies
    (TGRS.cxx uses MSimEvent::GetREnergy(), which is exactly that sum).
  * Each simulation component (photons, protons, ...) is normalised by ITS OWN
    observation time and bin width, and the components are then added together,
    like the "Sim" histogram in TGRS.cxx.
  * A component can consist of many files (mcosima writes one per instance).
    Counts are summed over the files and so are their observation times.
  * Observation time per file = last event time minus the start time, where the
    start is the "TB" header value if present, otherwise the first event time.
    (This is how MFileEvents::GetObservationTime works.)

Not done here: detector energy resolution. TGRS.cxx applies it unless it is run
with -p; this script uses the raw simulated hit energies. This does not matter
for A-versus-B, but line shapes against the measurements will differ.

Usage:
    python3 TGRS_Comparator.py --tgrs TGRS.dat \\
        --ref  "photons=a/CosmicPhotons.p1.inc*.id1.sim.gz" \\
        --ref  "protons=a/CosmicProtons.Decays.p1.inc*.id1.sim.gz" \\
        --test "photons=b/CosmicPhotons.p1.inc*.id1.sim.gz" \\
        --test "protons=b/CosmicProtons.Decays.p1.inc*.id1.sim.gz" \\
        --ref-label "Geant4 10.2 (main)" --test-label "Geant4 11.4 (dc6)" \\
        --out results --config steering_config.json --reports reports

Each --ref / --test is NAME=GLOB. Use --ref/--test as often as you have
components (a photons-only run is fine). Quote the globs.
"""

import argparse
import glob
import gzip
import json
import zlib
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")   # no display needed (works on Savio)
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_dist

DETECTOR_AREA_CM2 = 35.0     # TGRS.cxx: TGRSHist->Scale(35)
SIGMA_THRESHOLD = 3.0
P_THRESHOLD = 0.05


# ------------------------------------------------------------------ measurements

def load_tgrs(filepath):
    """
    Reads TGRS.dat the way TGRS.cxx does: everything after the line starting
    with 'DATA', keeping lines with exactly 5 tokens: ch emin emax flux +-.
    Flux is in ph/cm2/s/keV and is multiplied by the 35 cm2 area to give
    cts/keV/s; the uncertainty is scaled the same way.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"TGRS data file not found: {filepath}")

    rows = []
    started = False
    with open(filepath) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "DATA":
                started = True
            if started and len(parts) == 5:
                try:
                    rows.append([float(parts[1]), float(parts[2]),
                                 float(parts[3]), float(parts[4])])
                except ValueError:
                    continue

    if not rows:
        raise ValueError(f"No data rows found in {filepath}")

    arr = np.array(rows)
    emin, emax = arr[:, 0], arr[:, 1]
    if np.any(emin[1:] != emax[:-1]):
        print("[TGRS] Warning: energy bins are not contiguous; "
              "bin edges are taken from the lower edges (as TGRS.cxx does).")
    edges = np.append(emin, emax[-1])
    rate = arr[:, 2] * DETECTOR_AREA_CM2
    err = arr[:, 3] * DETECTOR_AREA_CM2
    return edges, rate, err


# ------------------------------------------------------------------ simulations

def read_sim_file(path):
    """
    Streams one (gzipped) .sim file.
    Returns (energies, duration, truncated):
      energies  - array, one entry per event with >= 1 HTsim hit = sum of hit energies [keV]
      duration  - observation time in seconds, or None if the file has no events
      truncated - True if the gzip stream ended early (killed job); what was read is kept
    """
    path = str(path)
    opener = gzip.open if path.endswith(".gz") else open

    energies = []
    header_start = None
    first_time = None
    last_time = None

    in_event = False
    event_energy = 0.0
    event_has_hit = False
    truncated = False

    try:
        with opener(path, "rt") as f:
            for line in f:
                parts = line.split(None, 1)
                if not parts:
                    continue
                tag = parts[0]

                if tag == "SE":
                    if in_event and event_has_hit:
                        energies.append(event_energy)
                    in_event = True
                    event_energy = 0.0
                    event_has_hit = False

                elif tag == "TI":
                    t = float(parts[1].split()[0])
                    if in_event:
                        if first_time is None:
                            first_time = t
                        last_time = t

                elif tag == "TB" and not in_event and len(parts) > 1:
                    header_start = float(parts[1].split()[0])

                elif tag == "HTsim" and in_event and len(parts) > 1:
                    # "HTsim 4; x; y; z; energy; ..." -> energy is field 4
                    fields = parts[1].split(";")
                    event_energy += float(fields[4])
                    event_has_hit = True

                elif tag == "EN":
                    if in_event and event_has_hit:
                        energies.append(event_energy)
                    in_event = False
                    break
            else:
                # reached end of file without an EN line
                if in_event and event_has_hit:
                    energies.append(event_energy)
    except (EOFError, OSError, zlib.error):
        # Truncated gzip (job killed). The event being read is incomplete: drop it.
        truncated = True

    if first_time is None or last_time is None:
        return np.array([]), None, truncated

    start = header_start if header_start is not None else first_time
    duration = last_time - start
    if duration <= 0:
        return np.array(energies, dtype=float), None, truncated
    return np.array(energies, dtype=float), duration, truncated


def component_rate(pattern, edges):
    """
    Reads every file matching `pattern` (a glob), sums counts and observation
    times, and returns (rate, error, info) in cts/keV/s.
    The error is Poisson, with a floor of one count so empty bins keep a
    non-zero uncertainty.
    """
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files match: {pattern}")

    counts = np.zeros(len(edges) - 1)
    total_time = 0.0
    used, skipped, truncated_files = 0, [], []

    for p in files:
        energies, duration, truncated = read_sim_file(p)
        if duration is None:
            skipped.append(p)       # empty file, include-stub, or no events
            continue
        counts += np.histogram(energies, bins=edges)[0]
        total_time += duration
        used += 1
        if truncated:
            truncated_files.append(p)

    if used == 0 or total_time <= 0:
        raise ValueError(f"No usable events in files matching: {pattern}")

    if skipped:
        print(f"[Sim] {len(skipped)} file(s) skipped (no events): {skipped[:3]}"
              f"{' ...' if len(skipped) > 3 else ''}")
    if truncated_files:
        print(f"[Sim] WARNING: {len(truncated_files)} file(s) were truncated "
              f"(job killed?); partial data used.")

    width = np.diff(edges)
    rate = counts / (total_time * width)
    err = np.sqrt(np.maximum(counts, 1.0)) / (total_time * width)
    info = {
        "pattern": pattern,
        "files_used": used,
        "files_skipped": len(skipped),
        "files_truncated": len(truncated_files),
        "observation_time_s": total_time,
        "events_in_range": int(counts.sum()),
    }
    return rate, err, info


def side_rate(components, edges):
    """Adds the components of one side (e.g. photons + protons), like TGRS.cxx."""
    total = np.zeros(len(edges) - 1)
    var = np.zeros(len(edges) - 1)
    infos = {}
    for name, pattern in components.items():
        rate, err, info = component_rate(pattern, edges)
        total += rate
        var += err ** 2
        infos[name] = info
        print(f"[Sim] {name}: {info['files_used']} file(s), "
              f"{info['observation_time_s']:.1f} s, {info['events_in_range']} events")
    return total, np.sqrt(var), infos


# ------------------------------------------------------------------ statistics

def _chi2_block(pulls, used, prefix=""):
    chi2_val = float(np.sum(pulls[used] ** 2))
    ndf = int(np.sum(used))
    p = float(chi2_dist.sf(chi2_val, ndf)) if ndf > 0 else float("nan")
    return {
        f"{prefix}Chi2": chi2_val,
        f"{prefix}NDF": ndf,
        f"{prefix}Chi2/NDF": chi2_val / ndf if ndf > 0 else float("nan"),
        f"{prefix}Chi2 p-value": p,
    }, p


def compare_rates(x, ex, y, ey, edges, name_x, name_y):
    """Per-bin pulls and chi2 between two rate spectra, with propagated errors."""
    sigma = np.sqrt(ex ** 2 + ey ** 2)
    used = (sigma > 0) & ~((x == 0) & (y == 0))
    pulls = np.zeros_like(x)
    pulls[used] = (x[used] - y[used]) / sigma[used]

    block, p = _chi2_block(pulls, used)

    worst = int(np.argmax(np.abs(pulls)))
    width = np.diff(edges)
    cx, cy = np.cumsum(x * width), np.cumsum(y * width)
    ks = (float(np.max(np.abs(cx / cx[-1] - cy / cy[-1])))
          if cx[-1] > 0 and cy[-1] > 0 else float("nan"))

    stats = {
        "Compared": f"{name_x} vs {name_y}",
        "Bins compared": int(np.sum(used)),
        f"Mean rate, {name_x} (cts/keV/s)": float(np.mean(x)),
        f"Mean rate, {name_y} (cts/keV/s)": float(np.mean(y)),
        **block,
        f"Bins beyond {SIGMA_THRESHOLD:g} sigma": int(np.sum(np.abs(pulls) > SIGMA_THRESHOLD)),
        "Worst bin energy range (keV)": [float(edges[worst]), float(edges[worst + 1])],
        "Worst bin pull": float(pulls[worst]),
        "KS distance (cumulative spectra)": ks,
        "pass": bool(p > P_THRESHOLD),
    }
    return stats, pulls


def compare_ratio(meas, emeas, sim, esim, edges, name_sim):
    """Statistics of measurement/simulation; perfect agreement is a ratio of 1."""
    ok = (sim > 0) & (meas > 0)
    ratio = np.full_like(meas, np.nan)
    ratio_err = np.full_like(meas, np.nan)
    ratio[ok] = meas[ok] / sim[ok]
    ratio_err[ok] = ratio[ok] * np.sqrt((emeas[ok] / meas[ok]) ** 2 + (esim[ok] / sim[ok]) ** 2)

    used = ok & (ratio_err > 0)
    pulls = np.zeros_like(meas)
    pulls[used] = (ratio[used] - 1.0) / ratio_err[used]
    block, p = _chi2_block(pulls, used, prefix="Chi2 about 1: ")
    # keep the usual key names for the reporter
    block = {k.replace("Chi2 about 1: ", ""): v for k, v in block.items()}

    r = ratio[ok]
    stats = {
        "Compared": f"Measurements / {name_sim}",
        "Bins with a ratio": int(np.sum(ok)),
        "Median ratio": float(np.median(r)) if r.size else float("nan"),
        "Mean ratio": float(np.mean(r)) if r.size else float("nan"),
        "Std of ratio": float(np.std(r)) if r.size else float("nan"),
        "Fraction of bins within 20% of 1": float(np.mean(np.abs(r - 1.0) < 0.2)) if r.size else float("nan"),
        **block,
        "pass": bool(p > P_THRESHOLD),
    }
    return stats, ratio


# ------------------------------------------------------------------ plotting

def _floor(*arrays):
    vals = np.concatenate([a[a > 0] for a in arrays if np.any(a > 0)] or [np.array([1e-7])])
    return float(vals.min()) / 2.0


def plot_spectra(path, edges, meas, emeas, series, title):
    """series: list of (label, rate, color)."""
    emid = 0.5 * (edges[:-1] + edges[1:])
    plt.figure(figsize=(12, 6))
    plt.step(emid, meas, where="mid", color="black", linewidth=1.5, label="Measurements (x35 cm2)")
    plt.errorbar(emid, meas, yerr=emeas, fmt="none", color="black", alpha=0.3)
    for label, rate, color in series:
        plt.step(emid, rate, where="mid", color=color, linewidth=1, alpha=0.8, label=label)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Hits (cts/keV/s)")
    plt.ylim(bottom=_floor(meas, *[s[1] for s in series]))
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"Saved {path}")


def plot_ratio(path, edges, ratio, label, color, title):
    emid = 0.5 * (edges[:-1] + edges[1:])
    plt.figure(figsize=(12, 6))
    plt.step(emid, ratio, where="mid", color=color, linewidth=1.5, label=label)
    plt.axhline(1.0, color="black", linestyle="--", linewidth=1, label="Perfect agreement")
    plt.xscale("log")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Ratio (Measurements / Sim)")
    plt.title(title)
    plt.ylim(0, 3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"Saved {path}")


# ------------------------------------------------------------------ driver

def run_comparison(tgrs_file, ref_components, test_components, output_dir,
                   ref_label="Reference (A)", test_label="Test (B)"):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    edges, meas, emeas = load_tgrs(tgrs_file)
    print(f"[TGRS] {len(meas)} bins, {edges[0]:g} - {edges[-1]:g} keV")

    print(f"\n[{ref_label}]")
    rate_a, err_a, info_a = side_rate(ref_components, edges)
    print(f"\n[{test_label}]")
    rate_b, err_b, info_b = side_rate(test_components, edges)

    def save(name, stats, plot):
        stats["overlay_plot"] = str(plot)
        with open(output_dir / name, "w") as f:
            json.dump(stats, f, indent=2)
        print(f"Saved {name}")

    # 1/2: each simulation against the measurements
    p1 = output_dir / "plot1_sim_a.png"
    plot_spectra(p1, edges, meas, emeas, [(ref_label, rate_a, "blue")],
                 f"{ref_label} vs Measurements")
    s1, _ = compare_rates(rate_a, err_a, meas, emeas, edges, ref_label, "Measurements")
    save("reference_stats.json", s1, p1)

    p2 = output_dir / "plot2_sim_b.png"
    plot_spectra(p2, edges, meas, emeas, [(test_label, rate_b, "red")],
                 f"{test_label} vs Measurements")
    s2, _ = compare_rates(rate_b, err_b, meas, emeas, edges, test_label, "Measurements")
    save("test_stats.json", s2, p2)

    # 3/4: measurement / simulation ratios
    s3, ratio_a = compare_ratio(meas, emeas, rate_a, err_a, edges, ref_label)
    p3 = output_dir / "plot3_ratio_a.png"
    plot_ratio(p3, edges, ratio_a, f"Measurements / {ref_label}", "blue",
               f"Ratio - Measurements vs {ref_label}")
    save("reference_ratio_stats.json", s3, p3)

    s4, ratio_b = compare_ratio(meas, emeas, rate_b, err_b, edges, test_label)
    p4 = output_dir / "plot4_ratio_b.png"
    plot_ratio(p4, edges, ratio_b, f"Measurements / {test_label}", "red",
               f"Ratio - Measurements vs {test_label}")
    save("test_ratio_stats.json", s4, p4)

    # 5: everything together, and A versus B (the actual build-to-build check)
    p5 = output_dir / "plot5_all.png"
    plot_spectra(p5, edges, meas, emeas,
                 [(ref_label, rate_a, "blue"), (test_label, rate_b, "red")],
                 f"{ref_label} & {test_label} vs Measurements")
    s5, _ = compare_rates(rate_a, err_a, rate_b, err_b, edges, ref_label, test_label)
    save("combined_stats.json", s5, p5)

    emid = 0.5 * (edges[:-1] + edges[1:])
    with open(output_dir / "binned_rates.json", "w") as f:
        json.dump({
            "energy_mid_keV": emid.tolist(),
            "measurements": meas.tolist(), "measurements_err": emeas.tolist(),
            "reference": rate_a.tolist(), "reference_err": err_a.tolist(),
            "test": rate_b.tolist(), "test_err": err_b.tolist(),
            "reference_info": info_a, "test_info": info_b,
        }, f, indent=2)
    print("Saved binned_rates.json")

    return {"reference": s1, "test": s2, "ref_ratio": s3, "test_ratio": s4, "ref_vs_test": s5}


def _parse_components(items, flag):
    comps = {}
    for item in items:
        if "=" not in item:
            raise SystemExit(f"{flag} needs NAME=GLOB, got: {item}")
        name, pattern = item.split("=", 1)
        comps[name.strip()] = pattern.strip()
    if not comps:
        raise SystemExit(f"At least one {flag} NAME=GLOB is required")
    return comps


def main():
    ap = argparse.ArgumentParser(description="Compare two MEGAlib builds against TGRS data")
    ap.add_argument("--tgrs", default="TGRS.dat", help="TGRS.dat measurement file")
    ap.add_argument("--ref", action="append", default=[], metavar="NAME=GLOB",
                    help="Reference-build component, e.g. photons=a/CosmicPhotons.p1.inc*.id1.sim.gz")
    ap.add_argument("--test", action="append", default=[], metavar="NAME=GLOB",
                    help="Test-build component")
    ap.add_argument("--ref-label", default="Reference (A)")
    ap.add_argument("--test-label", default="Test (B)")
    ap.add_argument("--out", default="results", help="Directory for plots and stats")
    ap.add_argument("--config", default="steering_config.json",
                    help="Config JSON shown in the PDF (optional)")
    ap.add_argument("--reports", default="./reports", help="Directory for the PDF")
    ap.add_argument("--no-report", action="store_true", help="Skip the PDF")
    args = ap.parse_args()

    run_comparison(args.tgrs,
                   _parse_components(args.ref, "--ref"),
                   _parse_components(args.test, "--test"),
                   args.out, args.ref_label, args.test_label)

    if not args.no_report:
        from TGRS_Reporter import TGRS_Reporter   # imported late: needs reportlab
        reporter = TGRS_Reporter(config_json=args.config, output_dir=args.reports)
        reporter.generate_pdf(stats_dir=args.out)


if __name__ == "__main__":
    main()
