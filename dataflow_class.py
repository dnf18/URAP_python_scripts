import os
import json
import subprocess
import numpy as np
from pathlib import Path


class DataFlow:
    def __init__(self, json_path: str):

        self.json_path = Path(json_path)
        if not self.json_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.json_path}")

        # load config
        with open(self.json_path, "r") as f:
            self.config = json.load(f)

        # extract fields
        self.cosima_file = self.config["cosima_file"]
        self.geometry_file = self.config["geometry_file"]
        self.revan_cfg = self.config["revan_output"]
        self.mimrec_cfg = self.config["mimrec_output"]
        self.energy_cut = self.config.get("energy_cut", [10, 2000])
        self.max_events = self.config.get("max_events", 100000)

        # directory for this run
        self.dir = self.json_path.parent

        # output directory
        self.output_dir = self.dir / "results"
        self.output_dir.mkdir(exist_ok=True)

    # --------------------------------------------------------------

    def run_command(self, cmd):
        print(f">>> {' '.join(cmd)}\n")
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        for line in process.stdout:
            print(line, end="")

        process.wait()

        if process.returncode != 0:
            raise RuntimeError(f"Command failed: {' '.join(cmd)}")

    # --------------------------------------------------------------

    def run_simulation(self):
        print("==== Step 1: Simulation ====")
        self.run_command(["cosima", self.cosima_file])

    # --------------------------------------------------------------

    def run_reconstruction(self):
        print("==== Step 2: Reconstruction ====")

        base = os.path.splitext(self.cosima_file)[0]
        sim_file = base + ".inc1.id1.sim.gz"

        self.run_command([
            "revan",
            "-c", self.revan_cfg,
            "-g", self.geometry_file,
            "-f", sim_file,
            "-a", "-n"
        ])

    # --------------------------------------------------------------

    def run_spectrum(self):
        print("==== Step 3: Spectrum Creation ====")

        base = os.path.splitext(self.cosima_file)[0]
        tra_file = base + ".inc1.id1.tra.gz"
        output_spectrum = self.output_dir / "spectrum.C"

        self.run_command([
            "mimrec",
            "-c", self.mimrec_cfg,
            "-g", self.geometry_file,
            "-f", tra_file,
            "-s",
            "-o", str(output_spectrum)
        ])

    # --------------------------------------------------------------

    def extract_energy_list(self):
        """
        Reads spectrum.C and extracts energies, saving to {reference,test}_energy.txt
        """

        spectrum_file = self.output_dir / "spectrum.C"
        if not spectrum_file.exists():
            print(f"[DataFlow] No spectrum.C found at {spectrum_file}")
            return

        # choose file name
        if "run_ref" in str(self.output_dir):
            out_file = self.output_dir / "reference_energy.txt"
        else:
            out_file = self.output_dir / "test_energy.txt"

        energies = []

        with open(spectrum_file, "r") as f:
            for line in f:
                line = line.strip()

                # Energy[i] = 215.0;
                if "Energy" in line and "=" in line:
                    try:
                        energies.append(float(line.split("=")[1].replace(";", "").strip()))
                        continue
                    except:
                        pass

                # hEnergy->Fill(215.0);
                if "Fill" in line:
                    try:
                        val = float(line.split("(")[1].split(")")[0])
                        energies.append(val)
                        continue
                    except:
                        pass

                # push_back(215.0)
                if "push_back" in line:
                    try:
                        val = float(line.split("(")[1].split(")")[0])
                        energies.append(val)
                        continue
                    except:
                        pass

        with open(out_file, "w") as f:
            for e in energies:
                f.write(f"{e}\n")

        print(f"[DataFlow] Extracted {len(energies)} energies → {out_file}")

        return energies

    # --------------------------------------------------------------
    # NEW: compute histogram (fixes bin mismatch problems)
    # --------------------------------------------------------------

    def save_histogram(self, energies):
        """
        Converts energy list into histogram + bins and saves them:
        - energy_hist.npy
        - energy_bins.npy
        """

        energies = np.array(energies, dtype=float)

        hist, bins = np.histogram(energies, bins=100)

        np.save(self.output_dir / "energy_hist.npy", hist)
        np.save(self.output_dir / "energy_bins.npy", bins)

        print(f"[DataFlow] Saved histogram + bins in {self.output_dir}")

    # --------------------------------------------------------------

    def organize_results(self):
        base_prefix = os.path.splitext(os.path.basename(self.cosima_file))[0]
        safe_ext = (".sim.gz", ".tra.gz", ".root", ".C", ".txt", ".dat")

        for file in os.listdir('.'):
            if file.startswith(base_prefix) and file.endswith(safe_ext):
                os.rename(file, self.output_dir / file)

    # --------------------------------------------------------------

    def run_full_pipeline(self):
        os.chdir(self.dir)

        self.run_simulation()
        self.run_reconstruction()
        self.run_spectrum()
        self.organize_results()

        energies = self.extract_energy_list()

        if energies:
            self.save_histogram(energies)
