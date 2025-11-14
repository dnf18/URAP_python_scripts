#!/usr/bin/env python3

import os
import json
import subprocess
from pathlib import Path


class DataFlow:
    def __init__(self, json_path: str):

        """
        Initialize DataFlow with a path to a steering_config.json file.
        Loads all necessary simulation parameters from the JSON.
        """

        self.json_path = Path(json_path)

        if not self.json_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.json_path}")

        with open(self.json_path, "r") as f:
            self.config = json.load(f)

        # Extract everything we need
        self.cosima_file = self.config["cosima_file"]
        self.geometry_file = self.config["geometry_file"]
        self.revan_cfg = self.config["revan_output"]
        self.mimrec_cfg = self.config["mimrec_output"]
        self.energy_cut = self.config.get("energy_cut", [10,2000])
        self.max_events = self.config.get("max_events", 100000)

        # directory of THIS run
        self.dir = self.json_path.parent

        # where ALL output files go (per run)
        self.output_dir = self.dir / "results"
        os.makedirs(self.output_dir, exist_ok=True)

        # keep static copies inside /home/linusb/algo
        algo_dir = Path("/home/linusb/algo")
        for fpath in [self.cosima_file, self.geometry_file, self.revan_cfg, self.mimrec_cfg]:
            fpath = Path(fpath)
            if fpath.exists():
                dest = algo_dir / fpath.name
                if not dest.exists():
                    os.system(f"cp {fpath} {dest}")


    # --------------------------------------------------------------
    # Utility for running commands
    # --------------------------------------------------------------

    def run_command(self, cmd):
        print(f">>> {' '.join(cmd)}\n")
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        for line in process.stdout:
            print(line, end="")

        process.wait()

        if process.returncode != 0:
            raise RuntimeError(f"Command failed with exit code {process.returncode}")


    # --------------------------------------------------------------
    # Step 1: cosima
    # --------------------------------------------------------------

    def run_simulation(self):
        print("==== Step 1: Simulation ====")
        self.run_command(["cosima", self.cosima_file])


    # --------------------------------------------------------------
    # Step 2: revan
    # --------------------------------------------------------------

    def run_reconstruction(self):
        print("==== Step 2: Reconstruction ====")

        base = os.path.splitext(self.cosima_file)[0]
        sim_file = base + ".inc1.id1.sim.gz"

        if not os.path.exists(sim_file):
            print(f"Warning: Simulation file '{sim_file}' not found after cosima step.")

        self.run_command([
            "revan",
            "-c", self.revan_cfg,
            "-g", self.geometry_file,
            "-f", sim_file,
            "-a", "-n"
        ])


    # --------------------------------------------------------------
    # Step 3: mimrec
    # --------------------------------------------------------------

    def run_spectrum(self):
        print("==== Step 3: Spectrum Creation ====")

        base = os.path.splitext(self.cosima_file)[0]
        tra_file = base + ".inc1.id1.tra.gz"

        if not os.path.exists(tra_file):
            print(f"Warning: Tracking file '{tra_file}' not found after revan step.")

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
    # NEW: extract energies from .tra.gz so Comparator can use them
    # --------------------------------------------------------------

    def extract_energy(self):
        """
        Reads the .tra.gz file and pulls out all energy values.
        Writes them into reference_energy.txt or test_energy.txt.
        """

        base = os.path.splitext(self.cosima_file)[0]
        tra_file = base + ".inc1.id1.tra.gz"

        # Detect which run this is
        run_name = "reference_energy.txt" if "run_ref" in str(self.dir) else "test_energy.txt"
        out_file = self.output_dir / run_name

        energies = []

        # Use zcat to avoid python gzip overhead
        with os.popen(f"zcat {tra_file}") as f:
            for line in f:
                if "E=" in line or "E =" in line:
                    try:
                        val = float(line.replace("keV","").split("=")[1])
                        energies.append(val)
                    except:
                        pass

        with open(out_file, "w") as f:
            for e in energies:
                f.write(f"{e}\n")

        print(f"Extracted {len(energies)} energies → {out_file}")


    # --------------------------------------------------------------
    # organize results, BUT avoid moving any .source/.cfg/.setup
    # --------------------------------------------------------------

    def organize_results(self):

        base_prefix = os.path.splitext(os.path.basename(self.cosima_file))[0]

        safe_ext = (".sim.gz", ".tra.gz", ".root", ".C", ".txt", ".dat")

        # Only move files that match the simulation prefix & safe extensions
        for file in os.listdir('.'):
            if file.startswith(base_prefix) and file.endswith(safe_ext):
                os.rename(file, self.output_dir / file)

    def extract_energy_list(self):
        spectrum_file = self.output_dir / "spectrum.C"
        if not spectrum_file.exists():
            print(f"[DataFlow] No spectrum.C found at {spectrum_file}")
            return

        # determine filename
        if "run_ref" in str(self.output_dir):
            out_file = self.output_dir / "reference_energy.txt"
        else:
            out_file = self.output_dir / "test_energy.txt"

        energies = []

        with open(spectrum_file, "r") as f:
            for line in f:
                line = line.strip()

                # type A: Energy[i] = 215.0;
                if "Energy" in line and "=" in line:
                    try:
                        val = float(line.split("=")[1].replace(";", "").strip())
                        energies.append(val)
                        continue
                    except:
                        pass

                # type B: hEnergy->Fill(215.0);
                if "Fill" in line and "(" in line and ")" in line:
                    try:
                        val = float(line.split("(")[1].split(")")[0])
                        energies.append(val)
                        continue
                    except:
                        pass

                # type C: push_back(215.0)
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

    # --------------------------------------------------------------
    # Full pipeline
    # --------------------------------------------------------------

    def run_full_pipeline(self):
        os.chdir(self.dir)

        self.run_simulation()
        self.run_reconstruction()
        self.run_spectrum()
        self.organize_results()
        self.extract_energy_list()