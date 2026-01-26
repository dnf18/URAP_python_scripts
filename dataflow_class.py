import os
import json
import subprocess
from pathlib import Path
import numpy as np

class DataFlow:
    def __init__(self, json_path: str):
        self.json_path = Path(json_path)

        if not self.json_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.json_path}")

        # Load steering config
        with open(self.json_path, "r") as f:
            self.config = json.load(f)

        self.cosima_file = self.config["cosima_file"]
        self.geometry_file = self.config["geometry_file"]
        self.revan_cfg = self.config["revan_output"]
        self.mimrec_cfg = self.config["mimrec_output"]
        self.energy_cut = self.config.get("energy_cut", [10, 2000])
        self.max_events = self.config.get("max_events", 100000)

        # Directory for this run
        self.dir = self.json_path.parent

        # Output directory
        self.output_dir = self.dir / "results"
        self.output_dir.mkdir(exist_ok=True, parents=True)

    # -----------------------------
    # Run shell commands
    # -----------------------------
    def run_command(self, cmd):
        print(f">>> {' '.join(cmd)}\n")
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end="")
        process.wait()
        if process.returncode != 0:
            raise RuntimeError(f"Command failed: {cmd}")

    # -----------------------------
    # Step 1: Simulation
    # -----------------------------
    def run_simulation(self):
        print("\n==== Step 1: Simulation ====")
        self.run_command(["cosima", self.cosima_file])

    # -----------------------------
    # Step 2: Revan
    # -----------------------------
    def run_reconstruction(self):
        print("\n==== Step 2: Reconstruction ====")
        base = os.path.splitext(self.cosima_file)[0]
        sim_file = base + ".inc1.id1.sim.gz"
        if not os.path.exists(sim_file):
            print(f"Warning: Missing simulation file: {sim_file}")

        self.run_command([
            "revan",
            "-c", self.revan_cfg,
            "-g", self.geometry_file,
            "-f", sim_file,
            "-a", "-n"
        ])

    # -----------------------------
    # Step 3: Mimrec
    # -----------------------------
    def run_spectrum(self):
        print("\n==== Step 3: Spectrum Creation ====")
        base = os.path.splitext(self.cosima_file)[0]
        tra_file = base + ".inc1.id1.tra.gz"
        if not os.path.exists(tra_file):
            print(f"Warning: Missing tracking file: {tra_file}")

        output_spectrum = self.output_dir / "spectrum.C"
        self.run_command([
            "mimrec",
            "-c", self.mimrec_cfg,
            "-g", self.geometry_file,
            "-f", tra_file,
            "-s",
            "-o", str(output_spectrum)
        ])

    # -----------------------------
    # Extract energies
    # -----------------------------
    def extract_energy_list(self):
        spectrum_file = self.output_dir / "spectrum.C"

        run_type = "reference" if "run_ref" in str(self.dir) else "test"
        out_file = self.output_dir / f"{run_type}_energy.txt"

        if not spectrum_file.exists():
            print(f"[DataFlow] spectrum.C not found: {spectrum_file}")
            return

        energies = []
        with open(spectrum_file, "r") as f:
            for line in f:
                line = line.strip()
                if "Fill(" in line and ")" in line:
                    try:
                        val = float(line.split("(")[1].split(")")[0])
                        energies.append(val)
                    except:
                        pass

        with open(out_file, "w") as f:
            for e in energies:
                f.write(f"{e}\n")

        print(f"[DataFlow] Extracted {len(energies)} energies → {out_file}")

        # Create histogram JSON from extracted energies
        self.generate_histogram(energies)

    # -----------------------------
    # Generate histogram JSON
    # -----------------------------
    def generate_histogram(self, energies=None):
        if energies is None:
            # fallback to energy file
            run_type = "reference" if "run_ref" in str(self.dir) else "test"
            energy_txt = self.output_dir / f"{run_type}_energy.txt"
            if not energy_txt.exists():
                raise FileNotFoundError(f"Energy file not found: {energy_txt}")
            with open(energy_txt, "r") as f:
                energies = [float(line.strip()) for line in f if line.strip()]

        hist_json = self.output_dir / "energy_hist.json"
        counts, edges = np.histogram(energies, bins=50)
        hist_data = {"bins": counts.tolist(), "edges": edges.tolist()}

        with open(hist_json, "w") as f:
            json.dump(hist_data, f, indent=4)

        print(f"[DataFlow] Histogram saved → {hist_json}")

    # -----------------------------
    # Organize result files
    # -----------------------------
    def organize_results(self):
        base_prefix = os.path.splitext(os.path.basename(self.cosima_file))[0]
        safe_ext = (".sim.gz", ".tra.gz", ".root", ".C", ".txt", ".dat")

        for file in os.listdir('.'):
            if file.startswith(base_prefix) and file.endswith(safe_ext):
                os.rename(file, self.output_dir / file)

    # -----------------------------
    # Full pipeline
    # -----------------------------
    def run_full_pipeline(self):
        os.chdir(self.dir)
        self.run_simulation()
        self.run_reconstruction()
        self.run_spectrum()
        self.organize_results()
        self.extract_energy_list()
