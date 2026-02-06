import os
import json
import subprocess
from pathlib import Path
import numpy as np


"""""
DataFlow Class

This class is responsible for running the MEGAlib pipeline for ONE run directory.

It takes in a steering_config.json and then:

1) Runs cosima (simulation) -> produces *.sim.gz
2) Runs revan  (reconstruction) -> produces *.tra.gz
3) Runs mimrec -s to produce a ROOT macro: results/spectrum.C
4) Parses spectrum.C to extract Fill(E) values and bins them -> results/energy_hist.json (for Comparator)
5) Patches spectrum.C so it can SaveAs(...) and runs ROOT to make a PNG -> results/*_spectrum.png (for Reporter)
6) Writes a small spectrum_meta.json so Supervisor/Reporter can grab the PNG path easily

"""""


class DataFlow:
    def __init__(self, json_path: str):
        self.json_path = Path(json_path)

        if not self.json_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.json_path}")

        #steering config
        with open(self.json_path, "r") as f:
            self.config = json.load(f)

        self.cosima_file = self.config["cosima_file"]
        self.geometry_file = self.config["geometry_file"]
        self.revan_cfg = self.config["revan_output"]
        self.mimrec_cfg = self.config["mimrec_output"]

        #telling where JSON lives
        self.dir = self.json_path.parent

        #output directory
        self.output_dir = self.dir / "results"
        self.output_dir.mkdir(exist_ok=True, parents=True)

    #runs shell command
    def run_command(self, cmd, cwd=None):
        cwd = str(cwd) if cwd is not None else None
        print(f">>> {' '.join(map(str, cmd))}\n")

        process = subprocess.Popen(
            list(map(str, cmd)),
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        for line in process.stdout:
            print(line, end="")

        process.wait()

        if process.returncode != 0:
            raise RuntimeError(f"Command failed (code {process.returncode}): {cmd}")

    #cosima Simulation

    def run_simulation(self):
        print("\n==== Step 1: Simulation (cosima) ====")
        self.run_command(["cosima", self.cosima_file], cwd=self.dir)

    #revan

    def run_reconstruction(self):
        print("\n==== Step 2: Reconstruction (revan) ====")

        base = os.path.splitext(self.cosima_file)[0]
        sim_file = base + ".inc1.id1.sim.gz"

        sim_path = self.dir / sim_file
        if not sim_path.exists():
            print(f"Warning: Missing simulation file: {sim_path}")

        self.run_command(
            [
                "revan",
                "-c", self.revan_cfg,
                "-g", self.geometry_file,
                "-f", sim_file,
                "-a", "-n"
            ],
            cwd=self.dir
        )

    #Mimrec (makes spectrum.C)

    def run_spectrum_macro(self):
        print("\n==== Step 3: Spectrum macro (mimrec) ====")

        base = os.path.splitext(self.cosima_file)[0]
        tra_file = base + ".inc1.id1.tra.gz"

        tra_path = self.dir / tra_file
        if not tra_path.exists():
            print(f"Warning: Missing tracking file: {tra_path}")

        macro_out = self.output_dir / "spectrum.C"

        self.run_command(
            [
                "mimrec",
                "-c", self.mimrec_cfg,
                "-g", self.geometry_file,
                "-f", tra_file,
                "-s",
                "-o", str(macro_out)
            ],
            cwd=self.dir
        )

        return macro_out

    #Step 3.5: Extract energies + build histogram JSON (for Comparator)

    def extract_energy_list(self):
        spectrum_file = self.output_dir / "spectrum.C"

        run_type = "reference" if "run_ref" in str(self.dir) else "test"
        out_file = self.output_dir / f"{run_type}_energy.txt"

        if not spectrum_file.exists():
            print(f"[DataFlow] spectrum.C not found: {spectrum_file}")
            return

        #spectrum.C from mimrec is binned already (SetBinContent + vector edges) so we take the binned values and compare those


        text = spectrum_file.read_text()

        # 1) Extract edges from the vector definition
        start_key = "std::vector<Double_t>"
        vec_start = text.find(start_key)

        if vec_start == -1:
            print("[DataFlow] Could not find x-axis vector (std::vector<Double_t> ...) in spectrum.C")
            return

        brace_open = text.find("{", vec_start)
        brace_close = text.find("};", brace_open)

        if brace_open == -1 or brace_close == -1:
            print("[DataFlow] Could not parse x-axis vector braces in spectrum.C")
            return

        vec_body = text[brace_open + 1:brace_close]

        edges = []
        for token in vec_body.split(","):
            token = token.strip()
            if token:
                try:
                    edges.append(float(token))
                except:
                    pass

        if len(edges) < 2:
            print("[DataFlow] Parsed edges list is too small, something went wrong")
            return

        #extracting SetBinContent(i, value)
        bins_count = len(edges) - 1
        bins = [0.0] * bins_count

        for line in text.splitlines():
            line = line.strip()
            if "SetBinContent(" not in line:
                continue
            try:
                inside = line.split("SetBinContent(", 1)[1].split(")", 1)[0]
                idx_str, val_str = inside.split(",", 1)
                idx = int(idx_str.strip())
                val = float(val_str.strip())

                #ROOT bins are 1..nbins (0 and nbins+1 are under/overflow)
                if 1 <= idx <= bins_count:
                    bins[idx - 1] = val
            except:
                pass

        #Save "energy" list file just as a debug artifact
        with open(out_file, "w") as f:
            for b in bins:
                f.write(f"{b}\n")

        print(f"[DataFlow] Extracted {len(bins)} bin contents → {out_file}")

        #create a histogram JSON directly from bins + edges
        self.generate_histogram((bins, edges))

    #generate histogram JSON

    def generate_histogram(self, energies):
        hist_json = self.output_dir / "energy_hist.json"

        # energies is now (bins, edges) coming from spectrum.C
        bins, edges = energies

        hist_data = {"bins": list(bins), "edges": list(edges)}

        with open(hist_json, "w") as f:
            json.dump(hist_data, f, indent=4)

        print(f"[DataFlow] Histogram saved → {hist_json}")

    #step 4: patch spectrum.C so ROOT can save a PNG

    def _ensure_includes(self, text: str) -> str:
        includes = [
            "#include <TCanvas.h>",
            "#include <TStyle.h>",
            "#include <TColor.h>",
        ]

        if any(inc in text for inc in includes):
            return text

        lines = text.splitlines(True)
        insert_at = 0

        #if a __CLING__ pragma block exists, insert after it.
        for i, line in enumerate(lines):
            if line.strip() == "#endif":
                insert_at = i + 1
                break

        include_block = "".join(inc + "\n" for inc in includes) + "\n"
        lines.insert(insert_at, include_block)
        return "".join(lines)

    def patch_macro_for_png(self, macro_path: Path, png_path: Path) -> Path:
        if not macro_path.exists():
            raise FileNotFoundError(f"Macro not found: {macro_path}")

        text = macro_path.read_text()

        #fix illegal function names like void MaxObservingCrab.spectrum()
        #illegal syntax will make it so your not able to save the PNG or use ROOT
        text = text.replace("void MaxObservingCrab.spectrum()", "void MaxObservingCrab_spectrum()")

        #add includes if missing
        text = self._ensure_includes(text)

        #makes it easier to find
        save_line = f'  c1->SaveAs("{png_path.as_posix()}");'

        if "SaveAs(" not in text:
            idx = text.rfind("}")
            if idx == -1:
                raise ValueError(f"Could not find closing brace in macro: {macro_path}")
            text = text[:idx] + save_line + "\n" + text[idx:]

        patched_path = macro_path.parent / (macro_path.stem + "_png.C")
        patched_path.write_text(text)
        return patched_path

    #Step 5: Run ROOT to produce PNG

    def run_root_macro(self, macro_path: Path):
        print("\n==== Step 4: ROOT render (batch) ====")

        #find the function name inside the macro (first "void NAME(" line)
        func_name = None
        with open(macro_path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("void "):
                    func_name = line.split()[1].split("(")[0]
                    break

        if func_name is None:
            raise RuntimeError(f"Could not find macro function name in: {macro_path}")

        #run ROOT from the macro directory, load then call
        self.run_command(
            [
                "root", "-l", "-b", "-q",
                "-e", f".L {macro_path.name}",
                "-e", f"{func_name}()"
            ],
            cwd=macro_path.parent
        )



    def make_spectrum_png(self) -> Path:
        macro_path = self.output_dir / "spectrum.C"

        run_type = "reference" if "run_ref" in str(self.dir) else "test"
        png_path = self.output_dir / f"{run_type}_spectrum.png"

        patched = self.patch_macro_for_png(macro_path, png_path)
        self.run_root_macro(patched)

        if not png_path.exists():
            raise RuntimeError(
                f"ROOT macro ran but PNG not found at {png_path}.\n"
                f"Check the patched macro: {patched}"
            )

        print(f"[DataFlow] Spectrum PNG saved → {png_path}")
        return png_path

    # Full pipeline

    def run_full_pipeline(self):
        self.run_simulation()
        self.run_reconstruction()
        self.run_spectrum_macro()

        # Make energy_hist.json for Comparator
        self.extract_energy_list()

        # Make spectrum PNG for Reporter
        png_path = self.make_spectrum_png()

        # Save spectrum path into a small JSON for Supervisor/Reporter
        meta_json = self.output_dir / "spectrum_meta.json"
        with open(meta_json, "w") as f:
            json.dump({"spectrum_png": str(png_path)}, f, indent=4)

        return png_path
