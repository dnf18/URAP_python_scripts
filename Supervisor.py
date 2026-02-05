from dataflow_class import DataFlow #DataFLow Class
from Steering import Steering #Steering Class
from comparator import Comparator #comparator Class
from Reporter import Reporter #Reporting Class
import argparse
import os
from pathlib import Path




"""""

The Supervisor Class essentially coordinates all the pipelines.

Then the Steering Class uses the Directory the Supervisor class calls upon to get the specific Geometries

It makes the Dataflow run all the simulations and create the Spectrum PNGs

Then it makes the sim go into the Comparison 

Then it makes the Comparison and those results, then get passed along to the Reporter class

The Reporter Class writes the pdf and attatches all the graphs and comparisons and data analysis

"""""




class Supervisor:
    def __init__(self, dir: str):
        self.dir = Path(dir)

    def run(self) -> int:
        #create run directories
        run_ref = self.dir / "run_ref"
        run_test = self.dir / "run_test"
        run_ref.mkdir(parents=True, exist_ok=True)
        run_test.mkdir(parents=True, exist_ok=True)

        print(f"Created run directories:\n- {run_ref}\n- {run_test}\n")

        #prompt once
        print("RUN CONFIGURATION (applies to both runs)")
        steering = Steering.user_input()

        #now save configs into BOTH run folders
        json_ref = (run_ref / "steering_config.json").resolve()
        json_test = (run_test / "steering_config.json").resolve()

        steering.save(json_ref)
        steering.save(json_test)

        print("Saved steering to both run_ref and run_test.\n")

        #run reference pipeline
        print("RUN 1: Reference")
        DataFlow(str(json_ref)).run_full_pipeline()

        #run test pipeline
        print("RUN 2: Test")
        DataFlow(str(json_test)).run_full_pipeline()

        #compare results
        ref_hist = run_ref / "results" / "energy_hist.json"
        test_hist = run_test / "results" / "energy_hist.json"

        comparison_json = self.dir / "comparison_results.json"
        histogram_dir = self.dir / "histograms"
        histogram_dir.mkdir(exist_ok=True)

        comp = Comparator(
            ref_json=str(ref_hist),
            test_json=str(test_hist),
            output_json=str(comparison_json),
            histogram_output_dir=str(histogram_dir),
            sigma_threshold=3.0
        )

        results = comp.compare()
        overlay_plot = comp.plot_overlay()

        #spectrum png paths (for Reporter)
        ref_meta = run_ref / "results" / "spectrum_meta.json"
        test_meta = run_test / "results" / "spectrum_meta.json"

        import json
        with open(ref_meta) as f:
            ref_png = json.load(f)["spectrum_png"]
        with open(test_meta) as f:
            test_png = json.load(f)["spectrum_png"]

        Reporter(config_json=str(json_ref)).generate_pdf(
            results=results,
            histograms={
                "Reference Spectrum": ref_png,
                "Test Spectrum": test_png,
                "Overlay Comparison": overlay_plot
            }
        )

        print("Done.")
        return 0



def main() -> int: 
    parser = argparse.ArgumentParser(prog='MEGAlib end2end DualRun', description='Runs two MEGAlib versions and compares outputs')
    parser.add_argument("path", nargs="?", default=".", help='Finding path towards Existing Directory')
    args = parser.parse_args()

    if not os.path.isdir(args.path):
        print(f"Error: {args.path} is not a directory")
        return 2
    else:
        print(f"Directory exists and the program will continue : {os.path.abspath(args.path)}")
        return Supervisor(args.path).run()


if __name__ == "__main__": #will exit the code if it doesnt work
    raise SystemExit(main())
