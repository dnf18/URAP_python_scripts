from dataflow_class import DataFlow
from Steering import Steering
from comparator import Comparator
from Reporter import Reporter
from pathlib import Path
import argparse, os, json


class Supervisor:
    def __init__(self, dir: str):
        self.dir = Path(dir)

    def run(self):

        # Output directories

        run_ref = self.dir / "run_ref"
        run_test = self.dir / "run_test"

        run_ref.mkdir(parents=True, exist_ok=True)
        run_test.mkdir(parents=True, exist_ok=True)

        # Get steering config

        steering = Steering.user_input()

        json_ref  = run_ref  / "steering_config.json"
        json_test = run_test / "steering_config.json"

        steering.save(json_ref)
        steering.save(json_test)

        # Run DataFlow

        DataFlow(str(json_ref)).run_full_pipeline()
        DataFlow(str(json_test)).run_full_pipeline()

        # Histogram paths
     
        ref_hist = run_ref / "results" / "energy_hist.json"
        test_hist = run_test / "results" / "energy_hist.json"

        comparison_json = self.dir / "comparison_results.json"
        histogram_dir = self.dir / "histograms"
        histogram_dir.mkdir(exist_ok=True)

        # Spectrum paths

        ref_meta = run_ref / "results" / "spectrum_meta.json"
        test_meta = run_test / "results" / "spectrum_meta.json"

        with open(ref_meta) as f:
            ref_png = json.load(f)["spectrum_png"]

        with open(test_meta) as f:
            test_png = json.load(f)["spectrum_png"]

        # Run Comparator
      
        comp = Comparator(
            ref_json=str(ref_hist),
            test_json=str(test_hist),
            output_json=str(comparison_json),
            histogram_output_dir=str(histogram_dir),
            sigma_threshold=3.0
        )

        results = comp.compare()
        overlay_plot = comp.plot_overlay()

        # Reporter Output
   
        reporter = Reporter(config_json=str(json_ref))
        reporter.generate_pdf(
            results=results,
            histograms={
                "Reference Spectrum": ref_png,
                "Test Spectrum": test_png,
                "Overlay Comparison": overlay_plot
            }
        )

        print("\n[Supervisor] COMPLETED SUCCESSFULLY\n")
