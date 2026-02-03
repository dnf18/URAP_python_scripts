from dataflow_class import DataFlow
from Steering import Steering
from comparator import Comparator
from Reporter import Reporter
from pathlib import Path
import argparse, os, json


class Supervisor:
    def __init__(self, dir: str):
        self.dir = Path(dir).resolve()

    def run(self):

        # Output directories

        run_ref = self.dir / "run_ref"
        run_test = self.dir / "run_test"

        run_ref.mkdir(parents=True, exist_ok=True)
        run_test.mkdir(parents=True, exist_ok=True)

        # Get steering config

        steering = Steering.user_input()

        ref_cfg  = run_ref  / "steering_config.json"
        test_cfg = run_test / "steering_config.json"

        steering.save(ref_cfg)
        steering.save(test_cfg)

        # Run DataFlow (reference)
        df_ref = DataFlow(str(ref_cfg))
        df_ref.run_full_pipeline()
        df_ref.generate_histogram()
		
        #Run Dataflow (test)
        df_test = DataFlow(str(test_cfg))
        df_test.run_full_pipeline()
        df_test.generate_histogram()

        # Histogram paths
     
        ref_hist = run_ref / "results" / "energy_hist.json"
        test_hist = run_test / "results" / "energy_hist.json"

        comparison_json = self.dir / "comparison_results.json"
        histogram_dir = self.dir / "histograms"
        histogram_dir.mkdir(exist_ok=True)


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
        
        # Spectrum paths
        
        histograms = {"Overlay Comparison": overlay_plot}
        
        ref_meta = run_ref / "results" / "spectrum_meta.json"
        test_meta = run_test / "results" / "spectrum_meta.json"

        with open(ref_meta) as f:
            histograms["Reference Spectrum"] = json.load(f)["spectrum_png"]

        with open(test_meta) as f:
            histograms["Test Spectrum"] = json.load(f)["spectrum_png"]

        # Reporter Output
   
        reporter = Reporter(config_json=str(ref_cfg))
        reporter.generate_pdf(
            results=results,
            histograms = histograms
        )

        print("\n[Supervisor] COMPLETED SUCCESSFULLY\n")
