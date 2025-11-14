from dataflow_class import DataFlow #DataFLow Class
from Steering import Steering #Steering Class
from Comparator import Comparator #comparator Class
from Reporter import Reporter #Reporting Class
import argparse
import os
from pathlib import Path




"""""
Note: We can also write it so the sim files are written into the 

The Supervisor Class essentially coordinates all the pipelines.

Then the Steering Class uses the Directory the Supervisor class calls upon to get the specific Geometries and maybe what tests you want to take

It makes the Dataflow run all the simulations

Then it makes the sim go into the Comparison 

Then it makes the Comparison and those results, then get passed along to the Reporter class

The Reporter Class writes the pdf and attatches all the graphs and comparisons and data analysis

"""""




class Supervisor:
    def __init__(self, dir: str):
        self.dir = Path(dir)

    def run(self) -> int:
        # create run directories
        run_ref = self.dir / "run_ref"
        run_test = self.dir / "run_test"
        run_ref.mkdir(parents=True, exist_ok=True)
        run_test.mkdir(parents=True, exist_ok=True)

        print(f"Created run directories:\n- {run_ref}\n- {run_test}\n")

        # prompt once
        print("RUN CONFIGURATION (applies to both runs)")
        steering = Steering.user_input()

        # now save configs into BOTH run folders
        json_ref = (run_ref / "steering_config.json").resolve()
        json_test = (run_test / "steering_config.json").resolve()

        steering.save(json_ref)
        steering.save(json_test)

        print("Saved steering to both run_ref and run_test.\n")

        # run reference pipeline
        print("RUN 1: Reference")
        DataFlow(json_ref).run_full_pipeline()

        # run test pipeline
        print("RUN 2: Test")
        DataFlow(json_test).run_full_pipeline()

        # compare results
        comparison_json = Comparator(self.dir).compare_energy_hist()
        Reporter(comparison_json).generate_pdf(
            results={"pass": True},
            histograms={}
        )

        print("Done.")



def main() -> int: 
    parser = argparse.ArgumentParser(prog='MEGAlib end2end DualRun', description='Runs two MEGAlib versions and compares outputs')
    parser.add_argument("path", nargs="?", default=".", help='Finding path towards Existing Directory')
    args = parser.parse_args()

    if not os.path.isdir(args.path):
        print(f"Error: {args.path} is not a directory")
    else:
        print(f"Directory exists and the program will continue : {os.path.abspath(args.path)}")
        Supervisor(args.path).run()


if __name__ == "__main__": #will exit the code if it doesnt work
    raise SystemExit(main())