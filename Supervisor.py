#!/usr/bin/env python3 

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
        self.dir = dir          #This is basically just the path to the test folder (i.e: smth like MaxObservingCrabs)
#self. is basically just saying that it belongs to this specific file
        

    def run(self) -> int:
        steering = Steering.user_input()
        json_path = Path(self.dir) / "steering_config.json"
        steering.save(json_path)
 #Information is basically just the list of commands thats your telling Data Flow to run
        #information then goes to .simulations which tells the instructions

        DataFlow(json_path).run_full_pipeline()
     #  runs the sims like cosima, revena, geomega
        
        
        comparison_json = Comparator(self.dir).compare_energy_hist()
        reporter = Reporter(config_json=comparison_json)  #load that fake JSON
        reporter.generate_pdf(results={"pass": True, "details": "Fake results from Comparator"})

         #run the comparisons, the more comparisons, the more functions that will get held 
         #.write() makes the pdf


#doesn't have to be a "pass" could just be a 1 or 0, we shall decide when we write more of the code
        
        
      #  print("Overall status: ", "PASS" if results.get("pass") else "FAIL") 
       # return 0 if results.get("pass") else 1



def main() -> int: 

    #creates an argument and basically just says if this Path isn't valid, don't run
    parser = argparse.ArgumentParser(prog='MEGAlib end2end Tests', description='Simulates, compares, and reports differences between 2 MEGAlib Versions')
    parser.add_argument("path", nargs="?", default=".", help='Finding path towards Existing Directory')
    args = parser.parse_args()

    if not os.path.isdir(args.path):
        print(f"Error: {args.path} is not a directory")
    else:
        print(f"Directory exists and the program will continue : {os.path.abspath(args.path)}")

        Supervisor(args.path).run()


if __name__ == "__main__": #will exit the code if it doesnt work
    raise SystemExit(main())




#parallelizaiton
#want to run multiple tests in parallel
#creates a supervisor class, which will either return a 1 or 0


#steering class just tells what files to run, and certain configurations, Just keeping the Information
#maybe just get rid of the DataFlow Class and have the specific runs in the Supervisor Class
#Use JSON Xfile
#check if the files exist and sanity check

