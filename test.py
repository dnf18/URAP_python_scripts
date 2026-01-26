from pathlib import Path
from dataflow_class import DataFlow
from Steering import Steering
from comparator import Comparator
from Reporter import Reporter
from Supervisor import Supervisor

import argparse
import os
from pathlib import Path

#!/usr/bin/env python3

def main():
    parser = argparse.ArgumentParser(
        prog="MEGAlib end2end DualRun",
        description="Runs two MEGAlib pipelines (reference vs test) and compares outputs"
    )
    parser.add_argument(
        "path", nargs="?", default=".", help="Path to the base directory for runs"
    )
    args = parser.parse_args()

    base_dir = Path(args.path).resolve()

    if not base_dir.is_dir():
        print(f"Error: {base_dir} is not a directory")
        return

    print(f"Directory exists, running validation in: {base_dir}\n")

    # Run the Supervisor
    supervisor = Supervisor(base_dir)
    supervisor.run()


if __name__ == "__main__":
    main()
