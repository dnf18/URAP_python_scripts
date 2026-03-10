import json
from pathlib import Path

"""""
This is the Steering Class, basically everything in this class is written on a JSON file

A JSON file is basically like a python dicionary that holds a string or in this case File Paths

We're saving the file paths to cosima, revan, mimrec, geometry, the energy

"""""

class Steering:

    def __init__(
        self,
        cosima_file: str,
        revan_output: str,
        mimrec_output: str,
        geometry_file: str,
        reference_source: str,
        test_source: str,
        energy_cut=(10, 2000),
        algorithm="Standard",
        max_events=100000
    ):

        self.config = {
            "tool": "revan",
            "cosima_file": str(Path(cosima_file).resolve()),
            "revan_output": str(Path(revan_output).resolve()),
            "mimrec_output": str(Path(mimrec_output).resolve()),
            "geometry_file": str(Path(geometry_file).resolve()),
            "energy_cut": list(energy_cut),
            "reconstruction_algorithm": algorithm,
            "max_events": max_events,
            "log_level": "INFO",
            "Megalib Main Branch Source": str(Path(reference_source).resolve()),
            "Megalib V2 Source": str(Path(test_source).resolve())
        }


    @classmethod
    def user_input(cls):

        cwd = Path.cwd()

        auto_source = next(cwd.glob("*.source"), None)
        auto_geo = next(cwd.glob("*.geo.setup"), None)
        auto_revan = next(cwd.glob("*.revan.cfg"), None)
        auto_mimrec = next(cwd.glob("*.mimrec.cfg"), None)

        if auto_source:
            print(f"Auto-detected source file: {auto_source}")
            cosima = str(auto_source)
        else:
            while True:
                cosima = input("Enter path to Cosima .source file: ").strip()
                if not cosima.endswith(".source"):
                    print("Invalid file type must end with '.source'")
                    continue
                if not Path(cosima).exists():
                    print("File not found at that path")
                    continue
                break

        if auto_geo:
            print(f"Auto-detected geometry file: {auto_geo}")
            geometry = str(auto_geo)
        else:
            while True:
                geometry = input("Enter path to Geometry .geo.setup file: ").strip()
                if not geometry.endswith(".geo.setup"):
                    print("Invalid file type must end with '.geo.setup'")
                    continue
                if not Path(geometry).exists():
                    print("File not found at that path")
                    continue
                break

        if auto_revan:
            print(f"Auto-detected Revan config: {auto_revan}")
            revan = str(auto_revan)
        else:
            while True:
                revan = input("Enter path to revan .revan.cfg file: ").strip()
                if not revan.endswith(".revan.cfg"):
                    print("Invalid file type must end with '.revan.cfg'")
                    continue
                if not Path(revan).exists():
                    print("File not found at that path")
                    continue
                break

        if auto_mimrec:
            print(f"Auto-detected Mimrec config: {auto_mimrec}")
            mimrec = str(auto_mimrec)
        else:
            while True:
                mimrec = input("Enter path to mimrec .mimrec.cfg file: ").strip()
                if not mimrec.endswith(".mimrec.cfg"):
                    print("Invalid file type must end with '.mimrec.cfg'")
                    continue
                if not Path(mimrec).exists():
                    print("File not found at that path")
                    continue
                break


        try:
            e_min = float(input("Enter minimum energy cut (keV) [default 10]: ") or 10)
            e_max = float(input("Enter maximum energy cut (keV) [default 2000]: ") or 2000)
        except ValueError:
            e_min, e_max = 10, 2000


        algo = input("Enter reconstruction algorithm [default Standard]: ").strip() or "Standard"

        max_ev = input("Enter max events [default 100000]: ").strip()
        max_ev = int(max_ev) if max_ev else 100000


        reference_source = input("Input your first path for MEGAlib: ").strip()
        test_source = input("Input your 2nd path for MEGAlib: ").strip()


        return cls(
            cosima_file=cosima,
            revan_output=revan,
            mimrec_output=mimrec,
            geometry_file=geometry,
            reference_source=reference_source,
            test_source=test_source,
            energy_cut=(e_min, e_max),
            algorithm=algo,
            max_events=max_ev
        )


    def save(self, path="./steering_config.json"):

        path = Path(path)
        run_dir = path.parent

        run_dir.mkdir(parents=True, exist_ok=True)

        src_cosima = Path(self.config["cosima_file"])
        src_geometry = Path(self.config["geometry_file"])
        src_revan = Path(self.config["revan_output"])
        src_mimrec = Path(self.config["mimrec_output"])

        dst_cosima = run_dir / src_cosima.name
        dst_geometry = run_dir / src_geometry.name
        dst_revan = run_dir / src_revan.name
        dst_mimrec = run_dir / src_mimrec.name


        # manual copy without shutil
        dst_cosima.write_bytes(src_cosima.read_bytes())
        dst_geometry.write_bytes(src_geometry.read_bytes())
        dst_revan.write_bytes(src_revan.read_bytes())
        dst_mimrec.write_bytes(src_mimrec.read_bytes())


        self.config["cosima_file"] = str(dst_cosima.resolve())
        self.config["geometry_file"] = str(dst_geometry.resolve())
        self.config["revan_output"] = str(dst_revan.resolve())
        self.config["mimrec_output"] = str(dst_mimrec.resolve())


        with open(path, "w") as f:
            json.dump(self.config, f, indent=4)


        print(f"Saved configuration to {path.resolve()}")


    def show(self):

        print("\n=== Current Steering Configuration ===")
        for k, v in self.config.items():
            print(f"  {k}: {v}")
        print()


if __name__ == "__main__":
    steering = Steering.user_input()
    steering.show()
    save_path = input("\nPath to save JSON config [./steering_config.json]: ").strip() or "./steering_config.json"
    steering.save(save_path)    #this is basically defining whats getting saved etc
