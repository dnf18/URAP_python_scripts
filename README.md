# URAP_python_scripts — `pipeline-only` branch

This branch holds only what is needed to run the TGRS validation pipeline: the scripts and inputs in
`tgrs_pipeline/` plus the three comparison/report modules it imports (`Comparator.py`,
`TGRS_Comparator.py`, `TGRS_Reporter.py`). **Do not merge this branch into `main`:** it deletes files that
`main` still needs.

## To Install Libraries from requirements.txt
Installing in python virtual env:

```
1. sudo apt install python3-pip
2. mkdir $HOME/pyenv
3, cd $HOME/pyenv
4. python3 -m venv myenv
5. source $HOME/pyenv/myenv/bin/activate
6. git clone this repo (URAP_python_scripts)
7. cd URAP_python_scripts  # (after git cloning repo)
8. pip3 install -r requirements.txt
```

On Savio there is no `sudo`: use `module load anaconda3` and install anything missing with
`pip install --user <package>`. The pipeline needs numpy, scipy, matplotlib, uproot, reportlab and pypdf.
Check that `requirements.txt` lists them.

---

# TGRS validation pipeline (MEGAlib regression test)

End-to-end test that runs the TGRS cosmic-ray activation benchmark (cosmic photons, prompt protons and
delayed decays) with two MEGAlib/Geant4 installs on Savio, compares the two simulated spectra with each
other and with the TGRS measurement, and writes a PDF validation report with a line-by-line comparison.

* **A (reference):** MEGAlib 10.2 with Geant4 10.2.3
* **B (test):** MEGAlib built against Geant4 11.4.3 ("dc6")

The scripts were written for the Savio cluster (account `fc_cosi`). Paths and the account name are specific to
that setup; edit them at the top of `tgrs_pipeline/scripts/submit_chain.sh` and `tgrs_pipeline/scripts/run_tgrs.sh`
before using the pipeline elsewhere.

## Layout

```
Comparator.py, TGRS_Comparator.py, TGRS_Reporter.py   modules imported by the report script
tgrs_pipeline/
  scripts/   simulation chain, TGRS run, report, line comparison, build records
  inputs/    old/ and dc6/ - the .source files, geometry setup and TGRS.dat for each install
```

Keep all files of `tgrs_pipeline/scripts/` in one folder: `submit_chain.sh` finds its helper scripts next to itself.

## 1. Run the simulation chain

```
cd tgrs_pipeline/scripts
bash submit_chain.sh <run_name> <install_dir> <inputs_dir> [total_seconds]
```

Example for the reference build, 100,000 s per component:

```
bash submit_chain.sh old_full /global/home/groups/fc_cosi/MEGAlib ../inputs/old 100000
```

The script submits, with Slurm dependencies: the photon array, a merge job, the Step 1 array (prompt protons and
isotope files), a merge job, Step 2 (one job that builds `CosmicProtons.activation.dat` from all Step 1 isotope
files), the Step 3 array (decays) and a final merge job. Each array task writes into its own folder
`<stage>_NNNN`; the merge jobs link the task outputs into the run folder under a global numbering and write a stub
`.sim.gz` file.

Run parameters can be changed with environment variables when calling the script: `P1_TIME`, `P1_CPUS`, `PHO_TIME`,
`PHO_CPUS`, `DEC_TIME`, `DEC_CPUS`, `TRIG_RATE`, `THROTTLE`, `WALL`. Note: `mcosima -t N` is the total number of runs,
and the simulated time applies per run.

If a stage fails or jobs are cancelled:

* `resume_chain.sh <run_name> <install_dir> [total]` resubmits the merge jobs, Step 2, Step 3 and the final merge for
  an existing run without rerunning the simulations.
* `rerun_step3.sh <run_name> <install_dir> [total]` reruns Step 3 only; earlier attempts are moved to `step3_attemptK/`.

Check each merge log. The merge scripts report skipped or truncated files, duplicated events and the total
observation time of the decay component, and stop with an error code when something is wrong.

## 2. Run the TGRS program

```
cd <folder of the runs>/old_full      # then repeat for dc6_full
MEGALIB_INSTALL=<install_dir> TGRS_BIN=<path to TGRS binary> bash <path to repo>/tgrs_pipeline/scripts/run_tgrs.sh
```

This writes `TGRS_output.root` (histograms `TGRS`, `Sim`, `Difference` and `Component_0..2`). Set `NO_ACTIVATION=1`
to leave out the activation component. See the header of `run_tgrs.sh` for details.

Both builds are analysed with the same TGRS binary, built against the Geant4 11.4.3 install (the TGRS binary of the
old install fails at library load). Build it with `build_tgrs.sh`; the TGRS example is patched to write the ROOT file
by `patch_tgrs.py`.

## 3. Validation report

Use the run names you gave `submit_chain.sh` in step 1; the examples use `old_full` (A) and `dc6_full` (B).
Replace the placeholders with the folders that contain your two run folders.

```
module load anaconda3
export PYTHONPATH=<path to repo>:$PYTHONPATH
RUN_A=<folder of the runs>/old_full
RUN_B=<folder of the runs>/dc6_full
python3 tgrs_pipeline/scripts/TGRS_ROOT_Report.py $RUN_A/TGRS_output.root $RUN_B/TGRS_output.root \
  --tgrs tgrs_pipeline/inputs/dc6/TGRS.dat --comparator-dir <path to repo> \
  --ref-label "MEGAlib 10.2 (Geant4 10.2.3)" --test-label "MEGAlib dc6 (Geant4 11.4.3)" \
  --out analysis/results --reports analysis/reports
```

Overall PASS/FAIL is decided by the A-versus-B comparison only, with a tolerance test: at least 95% of the merged bins
agree within 10% (or 3 sigma), and the integrated rate in each of three energy bands (below 100, 100-1000, above
1000 keV) agrees within 5% (or 3 sigma). The thresholds can be changed with `--bin-tol`, `--min-fraction` and
`--band-tol`; `--chi2-decides` uses the old chi2 p-value instead, which fails at this precision for any few-percent
difference. The comparisons with the measurements are information only. **The thresholds have not yet been approved
by the project lead.**

## 4. Line-by-line comparison

This uses `RUN_A` and `RUN_B` from step 3.

```
mkdir -p analysis/lines
python3 tgrs_pipeline/scripts/TGRS_Line_Compare.py $RUN_A/TGRS_output.root $RUN_B/TGRS_output.root \
  --lines tgrs_pipeline/scripts/tgrs_lines.csv --out analysis/lines --min-rate 0.02 \
  --ref-label "10.2" --test-label "dc6" --component 2 \
  --append-to analysis/reports/TGRS_Validation_Report.pdf
```

* `tgrs_lines.csv` is the line list of Weidenspointner et al. 2005 (ApJS 156, 69, Tables 5 and 6), 213 entries, parsed
  from the paper (first identification of each line only). Check a row against the paper before relying on it.
* For every line the script measures the net count rate (window of +-1.5 FWHM minus a sideband continuum) in the data
  and in both simulations. The assumed resolution (`--fwhm0`, `--fwhm1`) is approximate.
* `line_comparison.csv` holds all 213 lines (column `in_report` marks those in the PDF); the PDF table shows lines at or
  above `--min-rate` counts/s. Weak lines and lines next to a stronger line can give negative or unreliable rates.
* The combined report is written next to the original as `TGRS_Validation_Report_with_lines.pdf`.

## Build records

These document how the dc6 install and the TGRS binary were built on Savio (not needed to run the pipeline):
`build_megalib_dc6.sh`, `finish_megalib_dc6.sh`, `build_healpix_deps.sh`, `build_tgrs.sh`, `patch_tgrs.py`.

## 5. Validation report in the LaTeX layout (optional)

`TGRS_Latex_Reporter.py` builds the validation report in a LaTeX layout (table of contents, numbered sections, page
numbers, links) straight from the two `TGRS_output.root` files. It is independent of step 3: it does not use
`Comparator.py` or `TGRS_Reporter.py`, and needs only `TGRS_Line_Compare.py` and `tgrs_lines.csv` in the same folder.
The PASS/FAIL decision is the same A-versus-B tolerance test as in step 3.

Requirements: numpy, scipy, matplotlib, uproot, and a LaTeX installation (`pdflatex` with the packages booktabs,
longtable, pdflscape, hyperref). On Savio: `module load anaconda3 texlive`.

```
module load anaconda3 texlive
RUN_A=<folder of the runs>/old_full
RUN_B=<folder of the runs>/dc6_full
python3 tgrs_pipeline/scripts/TGRS_Latex_Reporter.py $RUN_A/TGRS_output.root $RUN_B/TGRS_output.root \
  --tgrs tgrs_pipeline/inputs/dc6/TGRS.dat --lines tgrs_pipeline/scripts/tgrs_lines.csv \
  --ref-label "MEGAlib 10.2 (Geant4 10.2.3)" --test-label "MEGAlib dc6 (Geant4 11.4.3)" \
  --out analysis/latex_report
```

The result is `analysis/latex_report/TGRS_Validation_Report.pdf`, next to the `.tex` file and a `figs/` folder. If
`pdflatex` is not found, the script still writes the `.tex` file and the figures; compile them elsewhere (for example
in Overleaf) or run `pdflatex TGRS_Validation_Report.tex` twice.

What the report contains:

1. Configuration summary, then spectra, goodness of fit and ratio plots of each build against the measurements
   (information only).
2. A versus B: combined plot, residuals, metrics and the tolerance test that decides PASS/FAIL.
3. Components (photons, prompt protons, activation), build against build, with integrated rates per energy band.
4. The energy ranges where A and B differ (|pull| > 3, more than 3% apart, integrated difference above 5 sigma), with
   the catalogue lines in each range and zoom plots of the 12 largest; and the ranges where each build differs from the
   data (both builds / only Geant4 A / only Geant4 B).
5. Lines of Tables 5 and 6 of the paper: summary, lines where A and B differ (with the same line fitted in the
   activation component alone), lines where each build differs from the data, and Appendix A with all lines in the
   layout of the paper.
6. Appendix B: group patterns of the line differences (median B/A per group with a bootstrap interval).

Each section starts with a short "How to read this" box.

Options:

* `--bin-tol`, `--min-fraction`, `--band-tol`: the tolerance-test thresholds (defaults 0.10, 0.95, 0.05). The thresholds
  have not yet been approved by the project lead.
* `--comp-labels` and `--activation-component`: the names of the simulation components and the index of the activation
  component used for the activation-only line columns (default 2).
* `--sim-merge N`: merge N simulation bins above 245 keV (default 1).
* `--fwhm0`, `--fwhm1`: the assumed line width FWHM = fwhm0 + fwhm1 * E in keV (approximate).
* `--no-pdf`: write only the `.tex` file and the figures.

Notes:

* Lines are measured as in step 4 (net rate in a window of +-1.5 FWHM minus a sideband continuum). A line is shown as an
  upper limit (net + 2 sigma, written `<`) if it is not detected at 2 sigma. Lines where the other fit methods differ
  (for example a fit of the measured line shape) will give somewhat different rates.
* The comparison covers 50 keV and above, using the data bins below 245 keV and the simulation bins above.
* The labels of the two builds should contain "Geant4 x.y"; the headings "Only Geant4 x.y differs" use it.
