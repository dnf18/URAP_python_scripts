#!/usr/bin/env python3
"""
patch_tgrs.py  <TGRS.cxx>  <TGRS_patched.cxx>

Patches the MEGAlib TGRS example so it can be used for the regression test:
  1. writes TGRS_output.root with the histograms TGRS, Sim, Difference and Component_i
     (Lindsey's output block, titles = the names given after -b)
  2. keeps Poisson errors: Sumw2() on every histogram and errors scaled together with the
     contents (without this ROOT reports sqrt(rate), which is wrong and breaks chi2)
  3. reads the measurement errors from TGRS.dat (5th column) and fixes the loop that wrote
     one element past the end of the Flux vector
  4. batch mode: no windows, no TApplication::Run() (the program hung there on Savio)

Every edit is an exact text replacement that must match exactly once, so the script works on the
old and the new MEGAlib version of TGRS.cxx and refuses to run on a file it does not understand.
"""
import sys

EDITS = [
 # includes
 ("#include <TLegend.h>\n",
  "#include <TLegend.h>\n#include <TFile.h>\n"),

 # measurement errors from the 5th column
 ("  vector<double> Flux;\n",
  "  vector<double> Flux;\n  vector<double> FluxErr;\n"),
 ("          Flux.push_back(T->GetTokenAtAsDouble(3));\n",
  "          Flux.push_back(T->GetTokenAtAsDouble(3));\n          FluxErr.push_back(T->GetTokenAtAsDouble(4));\n"),
 ("  for (int b = 1; b <= NBins; ++b) {\n    TGRSHist->SetBinContent(b, Flux[b-1]);\n  }\n",
  "  TGRSHist->Sumw2();\n  for (int b = 1; b <= NBins-1; ++b) {\n    TGRSHist->SetBinContent(b, Flux[b-1]);\n    TGRSHist->SetBinError(b, FluxErr[b-1]);\n  }\n"),

 # component histograms keep Poisson errors
 ("    Hists.push_back(new TH1D(m_BackgroundFileNames[i], m_BackgroundFileNames[i], NBinsSec-1, BinsSec));\n",
  "    Hists.push_back(new TH1D(m_BackgroundFileNames[i], m_BackgroundFileNames[i], NBinsSec-1, BinsSec));\n    Hists.back()->Sumw2();\n"),
 ("      Hists[i]->SetBinContent(b, Hists[i]->GetBinContent(b)/Hists[i]->GetBinWidth(b));\n",
  "      double BinErr = Hists[i]->GetBinError(b);\n"
  "      Hists[i]->SetBinContent(b, Hists[i]->GetBinContent(b)/Hists[i]->GetBinWidth(b));\n"
  "      Hists[i]->SetBinError(b, BinErr/Hists[i]->GetBinWidth(b));\n"),

 # combined simulation: errors added in quadrature
 ('  TH1D* SimHist = new TH1D("Sim", "Sim", NBinsSec-1, BinsSec);\n',
  '  TH1D* SimHist = new TH1D("Sim", "Sim", NBinsSec-1, BinsSec);\n  SimHist->Sumw2();\n'),
 ("      SimHist->SetBinContent(b, SimHist->GetBinContent(b) + Hists[i]->GetBinContent(b));\n",
  "      double SumErr = sqrt(SimHist->GetBinError(b)*SimHist->GetBinError(b) + Hists[i]->GetBinError(b)*Hists[i]->GetBinError(b));\n"
  "      SimHist->SetBinContent(b, SimHist->GetBinContent(b) + Hists[i]->GetBinContent(b));\n"
  "      SimHist->SetBinError(b, SumErr);\n"),

 # ROOT output (Lindsey's block)
 ("  delete Geometry;\n\n  return true;\n",
  '  TFile* OutFile = new TFile("TGRS_output.root", "RECREATE");\n'
  '  TGRSHist->Write();   // Measurements\n'
  '  SimHist->Write();    // Combined Simulations\n'
  '  DifferenceHist->Write();\n'
  '  for (unsigned int i = 0; i < Hists.size(); ++i) {\n'
  '    Hists[i]->SetName(Form("Component_%u", i));\n'
  '    Hists[i]->SetTitle(m_BackgroundTypeNames[i].Data());  // the name given after -b\n'
  '    Hists[i]->Write();\n'
  '  }\n'
  '  OutFile->Close();\n'
  '  cout<<"Wrote TGRS_output.root"<<endl;\n\n'
  '  delete Geometry;\n\n  return true;\n'),

 # batch mode
 ('  TApplication TGRSApp("TGRSApp", 0, 0);\n',
  '  gROOT->SetBatch(kTRUE);\n  TApplication TGRSApp("TGRSApp", 0, 0);\n'),
 ("  TGRSApp.Run();\n",
  "  // TGRSApp.Run();  // batch mode: no GUI loop (it never returns)\n"),
]

# Optional edits: applied only if the anchor is present (0 or 1 matches allowed).
# The old install's stock TGRS.cxx divides by the MTime object directly, which no longer
# compiles against its own MEGAlib library; the new install's file already uses GetAsDouble().
OPTIONAL_EDITS = [
 ("Hists[i]->Scale(1.0/SimFiles[i]->GetObservationTime());",
  "Hists[i]->Scale(1.0/SimFiles[i]->GetObservationTime().GetAsDouble());"),
]


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    src = open(sys.argv[1]).read()
    for old, new in EDITS:
        n = src.count(old)
        if n != 1:
            raise SystemExit(f"Edit anchor found {n} times (expected 1), file not patched:\n{old}")
        src = src.replace(old, new)
    applied = len(EDITS)
    for old, new in OPTIONAL_EDITS:
        n = src.count(old)
        if n > 1:
            raise SystemExit(f"Optional anchor found {n} times, file not patched:\n{old}")
        if n == 1:
            src = src.replace(old, new)
            applied += 1
            print("Applied optional edit: GetObservationTime() -> GetObservationTime().GetAsDouble()")
    open(sys.argv[2], "w").write(src)
    print(f"Patched file written to {sys.argv[2]} ({applied} edits)")


if __name__ == "__main__":
    main()
