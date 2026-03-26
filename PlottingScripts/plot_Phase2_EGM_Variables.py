#!/usr/bin/env python3
"""
🎯 ROOT Histogram Comparison Script for Multiple EGM Variables 🎯

This script takes multiple ROOT files and plots various EGM variables
as normalized histograms with both linear and log scale y-axis.

Usage: python plot_Phase2_EGM_Variables.py file1.root file2.root [file3.root ...] [output_name] [legend1] [legend2] [legend3] ...

Example: python3 plot_Phase2_EGM_Variables.py ../Ntuple_16_0_0_pre4_QCD.root ../Ntuple_16_0_0_pre4_SingleE.root ../Ntuple_16_0_0_pre4_ZEE.root ../Ntuple_16_0_0_pre4_ZpEE.root Comparison "QCD" "SingleE" "ZEE" "ZpEE"
"""

import sys
import os
import ROOT
from ROOT import TFile, TCanvas, TH1F, TLegend, gStyle, gROOT, TPaveText, TLatex
import argparse

# Minimum E_T (GeV) for candidates to be included in the distributions (pt > 30 GeV)
PT_MIN_GEV = 30.0

def setup_root_style():
    """Setup ROOT plot style for beautiful histograms"""
    gROOT.SetBatch(True)  # Run in batch mode for saving
    gStyle.SetOptStat(1111)  # Show all statistics
    gStyle.SetPalette(55)    # Beautiful color palette
    gStyle.SetTitleSize(0.04, "xyz")
    gStyle.SetLabelSize(0.03, "xyz")
    gStyle.SetTitleOffset(1.2, "y")
    gStyle.SetTitleOffset(1.1, "x")
    gStyle.SetPadLeftMargin(0.12)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadBottomMargin(0.12)

def plot_egm_variables(file_paths, output_name="egm_variables_comparison", legends=None):
    """Plot multiple variables from multiple ROOT files, creating separate plots for each"""
    
    # Set up default legends if not provided
    if legends is None or len(legends) == 0:
        legends = [f"File {i+1}" for i in range(len(file_paths))]
    elif len(legends) < len(file_paths):
        # Extend legends if not enough provided
        legends.extend([f"File {i+1}" for i in range(len(legends), len(file_paths))])
    
    print(f"🎯 Starting variable comparison plots! 🎯")
    print(f"📁 Processing {len(file_paths)} files:")
    for i, (file_path, legend) in enumerate(zip(file_paths, legends)):
        print(f"   {i+1}. {file_path} (Legend: {legend})")
    print(f"📂 Output base name: {output_name}")
    print(f"✂️  Cut: eg_et > {PT_MIN_GEV} GeV (pt > {PT_MIN_GEV} GeV)")
    print(f"🎨 Color scheme: Blue, Red, Green, Magenta, Orange, Cyan, Yellow, Pink, Violet, Teal")
    print("=" * 60)
    
    # Define variables to plot with their binning
    variables = {
        "eg_hcalHForHoverE": {"bins": 25, "xmin": 0.0, "xmax": 1, "title": "HCal H/E", "xlabel": "H/E"},
        "eg_et": {"bins": 100, "xmin": 0.0, "xmax": 3000.0, "title": "E_{T}", "xlabel": "E_{T} [GeV]"},
        "eg_energy": {"bins": 25, "xmin": 0.0, "xmax": 5000.0, "title": "Energy", "xlabel": "Energy [GeV]"},
        "eg_rawEnergy": {"bins": 25, "xmin": 0.0, "xmax": 5000.0, "title": "Raw Energy", "xlabel": "Raw Energy [GeV]"},
        "eg_nrClus": {"bins": 10, "xmin": 0.0, "xmax": 10.0, "title": "Number of Clusters", "xlabel": "Number of Clusters"},
        "eg_hgcaliso_layerclus": {"bins": 12, "xmin": 0.0, "xmax": 12.0, "title": "HgCal Isolated Clusters", "xlabel": "Layer Number"},
        "eg_phi": {"bins": 32, "xmin": -3.2, "xmax": 3.2, "title": "#phi", "xlabel": "#phi [rad]"},
        "eg_eta": {"bins": 30, "xmin": -3.0, "xmax": 3.0, "title": "#eta", "xlabel": "#eta"},
        "eg_sigmaIEtaIEta": {"bins": 50, "xmin": 0.0, "xmax": 0.03, "title": "#sigma_{i#eta i#eta}", "xlabel": "#sigma_{i#eta i#eta}"},
        "eg_sigma2uu": {"bins": 50, "xmin": 0.0, "xmax": 10.0, "title": "#sigma_{uu}", "xlabel": "#sigma^{2}_{uu}"},
        "eg_sigma2vv": {"bins": 50, "xmin": 0.0, "xmax": 3.0, "title": "#sigma_{vv}", "xlabel": "#sigma^{2}_{vv}"},
        "eg_sigma2ww": {"bins": 50, "xmin": 0.0, "xmax": 200.0, "title": "#sigma_{ww}", "xlabel": "#sigma^{2}_{ww}"},
        "eg_sigma2xx": {"bins": 50, "xmin": 0.0, "xmax": 10.0, "title": "#sigma_{xx}", "xlabel": "#sigma^{2}_{xx}"},
        "eg_sigma2xy": {"bins": 50, "xmin": -5.0, "xmax": 5.0, "title": "#sigma_{xy}", "xlabel": "#sigma^{2}_{xy}"},
        "eg_sigma2yy": {"bins": 50, "xmin": 0.0, "xmax": 10.0, "title": "#sigma_{yy}", "xlabel": "#sigma^{2}_{yy}"},
        "eg_sigma2yz": {"bins": 50, "xmin": -20.0, "xmax": 20.0, "title": "#sigma_{yz}", "xlabel": "#sigma^{2}_{yz}"},
        "eg_sigma2zx": {"bins": 50, "xmin": -20.0, "xmax": 20.0, "title": "#sigma_{zx}", "xlabel": "#sigma^{2}_{zx}"},
        "eg_sigma2zz": {"bins": 50, "xmin": 0.0, "xmax": 100.0, "title": "#sigma_{zz}", "xlabel": "#sigma^{2}_{zz}"},
        # Square root of sigma variables (only non-negative ones; ROOT draw expression in "draw")
        "sigmauu": {"bins": 50, "xmin": 0.0, "xmax": 3.16, "title": "#sigma_{uu}", "xlabel": "#sigma_{uu}", "draw": "sqrt(eg_sigma2uu)"},
        "sigmavv": {"bins": 50, "xmin": 0.0, "xmax": 1.73, "title": "#sigma_{vv}", "xlabel": "#sigma_{vv}", "draw": "sqrt(eg_sigma2vv)"},
        "sigmaww": {"bins": 50, "xmin": 0.0, "xmax": 14.14, "title": "#sigma_{ww}", "xlabel": "#sigma_{ww}", "draw": "sqrt(eg_sigma2ww)"},
        "sigmaxx": {"bins": 50, "xmin": 0.0, "xmax": 3.16, "title": "#sigma_{xx}", "xlabel": "#sigma_{xx}", "draw": "sqrt(eg_sigma2xx)"},
        "sigmayy": {"bins": 50, "xmin": 0.0, "xmax": 3.16, "title": "#sigma_{yy}", "xlabel": "#sigma_{yy}", "draw": "sqrt(eg_sigma2yy)"},
        "sigmazz": {"bins": 50, "xmin": 0.0, "xmax": 10.0, "title": "#sigma_{zz}", "xlabel": "#sigma_{zz}", "draw": "sqrt(eg_sigma2zz)"},
        "eg_invEInvP": {"bins": 30, "xmin": 0.0, "xmax": 0.01, "title": "1/E - 1/p", "xlabel": "1/E - 1/p"},
        "eg_invESeedInvP": {"bins": 30, "xmin": 0.0, "xmax": 0.1, "title": "1/E_{seed} - 1/p", "xlabel": "1/E_{seed} - 1/p"},
        "eg_trkDEta": {"bins": 30, "xmin": 0.0, "xmax": 0.08, "title": "#Delta#eta_{Track}", "xlabel": "#Delta#eta_{Track}"},
        #"eg_trkDEtaSeed": {"bins": 30, "xmin": 0.0, "xmax": 0.08, "title": "#Delta#eta_{Track}^{Seed}", "xlabel": "#Delta#eta_{Track}^{Seed}"},
        "eg_ecalPFIsol_default": {"bins": 25, "xmin": 0.0, "xmax": 500.0, "title": "ECAL PF Isolation", "xlabel": "ECAL PF Isolation [GeV]"},
        "eg_hcalPFIsol_default": {"bins": 25, "xmin": 0.0, "xmax": 100.0, "title": "HCAL PF Isolation", "xlabel": "HCAL PF Isolation [GeV]"},
        "eg_hgcalPFIsol_default": {"bins": 25, "xmin": -1.0, "xmax": 1.0, "title": "HgCal PF Isolation", "xlabel": "HgCal PF Isolation"},
        "eg_trkIsolV0_default": {"bins": 50, "xmin": 0.0, "xmax": 0.5, "title": "Track Isolation V0", "xlabel": "Track Isolation V0"},
        "eg_trkChi2_default": {"bins": 50, "xmin": 0.0, "xmax": 1.0, "title": "Track #chi^{2}", "xlabel": "Track #chi^{2}"},
        "eg_pms2_default": {"bins": 20, "xmin": 0.0, "xmax": 0.8, "title": "PMS2", "xlabel": "PMS2"},
        "eg_l1TrkIsoCMSSW": {"bins": 20, "xmin": 0.0, "xmax": 200.0, "title": "L1 Track Isolation CMSSW", "xlabel": "L1 Track Isolation CMSSW"}
    }
    
    # Open ROOT files
    try:
        files = []
        trees = []
        
        for i, file_path in enumerate(file_paths):
            file_obj = TFile(file_path, "READ")
            if file_obj.IsZombie():
                print(f"💀 Oops! File {i+1} ({file_path}) is a zombie! Check your file paths!")
                return
            files.append(file_obj)
            
            # Get the egHLTTree from each file
            tree = file_obj.Get("egHLTTree")
            if not tree:
                print(f"🌳 'egHLTTree' not found in file {i+1} ({file_path})!")
                return
            trees.append(tree)
            print(f"🌳 Tree {i+1}: {tree.GetName()} with {tree.GetEntries()} entries")
        
        # Define colors for different files (cycling through a nice palette)
        colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta+2, ROOT.kOrange+2, 
                 ROOT.kCyan+2, ROOT.kYellow+2, ROOT.kPink+2, ROOT.kViolet+2, ROOT.kTeal+2]
        
        # Loop over each variable and create separate plots
        for var_name, var_config in variables.items():
            print(f"\n📊 Creating plot for: {var_name}")
            print(f"   Bins: {var_config['bins']}, Range: {var_config['xmin']} to {var_config['xmax']}")
            
            # Create histograms for this variable for all files
            histograms = []
            # Selection: only candidates with pt (eg_et) > PT_MIN_GEV
            selection = f"eg_et > {PT_MIN_GEV}"
            draw_expr = var_config.get("draw", var_name)
            # For sigma variables, remove the spike at 0 (barrel-region candidates).
            sigma_zero_var = None
            var_lower = var_name.lower()
            if "sigma" in var_lower:
                if var_name.startswith("eg_sigma"):
                    # Underlying ROOT quantity is the histogram variable itself.
                    sigma_zero_var = var_name
                elif draw_expr.startswith("sqrt(") and draw_expr.endswith(")"):
                    # Underlying ROOT quantity comes from the sqrt(...) draw expression.
                    # Example: draw_expr = "sqrt(eg_sigma2uu)" -> sigma_zero_var = "eg_sigma2uu"
                    inner = draw_expr[len("sqrt("):-1]
                    if inner.startswith("eg_sigma"):
                        sigma_zero_var = inner
            if sigma_zero_var is not None:
                selection = f"{selection} && {sigma_zero_var} != 0"
            for i, tree in enumerate(trees):
                hist = TH1F(f"hist_{i}_{var_name}", "", var_config['bins'], var_config['xmin'], var_config['xmax'])
                tree.Draw(f"{draw_expr}>>hist_{i}_{var_name}", selection, "goff")
                print(f"   File {i+1}: {hist.GetEntries()} entries")
                
                # Normalize histogram to 1
                if hist.GetEntries() > 0:
                    hist.Scale(1.0 / hist.GetEntries())
                
                # Style the histogram
                hist.SetLineColor(colors[i % len(colors)])
                hist.SetLineWidth(2)
                hist.SetFillColor(0)
                hist.SetFillStyle(0)
                hist.SetStats(0)  # Hide statistics box
                
                histograms.append(hist)
            
            # Set up axis labels for the first histogram (they'll be the same for all)
            if histograms:
                histograms[0].GetXaxis().SetTitle(var_config['xlabel'])
                histograms[0].GetXaxis().SetTitleSize(0.05)
                histograms[0].GetXaxis().SetLabelSize(0.04)
                histograms[0].GetXaxis().SetTitleOffset(1.1)
                histograms[0].GetYaxis().SetTitle("a.u.")
                histograms[0].GetYaxis().SetLabelSize(0.045)
                histograms[0].GetYaxis().SetTitleSize(0.06)
                histograms[0].GetYaxis().SetTitleOffset(0.8)
            
            # Create canvas and draw - simplified without ratio panel
            # Draw both linear and log scale versions
            for ytag, logy in (("lin", False), ("log", True)):
                canvas = TCanvas(f"canvas_{var_name}_{ytag}", f"{var_config['title']} Comparison ({ytag})", 800, 800)
                canvas.SetLogy(logy)
                canvas.SetGridx(True)
                canvas.SetGridy(True)

                # Find the maximum to set proper scale
                max_val = max([hist.GetMaximum() for hist in histograms])
                histograms[0].SetMaximum(max_val * 1.5)

                # For log scale, avoid issues with bins at/near 0 by setting a sensible minimum.
                if logy:
                    min_positive = None
                    for hist in histograms:
                        nbins = hist.GetNbinsX()
                        for b in range(1, nbins + 1):
                            c = hist.GetBinContent(b)
                            if c > 0 and (min_positive is None or c < min_positive):
                                min_positive = c
                    if min_positive is not None:
                        histograms[0].SetMinimum(min_positive * 0.2)

                # Draw histograms (step-like line representation)
                for i, hist in enumerate(histograms):
                    if i == 0:
                        hist.Draw("HIST")
                    else:
                        hist.Draw("HIST SAME")

                # Add legend with custom labels
                legend = TLegend(0.62, 0.75, 0.93, 0.89)
                legend.SetBorderSize(1)
                legend.SetLineColor(1)
                legend.SetLineStyle(1)
                legend.SetLineWidth(1)
                legend.SetFillColor(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.045)

                # Add entries for all histograms
                for i, (hist, legend_label) in enumerate(zip(histograms, legends)):
                    legend.AddEntry(hist, legend_label, "l")
                legend.Draw()

                # Add CMS labels (keep existing text content)
                tex = TLatex()
                tex.SetTextFont(42)
                tex.SetTextSize(0.045)
                tex.SetLineWidth(2)
                tex.DrawLatexNDC(0.83, 0.96, "14 TeV")

                tex_cms = TLatex()
                tex_cms.SetTextSize(0.058)
                tex_cms.SetTextFont(42)
                tex_cms.DrawLatexNDC(0.12, 0.96, "#bf{CMS}")

                tex_private = TLatex()
                tex_private.SetTextSize(0.045)
                tex_private.SetTextFont(42)  # normal font
                tex_private.DrawLatexNDC(0.25, 0.96, "#it{Phase2 Simulation}")

                # Save the plot for this variable
                var_output_name = f"{output_name}_{var_name.replace('eg_', '')}"
                if ytag == "log":
                    # Backward-compatible output name for the log version
                    output_path = f"{var_output_name}.png"
                else:
                    output_path = f"{var_output_name}_lin.png"
                canvas.SaveAs(output_path)
                print(f"   💾 Saved: {output_path}")

                canvas.Close()

            # Clean up histograms for this variable
            for hist in histograms:
                hist.Delete()
            
            print(f"   ✅ Completed plot for {var_name}")
        
        print("\n" + "=" * 60)
        print("🎉 All variable comparisons complete! 🎉")
        print(f"📁 Check out your separate plots in the current directory!")
        
    except Exception as e:
        print(f"💥 Oops! Something went wrong: {e}")
    
    finally:
        # Clean up
        if 'files' in locals():
            for file_obj in files:
                file_obj.Close()

def main():
    """Main function to parse arguments and run the comparison"""
    
    # Show help if no arguments or --help
    if len(sys.argv) == 1 or "--help" in sys.argv or "-h" in sys.argv:
        print("""
🎯 Plot multiple EGM variables comparison from multiple ROOT files

Usage:
  python plot_Phase2_EGM_Variables.py file1.root file2.root [file3.root ...] [output_name] [legend1] [legend2] ...

Examples:
  python plot_Phase2_EGM_Variables.py file1.root file2.root file3.root
  python plot_Phase2_EGM_Variables.py file1.root file2.root file3.root my_comparison
  python plot_Phase2_EGM_Variables.py file1.root file2.root file3.root my_comparison "Ref" "Target1" "Target2"

Note: 
  - First file will be plotted in BLUE
  - Other files will be plotted in different colors (RED, GREEN, etc.)
  - You can provide up to 10 files with different colors
        """)
        sys.exit(0)
    
    # Parse arguments manually to handle the mixed file/string arguments
    all_args = sys.argv[1:]  # Get all command line arguments except script name
    
    if len(all_args) < 2:
        print("💀 Error: At least 2 ROOT files are required!")
        sys.exit(1)
    
    # First, collect all .root files (they come first)
    files = []
    output_name = "egm_variables_comparison"
    legends = []
    
    # Find where files end and output_name/legends begin
    # Collect all consecutive .root arguments
    i = 0
    while i < len(all_args) and all_args[i].endswith('.root'):
        files.append(all_args[i])
        i += 1
    
    # Check if we have at least 2 files
    if len(files) < 2:
        print("💀 Error: At least 2 ROOT files are required!")
        sys.exit(1)
    
    # The next argument (if exists) is the output_name
    if i < len(all_args):
        output_name = all_args[i]
        i += 1
        # Everything after output_name are legends
        legends = all_args[i:]
    
    # Check if files exist
    for file_path in files:
        if not os.path.exists(file_path):
            print(f"💀 File not found: {file_path}")
            sys.exit(1)
    
    # Setup ROOT style
    setup_root_style()
    
    # Let's plot this! 🚀
    plot_egm_variables(files, output_name, legends)

if __name__ == "__main__":
    main()
