#!/usr/bin/env python3
"""
Plot EGM sigma-variable distributions for a single sample in multiple pT ranges.

This script reads one ROOT file containing an `egHLTTree` and overlays the same
variable for three different pT (eg_et) regions on a single canvas.

pT ranges (GeV):
  - 30  to 400
  - 400 to 1200
  - 1200 to 3000

Usage:
  python3 plot_Phase2_EGM_Sigma_PtRanges.py input.root [output_base]

Example:
  python3 plot_Phase2_EGM_Sigma_PtRanges.py ../Ntuple.root SigmaPtBins
"""

import os
import sys
import ROOT
from ROOT import TFile, TCanvas, TH1F, TLegend, gStyle, gROOT, TLatex


PT_BINS_GEV = [
    (30.0, 400.0),
    (400.0, 1200.0),
    (1200.0, 3000.0),
]


def setup_root_style():
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    gStyle.SetTitleSize(0.04, "xyz")
    gStyle.SetLabelSize(0.03, "xyz")
    gStyle.SetTitleOffset(1.2, "y")
    gStyle.SetTitleOffset(1.1, "x")
    gStyle.SetPadLeftMargin(0.12)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadBottomMargin(0.12)


def pt_selection(pt_min, pt_max):
    # eg_et is in GeV in the existing plotting script
    return f"(eg_et > {pt_min}) && (eg_et <= {pt_max})"


def plot_sigma_pt_ranges(root_path, output_base="sigma_pt_ranges"):
    if not os.path.exists(root_path):
        raise FileNotFoundError(f"File not found: {root_path}")

    f = TFile(root_path, "READ")
    if f.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {root_path}")

    tree = f.Get("egHLTTree")
    if not tree:
        raise RuntimeError("TTree 'egHLTTree' not found in the file.")

    print(f"🌳 Tree: {tree.GetName()} with {tree.GetEntries()} entries")
    print(f"✂️  pT bins (eg_et): {', '.join([f'{a:g}-{b:g} GeV' for a,b in PT_BINS_GEV])}")

    # Sigma variables (base) + sqrt variants (where meaningful)
    # Note: signed covariance terms can be negative; no sqrt variants added for those.
    variables = {
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
        # sqrt variants (non-negative ranges)
        "sqrt_sigmaIEtaIEta": {"bins": 50, "xmin": 0.0, "xmax": 0.173, "title": "#sqrt{#sigma_{i#eta i#eta}}", "xlabel": "#sqrt{#sigma_{i#eta i#eta}}", "draw": "sqrt(eg_sigmaIEtaIEta)"},
        "sqrt_sigma2uu": {"bins": 50, "xmin": 0.0, "xmax": 3.16, "title": "#sigma_{uu}", "xlabel": "#sigma_{uu}", "draw": "sqrt(eg_sigma2uu)"},
        "sqrt_sigma2vv": {"bins": 50, "xmin": 0.0, "xmax": 1.73, "title": "#sigma_{vv}", "xlabel": "#sigma_{vv}", "draw": "sqrt(eg_sigma2vv)"},
        "sqrt_sigma2ww": {"bins": 50, "xmin": 0.0, "xmax": 14.14, "title": "#sigma_{ww}", "xlabel": "#sigma_{ww}", "draw": "sqrt(eg_sigma2ww)"},
        "sqrt_sigma2xx": {"bins": 50, "xmin": 0.0, "xmax": 3.16, "title": "#sigma_{xx}", "xlabel": "#sigma_{xx}", "draw": "sqrt(eg_sigma2xx)"},
        "sqrt_sigma2yy": {"bins": 50, "xmin": 0.0, "xmax": 3.16, "title": "#sigma_{yy}", "xlabel": "#sigma_{yy}", "draw": "sqrt(eg_sigma2yy)"},
        "sqrt_sigma2zz": {"bins": 50, "xmin": 0.0, "xmax": 10.0, "title": "#sigma_{zz}", "xlabel": "#sigma_{zz}", "draw": "sqrt(eg_sigma2zz)"},
    }

    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen + 2]
    labels = [f"{a:g} < p_{{T}} ≤ {b:g} GeV" for a, b in PT_BINS_GEV]

    for var_key, cfg in variables.items():
        draw_expr = cfg.get("draw", var_key)
        print(f"\n📊 Plotting {var_key} (draw: {draw_expr})")

        hists = []
        for idx, ((pt_min, pt_max), color) in enumerate(zip(PT_BINS_GEV, colors)):
            hname = f"hist_{var_key}_{idx}"
            hist = TH1F(hname, "", cfg["bins"], cfg["xmin"], cfg["xmax"])
            sel = pt_selection(pt_min, pt_max)
            tree.Draw(f"{draw_expr}>>{hname}", sel, "goff")

            mean = hist.GetMean() if hist.GetEntries() > 0 else float("nan")
            rms = hist.GetRMS() if hist.GetEntries() > 0 else float("nan")
            if hist.GetEntries() > 0:
                hist.Scale(1.0 / hist.GetEntries())

            hist.SetLineColor(color)
            hist.SetLineWidth(2)
            hist.SetFillStyle(0)
            hists.append(hist)

            print(
                f"   bin {idx+1}: {pt_min:g}-{pt_max:g} GeV -> {hist.GetEntries()} entries, mean = {mean:.6g}, RMS = {rms:.6g}"
            )

        # handle case where all are empty
        max_val = max((h.GetMaximum() for h in hists), default=0.0)
        if max_val <= 0:
            max_val = 1.0

        # For log-y plots, set a small positive minimum so ROOT has something to draw.
        # (With normalized histograms, empty bins are common.)
        for h in hists:
            h.SetMinimum(1e-6)

        def draw_common(canvas):
            canvas.SetGridx(True)
            canvas.SetGridy(True)

            hists[0].SetMaximum(max_val * 1.6)
            hists[0].GetXaxis().SetTitle(cfg["xlabel"])
            hists[0].GetYaxis().SetTitle("a.u.")

            for i, h in enumerate(hists):
                h.Draw("" if i == 0 else "same")

            leg = TLegend(0.58, 0.72, 0.93, 0.89)
            leg.SetBorderSize(1)
            leg.SetFillStyle(1001)
            leg.SetFillColor(0)
            leg.SetTextSize(0.040)
            for h, lab in zip(hists, labels):
                leg.AddEntry(h, lab, "l")
            leg.Draw()

            # CMS labels (match style used in the other script)
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
            tex_private.SetTextFont(42)
            tex_private.DrawLatexNDC(0.25, 0.96, "#it{Phase2 Simulation}")

            # Keep Python references alive until after SaveAs().
            # Otherwise, PyROOT objects like TLegend/TLatex can be garbage-collected
            # and disappear from the saved canvas.
            canvas._keepalive = [leg, tex, tex_cms, tex_private]

        # Linear-y plot
        c_lin = TCanvas(f"c_{var_key}_lin", cfg["title"], 800, 800)
        c_lin.SetLogy(False)
        draw_common(c_lin)
        out_lin = f"{output_base}_{var_key}_lin.png"
        c_lin.SaveAs(out_lin)
        c_lin.Close()
        print(f"   💾 Saved: {out_lin}")

        # Log-y plot
        c_log = TCanvas(f"c_{var_key}_log", cfg["title"], 800, 800)
        c_log.SetLogy(True)
        draw_common(c_log)
        out_log = f"{output_base}_{var_key}_log.png"
        c_log.SaveAs(out_log)
        c_log.Close()
        print(f"   💾 Saved: {out_log}")

        for h in hists:
            h.Delete()

    f.Close()
    print("\n🎉 Done. Plots written as PNG files.")


def main():
    if len(sys.argv) < 2 or "--help" in sys.argv or "-h" in sys.argv:
        print(__doc__)
        sys.exit(0)

    root_path = sys.argv[1]
    output_base = sys.argv[2] if len(sys.argv) > 2 else "sigma_pt_ranges"

    setup_root_style()
    plot_sigma_pt_ranges(root_path, output_base)


if __name__ == "__main__":
    main()

