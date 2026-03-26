#!/usr/bin/env python3
"""
Make "temperature" (colz) plots of sigma variables vs candidate pT and eta.

Implementation detail:
  - Uses TProfile2D so each (eta, pT) bin stores the MEAN of the sigma variable.
  - pT axis is binned in log-spaced bins (useful for 30–3000 GeV).
  - Drawn with colz and saved in both log and linear scale variants.

Usage:
  python3 plot_Phase2_EGM_Sigma_Colz_PtEta.py input.root [output_dir] [output_prefix]

Example:
  python3 plot_Phase2_EGM_Sigma_Colz_PtEta.py ../Ntuple.root Test_COLZ SigmaPtEta
"""

import os
import sys
from array import array

import ROOT
from ROOT import TFile, TCanvas, TLegend, gStyle, gROOT, TLatex, TProfile2D, TH2D, TF1


# Candidate selection
PT_MIN_GEV = 30.0
PT_MAX_GEV = 3000.0

# Axis binning
ETA_BINS = 60
ETA_MIN = -3.0
ETA_MAX = 3.0

PT_BINS = 60  # number of log-spaced bins in pT
SIGMA_BINS_1D = 50  # bins for sigma axis in sigma vs eta/pt 2D plots


def setup_root_style():
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    gStyle.SetTitleSize(0.04, "xyz")
    gStyle.SetLabelSize(0.03, "xyz")
    gStyle.SetTitleOffset(1.2, "y")
    gStyle.SetTitleOffset(1.1, "x")
    gStyle.SetPadLeftMargin(0.12)
    gStyle.SetPadRightMargin(0.14)  # room for z palette
    gStyle.SetPadTopMargin(0.06)
    gStyle.SetPadBottomMargin(0.12)


def log_edges(xmin, xmax, nbins):
    if xmin <= 0 or xmax <= 0:
        raise ValueError("log_edges requires xmin/xmax > 0")
    if xmax <= xmin:
        raise ValueError("log_edges requires xmax > xmin")
    if nbins < 1:
        raise ValueError("log_edges requires nbins >= 1")

    import math

    log_min = math.log10(xmin)
    log_max = math.log10(xmax)
    step = (log_max - log_min) / nbins
    edges = [10 ** (log_min + i * step) for i in range(nbins + 1)]
    return array("d", edges)


def selection():
    # pT branch in the existing scripts is eg_et (GeV)
    return f"(eg_et > {PT_MIN_GEV}) && (eg_et <= {PT_MAX_GEV})"


def _draw_cms_labels(canvas, prim=None):
    """Draw CMS-style labels on the current pad. Returns list of drawables to keep alive."""
    tex = TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.045)
    tex.SetLineWidth(2)
    tex.DrawLatexNDC(0.84, 0.96, "14 TeV")
    tex_cms = TLatex()
    tex_cms.SetTextSize(0.058)
    tex_cms.SetTextFont(42)
    tex_cms.DrawLatexNDC(0.12, 0.96, "#bf{CMS}")
    tex_private = TLatex()
    tex_private.SetTextSize(0.045)
    tex_private.SetTextFont(42)
    tex_private.DrawLatexNDC(0.25, 0.96, "#it{Phase2 Simulation}")
    tex_cut = TLatex()
    tex_cut.SetTextSize(0.040)
    tex_cut.SetTextFont(42)
    tex_cut.DrawLatexNDC(0.12, 0.90, f"p_{{T}} > {PT_MIN_GEV:g} GeV")
    keep = [tex, tex_cms, tex_private, tex_cut]
    if prim is not None:
        keep.insert(0, prim)
    existing = getattr(canvas, "_keepalive", [])
    canvas._keepalive = existing + keep


def _save_canvas_variants(canvas, out_base, logx=False, logy=False, logz=False):
    """Save both log and linear variants of the same canvas."""
    canvas.SetLogx(bool(logx))
    canvas.SetLogy(bool(logy))
    canvas.SetLogz(bool(logz))
    out_log = f"{out_base}_log.png"
    canvas.SaveAs(out_log)
    print(f"   💾 Saved: {out_log}")

    canvas.SetLogx(False)
    canvas.SetLogy(False)
    canvas.SetLogz(False)
    out_lin = f"{out_base}_lin.png"
    canvas.SaveAs(out_lin)
    print(f"   💾 Saved: {out_lin}")


def main():
    if len(sys.argv) < 2 or "--help" in sys.argv or "-h" in sys.argv:
        print(__doc__)
        sys.exit(0)

    in_root = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    out_prefix = sys.argv[3] if len(sys.argv) > 3 else "SigmaColzPtEta"

    if not os.path.exists(in_root):
        raise FileNotFoundError(f"File not found: {in_root}")
    os.makedirs(out_dir, exist_ok=True)

    setup_root_style()

    f = TFile(in_root, "READ")
    if f.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {in_root}")

    tree = f.Get("egHLTTree")
    if not tree:
        raise RuntimeError("TTree 'egHLTTree' not found in the file.")

    print(f"🌳 Tree: {tree.GetName()} with {tree.GetEntries()} entries")
    print(f"✂️  Selection: {selection()}")

    # Sigma variables (raw + sqrt variants). For 1D colz: sigma axis range (smin, smax).
    variables = {
        "eg_sigmaIEtaIEta": {"title": "#sigma_{i#eta i#eta}", "zlabel": "Mean #sigma_{i#eta i#eta}", "smin": 0.0, "smax": 0.03},
        "eg_sigma2uu": {"title": "#sigma^{2}_{uu}", "zlabel": "Mean #sigma^{2}_{uu}", "smin": 0.0, "smax": 10.0},
        "eg_sigma2vv": {"title": "#sigma^{2}_{vv}", "zlabel": "Mean #sigma^{2}_{vv}", "smin": 0.0, "smax": 3.0},
        "eg_sigma2ww": {"title": "#sigma^{2}_{ww}", "zlabel": "Mean #sigma^{2}_{ww}", "smin": 0.0, "smax": 200.0},
        "eg_sigma2xx": {"title": "#sigma^{2}_{xx}", "zlabel": "Mean #sigma^{2}_{xx}", "smin": 0.0, "smax": 10.0},
        "eg_sigma2xy": {"title": "#sigma^{2}_{xy}", "zlabel": "Mean #sigma^{2}_{xy}", "smin": -5.0, "smax": 5.0},
        "eg_sigma2yy": {"title": "#sigma^{2}_{yy}", "zlabel": "Mean #sigma^{2}_{yy}", "smin": 0.0, "smax": 10.0},
        "eg_sigma2yz": {"title": "#sigma^{2}_{yz}", "zlabel": "Mean #sigma^{2}_{yz}", "smin": -20.0, "smax": 20.0},
        "eg_sigma2zx": {"title": "#sigma^{2}_{zx}", "zlabel": "Mean #sigma^{2}_{zx}", "smin": -20.0, "smax": 20.0},
        "eg_sigma2zz": {"title": "#sigma^{2}_{zz}", "zlabel": "Mean #sigma^{2}_{zz}", "smin": 0.0, "smax": 100.0},
        "sigmauu": {"title": "#sigma_{uu}", "zlabel": "Mean #sigma_{uu}", "draw": "sqrt(eg_sigma2uu)", "smin": 0.0, "smax": 3.16},
        "sigmavv": {"title": "#sigma_{vv}", "zlabel": "Mean #sigma_{vv}", "draw": "sqrt(eg_sigma2vv)", "smin": 0.0, "smax": 1.73},
        "sigmaww": {"title": "#sigma_{ww}", "zlabel": "Mean #sigma_{ww}", "draw": "sqrt(eg_sigma2ww)", "smin": 0.0, "smax": 14.14},
        "sigmaxx": {"title": "#sigma_{xx}", "zlabel": "Mean #sigma_{xx}", "draw": "sqrt(eg_sigma2xx)", "smin": 0.0, "smax": 3.16},
        "sigmayy": {"title": "#sigma_{yy}", "zlabel": "Mean #sigma_{yy}", "draw": "sqrt(eg_sigma2yy)", "smin": 0.0, "smax": 3.16},
        "sigmazz": {"title": "#sigma_{zz}", "zlabel": "Mean #sigma_{zz}", "draw": "sqrt(eg_sigma2zz)", "smin": 0.0, "smax": 10.0},
    }

    pt_edges = log_edges(PT_MIN_GEV, PT_MAX_GEV, PT_BINS)

    for key, cfg in variables.items():
        draw_expr = cfg.get("draw", key)
        prof_name = f"p2_{key}"

        prof = TProfile2D(
            prof_name,
            "",
            ETA_BINS,
            ETA_MIN,
            ETA_MAX,
            PT_BINS,
            pt_edges,
        )

        # Fill: z : y : x
        tree.Draw(f"{draw_expr}:eg_et:eg_eta>>{prof_name}", selection(), "prof goff")

        c = TCanvas(f"c_{key}", key, 900, 800)
        c.SetGridx(True)
        c.SetGridy(True)

        prof.GetXaxis().SetTitle("#eta")
        prof.GetYaxis().SetTitle("p_{T} [GeV]")
        prof.GetZaxis().SetTitle(cfg.get("zlabel", "Mean value"))

        prof.Draw("colz")

        # CMS-style labels
        tex = TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.045)
        tex.SetLineWidth(2)
        tex.DrawLatexNDC(0.84, 0.96, "14 TeV")

        tex_cms = TLatex()
        tex_cms.SetTextSize(0.058)
        tex_cms.SetTextFont(42)
        tex_cms.DrawLatexNDC(0.12, 0.96, "#bf{CMS}")

        tex_private = TLatex()
        tex_private.SetTextSize(0.045)
        tex_private.SetTextFont(42)
        tex_private.DrawLatexNDC(0.25, 0.96, "#it{Phase2 Simulation}")

        tex_cut = TLatex()
        tex_cut.SetTextSize(0.040)
        tex_cut.SetTextFont(42)
        tex_cut.DrawLatexNDC(0.12, 0.90, f"p_{{T}} > {PT_MIN_GEV:g} GeV")

        # Keep references alive (PyROOT can GC these otherwise)
        c._keepalive = [prof, tex, tex_cms, tex_private, tex_cut]

        out_base = os.path.join(out_dir, f"{out_prefix}_{key}")
        _save_canvas_variants(c, out_base, logy=True, logz=True)
        c.Close()

        prof.Delete()

        # ---- Sigma vs eta (2D colz: x=eta, y=sigma) ----
        smin = cfg.get("smin", 0.0)
        smax = cfg.get("smax", 1.0)
        h2_eta_name = f"h2_eta_{key}"
        h2_eta = TH2D(h2_eta_name, "", ETA_BINS, ETA_MIN, ETA_MAX, SIGMA_BINS_1D, smin, smax)
        tree.Draw(f"{draw_expr}:eg_eta>>{h2_eta_name}", selection(), "goff")
        h2_eta.GetXaxis().SetTitle("#eta")
        h2_eta.GetYaxis().SetTitle(cfg.get("title", key))
        h2_eta.GetZaxis().SetTitle("Entries")

        c_eta = TCanvas(f"c_eta_{key}", f"{cfg['title']} vs #eta", 900, 700)
        c_eta.SetGridx(True)
        c_eta.SetGridy(True)
        h2_eta.Draw("colz")
        _draw_cms_labels(c_eta, h2_eta)
        out_eta_base = os.path.join(out_dir, f"{out_prefix}_{key}_vs_eta")
        _save_canvas_variants(c_eta, out_eta_base, logz=True)
        c_eta.Close()
        h2_eta.Delete()

        # ---- Sigma vs pT (2D colz: x=pT, y=sigma) ----
        h2_pt_name = f"h2_pt_{key}"
        h2_pt = TH2D(h2_pt_name, "", PT_BINS, pt_edges, SIGMA_BINS_1D, smin, smax)
        tree.Draw(f"{draw_expr}:eg_et>>{h2_pt_name}", selection(), "goff")
        h2_pt.GetXaxis().SetTitle("p_{T} [GeV]")
        h2_pt.GetYaxis().SetTitle(cfg.get("title", key))
        h2_pt.GetZaxis().SetTitle("Entries")

        c_pt = TCanvas(f"c_pt_{key}", f"{cfg['title']} vs p_{{T}}", 900, 700)
        c_pt.SetGridx(True)
        c_pt.SetGridy(True)
        h2_pt.Draw("colz")
        
        # Overlay requested reference parameterization for sigma2vv: 0.64 + 0.0008*pT
        if key == "eg_sigma2vv":
            ref_curve_1 = TF1(
                f"ref_curve1_{key}",
                "0.64 + 0.0008*x",
                PT_MIN_GEV,
                PT_MAX_GEV,
            )
            ref_curve_1.SetLineColor(ROOT.kRed + 1)
            ref_curve_1.SetLineStyle(2)
            ref_curve_1.SetLineWidth(3)
            ref_curve_1.Draw("same")
            
            ref_curve_2 = TF1(
                f"ref_curve2_{key}",
                "0.81 + 0.0002*x",
                PT_MIN_GEV,
                PT_MAX_GEV,
            )
            ref_curve_2.SetLineColor(ROOT.kAzure + 1)
            ref_curve_2.SetLineStyle(9)
            ref_curve_2.SetLineWidth(3)
            ref_curve_2.Draw("same")
            
            #ref_curve_3 = TF1(
            #    f"ref_curve3_{key}",
            #    "0.64 + 0.0002*x",
            #    PT_MIN_GEV,
            #    PT_MAX_GEV,
            #)
            #ref_curve_3.SetLineColor(ROOT.kGreen + 2)
            #ref_curve_3.SetLineStyle(7)
            #ref_curve_3.SetLineWidth(3)
            #ref_curve_3.Draw("same")
            
            leg = TLegend(0.46, 0.74, 0.88, 0.88)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetTextSize(0.030)
            leg.AddEntry(ref_curve_1, "0.64 + 0.0008 #times p_{T}", "l")
            leg.AddEntry(ref_curve_2, "0.81 + 0.0002 #times p_{T}", "l")
            #leg.AddEntry(ref_curve_3, "0.64 + 0.0002 #times p_{T}", "l")
            leg.Draw()
            #c_pt._keepalive = [h2_pt, ref_curve_1, ref_curve_2, ref_curve_3, leg]
            c_pt._keepalive = [h2_pt, ref_curve_1, ref_curve_2, leg]

        # Keep the existing reference overlays for sigmavv as well.
        if key == "sigmavv":
            ref_curve_1 = TF1(
                f"ref_curve1_{key}",
                "sqrt(0.64 + 0.0008*x)",
                PT_MIN_GEV,
                PT_MAX_GEV,
            )
            ref_curve_1.SetLineColor(ROOT.kRed + 1)
            ref_curve_1.SetLineStyle(2)
            ref_curve_1.SetLineWidth(3)
            ref_curve_1.Draw("same")
            
            ref_curve_2 = TF1(
                f"ref_curve2_{key}",
                "sqrt(0.81 + 0.0002*x)",
                PT_MIN_GEV,
                PT_MAX_GEV,
            )
            ref_curve_2.SetLineColor(ROOT.kAzure + 1)
            ref_curve_2.SetLineStyle(9)
            ref_curve_2.SetLineWidth(3)
            ref_curve_2.Draw("same")

            #ref_curve_3 = TF1(
            #    f"ref_curve3_{key}",
            #    "sqrt(0.64 + 0.0002*x)",
            #    PT_MIN_GEV,
            #    PT_MAX_GEV,
            #)
            #ref_curve_3.SetLineColor(ROOT.kGreen + 2)
            #ref_curve_3.SetLineStyle(7)
            #ref_curve_3.SetLineWidth(3)
            #ref_curve_3.Draw("same")
            
            leg = TLegend(0.44, 0.66, 0.88, 0.88)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetTextSize(0.030)
            leg.AddEntry(ref_curve_1, "#sqrt{0.64 + 0.0008 #times p_{T}}", "l")
            leg.AddEntry(ref_curve_2, "#sqrt{0.81 + 0.0002 #times p_{T}}", "l")
            #leg.AddEntry(ref_curve_3, "#sqrt{0.64 + 0.0002 #times p_{T}}", "l")
            leg.Draw()
            c_pt._keepalive = [h2_pt, ref_curve_1, ref_curve_2, leg]
            #c_pt._keepalive = [h2_pt, ref_curve_1, ref_curve_2, ref_curve_3,leg]

        # Overlay requested reference parameterization for sigma2ww: 64 + 0.04*pT and 42 + 0.02*pT
        if key == "eg_sigma2ww":
            ref_curve_1 = TF1(
                f"ref_curve1_{key}",
                "64 + 0.04*x",
                PT_MIN_GEV,
                PT_MAX_GEV,
            )
            ref_curve_1.SetLineColor(ROOT.kRed + 1)
            ref_curve_1.SetLineStyle(2)
            ref_curve_1.SetLineWidth(3)
            ref_curve_1.Draw("same")
            
            ref_curve_2 = TF1(
                f"ref_curve2_{key}",
                "81 + 0.02*x",
                PT_MIN_GEV,
                PT_MAX_GEV,
            )
            ref_curve_2.SetLineColor(ROOT.kAzure + 1)
            ref_curve_2.SetLineStyle(9)
            ref_curve_2.SetLineWidth(3)
            ref_curve_2.Draw("same")
            
            #ref_curve_3 = TF1(
            #    f"ref_curve3_{key}",
            #    "64 + 0.02*x",
            #    PT_MIN_GEV,
            #    PT_MAX_GEV,
            #)
            #ref_curve_3.SetLineColor(ROOT.kGreen + 2)
            #ref_curve_3.SetLineStyle(7)
            #ref_curve_3.SetLineWidth(3)
            #ref_curve_3.Draw("same")
            
            leg = TLegend(0.46, 0.74, 0.88, 0.88)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetTextSize(0.030)
            leg.AddEntry(ref_curve_1, "64 + 0.04 #times p_{T}", "l")
            leg.AddEntry(ref_curve_2, "81 + 0.02 #times p_{T}", "l")
            #leg.AddEntry(ref_curve_3, "64 + 0.02 #times p_{T}", "l")
            leg.Draw()
            #c_pt._keepalive = [h2_pt, ref_curve_1, ref_curve_2, ref_curve_3, leg]
            c_pt._keepalive = [h2_pt, ref_curve_1, ref_curve_2, leg]
        # Overlay requested reference parameterization for sigmaww using sqrt of sigma2ww references.
        if key == "sigmaww":
            ref_curve_1 = TF1(
                f"ref_curve1_{key}",
                "sqrt(64 + 0.04*x)",
                PT_MIN_GEV,
                PT_MAX_GEV,
            )
            ref_curve_1.SetLineColor(ROOT.kRed + 1)
            ref_curve_1.SetLineStyle(2)
            ref_curve_1.SetLineWidth(3)
            ref_curve_1.Draw("same")
            
            ref_curve_2 = TF1(
                f"ref_curve2_{key}",
                "sqrt(81 + 0.02*x)",
                PT_MIN_GEV,
                PT_MAX_GEV,
            )
            ref_curve_2.SetLineColor(ROOT.kAzure + 1)
            ref_curve_2.SetLineStyle(9)
            ref_curve_2.SetLineWidth(3)
            ref_curve_2.Draw("same")

            #ref_curve_3 = TF1(
            #    f"ref_curve3_{key}",
            #    "sqrt(64 + 0.02*x)",
            #    PT_MIN_GEV,
            #    PT_MAX_GEV,
            #)
            #ref_curve_3.SetLineColor(ROOT.kGreen + 2)
            #ref_curve_3.SetLineStyle(7)
            #ref_curve_3.SetLineWidth(3)
            #ref_curve_3.Draw("same")

            leg = TLegend(0.44, 0.66, 0.88, 0.88)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetTextSize(0.030)
            leg.AddEntry(ref_curve_1, "#sqrt{64 + 0.04 #times p_{T}}", "l")
            leg.AddEntry(ref_curve_2, "#sqrt{81 + 0.02 #times p_{T}}", "l")
            #leg.AddEntry(ref_curve_3, "#sqrt{64 + 0.02 #times p_{T}}", "l")
            leg.Draw()
            c_pt._keepalive = [h2_pt, ref_curve_1, ref_curve_2, leg]
            #c_pt._keepalive = [h2_pt, ref_curve_1, ref_curve_2, ref_curve_3, leg]
        
        _draw_cms_labels(c_pt, h2_pt)
        out_pt_base = os.path.join(out_dir, f"{out_prefix}_{key}_vs_pt")
        _save_canvas_variants(c_pt, out_pt_base, logx=True, logz=True)
        c_pt.Close()
        h2_pt.Delete()

    f.Close()
    print("\n🎉 Done. colz maps written (log + linear).")


if __name__ == "__main__":
    main()

