#!/usr/bin/env python3

"""
Plot per-filter COLZ temperature maps:
  - sigma2uu (y) vs pT (x)
  - sigma2vv (y) vs pT (x)
  - sigma2ww (y) vs pT (x)

Expected input histograms:
  - h2_sigma2uu_vs_pt_<filterName>
  - h2_sigma2vv_vs_pt_<filterName>
  - h2_sigma2ww_vs_pt_<filterName>
"""

import argparse
import os
import re
from typing import Dict, Tuple

import ROOT


def setup_root_style() -> None:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTitleSize(0.04, "xyz")
    ROOT.gStyle.SetLabelSize(0.03, "xyz")
    ROOT.gStyle.SetTitleOffset(1.2, "y")
    ROOT.gStyle.SetTitleOffset(1.1, "x")
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.14)
    ROOT.gStyle.SetPadTopMargin(0.06)
    ROOT.gStyle.SetPadBottomMargin(0.12)


def draw_cms_labels(canvas: ROOT.TCanvas) -> None:
    tex = ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.045)
    tex.DrawLatexNDC(0.84, 0.96, "14 TeV")
    tex_cms = ROOT.TLatex()
    tex_cms.SetTextSize(0.058)
    tex_cms.SetTextFont(42)
    tex_cms.DrawLatexNDC(0.12, 0.96, "#bf{CMS}")
    tex_private = ROOT.TLatex()
    tex_private.SetTextSize(0.045)
    tex_private.SetTextFont(42)
    tex_private.DrawLatexNDC(0.25, 0.96, "#it{Phase2 Simulation}")
    canvas._keepalive = getattr(canvas, "_keepalive", []) + [tex, tex_cms, tex_private]


def save_canvas_variants(canvas: ROOT.TCanvas, out_base: str) -> None:
    canvas.SetLogx(True)
    canvas.SetLogz(True)
    canvas.SaveAs(f"{out_base}_log.png")
    canvas.SetLogx(False)
    canvas.SetLogz(False)
    canvas.SaveAs(f"{out_base}_lin.png")


def sanitize_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def collect_h2_triplets(
    root_file: ROOT.TFile,
) -> Dict[str, Tuple[ROOT.TH2, ROOT.TH2, ROOT.TH2]]:
    uu = {}
    vv = {}
    ww = {}
    for key in root_file.GetListOfKeys():
        name = key.GetName()
        obj = root_file.Get(name)
        if not obj or not obj.InheritsFrom("TH2"):
            continue
        if name.startswith("h2_sigma2uu_vs_pt_"):
            filt = name.replace("h2_sigma2uu_vs_pt_", "", 1)
            uu[filt] = obj
        elif name.startswith("h2_sigma2vv_vs_pt_"):
            filt = name.replace("h2_sigma2vv_vs_pt_", "", 1)
            vv[filt] = obj
        elif name.startswith("h2_sigma2ww_vs_pt_"):
            filt = name.replace("h2_sigma2ww_vs_pt_", "", 1)
            ww[filt] = obj

    triplets = {}
    for filt in sorted(set(uu.keys()) | set(vv.keys()) | set(ww.keys())):
        if filt in uu and filt in vv and filt in ww:
            triplets[filt] = (uu[filt], vv[filt], ww[filt])
    return triplets


def draw_maps(
    triplets: Dict[str, Tuple[ROOT.TH2, ROOT.TH2, ROOT.TH2]],
    out_dir: str,
    out_prefix: str,
) -> None:
    for filt, (h2_uu_in, h2_vv_in, h2_ww_in) in triplets.items():
        h2_uu = h2_uu_in.Clone(f"{h2_uu_in.GetName()}__clone")
        h2_uu.SetDirectory(0)
        c_uu = ROOT.TCanvas(f"c_uu_{sanitize_filename(filt)}", filt, 900, 800)
        c_uu.SetGridx(True)
        c_uu.SetGridy(True)
        h2_uu.GetXaxis().SetTitle("p_{T} [GeV]")
        h2_uu.GetYaxis().SetTitle("#sigma^{2}_{uu}")
        h2_uu.GetZaxis().SetTitle("Entries")
        h2_uu.Draw("colz")
        txt_uu = ROOT.TLatex()
        txt_uu.SetNDC(True)
        txt_uu.SetTextSize(0.03)
        txt_uu.DrawLatex(0.12, 0.90, f"Filter: {filt}")
        draw_cms_labels(c_uu)
        c_uu._keepalive = getattr(c_uu, "_keepalive", []) + [txt_uu]
        save_canvas_variants(
            c_uu,
            os.path.join(out_dir, f"{out_prefix}_{sanitize_filename(filt)}_sigma2uu_vs_pt"),
        )
        c_uu.Close()

        h2_vv = h2_vv_in.Clone(f"{h2_vv_in.GetName()}__clone")
        h2_vv.SetDirectory(0)
        c_vv = ROOT.TCanvas(f"c_vv_{sanitize_filename(filt)}", filt, 900, 800)
        c_vv.SetGridx(True)
        c_vv.SetGridy(True)
        h2_vv.GetXaxis().SetTitle("p_{T} [GeV]")
        h2_vv.GetYaxis().SetTitle("#sigma^{2}_{vv}")
        h2_vv.GetZaxis().SetTitle("Entries")
        h2_vv.Draw("colz")
        txt_vv = ROOT.TLatex()
        txt_vv.SetNDC(True)
        txt_vv.SetTextSize(0.03)
        txt_vv.DrawLatex(0.12, 0.90, f"Filter: {filt}")
        draw_cms_labels(c_vv)
        c_vv._keepalive = getattr(c_vv, "_keepalive", []) + [txt_vv]
        save_canvas_variants(
            c_vv,
            os.path.join(out_dir, f"{out_prefix}_{sanitize_filename(filt)}_sigma2vv_vs_pt"),
        )
        c_vv.Close()

        h2_ww = h2_ww_in.Clone(f"{h2_ww_in.GetName()}__clone")
        h2_ww.SetDirectory(0)
        c_ww = ROOT.TCanvas(f"c_ww_{sanitize_filename(filt)}", filt, 900, 800)
        c_ww.SetGridx(True)
        c_ww.SetGridy(True)
        h2_ww.GetXaxis().SetTitle("p_{T} [GeV]")
        h2_ww.GetYaxis().SetTitle("#sigma^{2}_{ww}")
        h2_ww.GetZaxis().SetTitle("Entries")
        h2_ww.Draw("colz")
        txt_ww = ROOT.TLatex()
        txt_ww.SetNDC(True)
        txt_ww.SetTextSize(0.03)
        txt_ww.DrawLatex(0.12, 0.90, f"Filter: {filt}")
        draw_cms_labels(c_ww)
        c_ww._keepalive = getattr(c_ww, "_keepalive", []) + [txt_ww]
        save_canvas_variants(
            c_ww,
            os.path.join(out_dir, f"{out_prefix}_{sanitize_filename(filt)}_sigma2ww_vs_pt"),
        )
        c_ww.Close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot per-filter COLZ maps of sigma2uu/vv/ww vs pT."
    )
    parser.add_argument("--input-file", required=True, help="Input ROOT file")
    parser.add_argument("--out-dir", default="SigmaPlots", help="Output directory")
    parser.add_argument("--out-prefix", default="SigmaPerFilter", help="Output prefix")
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        raise FileNotFoundError(f"Input ROOT file not found: {args.input_file}")

    os.makedirs(args.out_dir, exist_ok=True)
    setup_root_style()

    f_in = ROOT.TFile.Open(args.input_file, "READ")
    if not f_in or f_in.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {args.input_file}")

    triplets = collect_h2_triplets(f_in)
    if not triplets:
        f_in.Close()
        raise RuntimeError(
            "No TH2 triplets found. Expected h2_sigma2uu/vv/ww_vs_pt_<filter>."
        )

    draw_maps(triplets, args.out_dir, args.out_prefix)
    f_in.Close()
    print(f"Saved sigma2uu/vv/ww vs pT COLZ maps for {len(triplets)} filters in: {args.out_dir}")


if __name__ == "__main__":
    main()

