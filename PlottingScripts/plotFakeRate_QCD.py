#!/usr/bin/env python3

import os
import re
import argparse
import logging
import math
import ROOT

from PlotTDRStyle import ModTDRStyle
from PlotCMSLumi import CMS_lumi


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


PLOT_TYPES = ("pt", "eta", "phi")
X_TITLES = {
    "pt": "p_{T} [GeV]",
    "eta": "#eta",
    "phi": "#phi [rad]",
}


def collect_histogram_names(root_file):
    names = set()
    for key in root_file.GetListOfKeys():
        names.add(key.GetName())
    return names


def discover_sequences_and_filters(hist_names):
    """
    Discover structure from histogram naming:
      den: <sequence>_den_ele_<plotType>
      num: <sequence>_num_ele_<plotType>_<filter>
    """
    structure = {}
    num_pattern = re.compile(r"^(?P<seq>.+)_num_ele_(?P<ptype>pt|eta|phi)_(?P<flt>.+)$")

    for hname in hist_names:
        match = num_pattern.match(hname)
        if not match:
            continue
        seq = match.group("seq")
        ptype = match.group("ptype")
        flt = match.group("flt")
        structure.setdefault(seq, {}).setdefault(ptype, set()).add(flt)

    return structure


def build_efficiency_graph(root_file, numerator_name, denominator_name):
    h_num = root_file.Get(numerator_name)
    h_den = root_file.Get(denominator_name)

    if not h_num or not h_den:
        return None, "missing_hist", -1.0, -1.0
    num_int = h_num.Integral()
    den_int = h_den.Integral()
    if h_den.Integral() <= 0:
        return None, "zero_denominator", num_int, den_int
    if ROOT.TEfficiency.CheckConsistency(h_num, h_den):
        eff = ROOT.TEfficiency(h_num, h_den)
        graph = eff.CreateGraph()
        graph.SetName(f"eff_{numerator_name}")
        return graph, "teff", num_int, den_int

    # Fallback for inconsistent bins (e.g. num > den in some bins):
    # build a robust graph using clipped ratio per bin.
    logger.warning(
        "TEfficiency consistency failed, using fallback ratio graph: num=%s den=%s",
        numerator_name,
        denominator_name,
    )
    graph = ROOT.TGraphAsymmErrors()
    graph.SetName(f"effFallback_{numerator_name}")

    point_idx = 0
    for ibin in range(1, h_den.GetNbinsX() + 1):
        den = h_den.GetBinContent(ibin)
        if den <= 0:
            continue
        num = h_num.GetBinContent(ibin)
        num_clipped = min(max(num, 0.0), den)
        p = num_clipped / den
        err = math.sqrt(max(p * (1.0 - p) / den, 0.0))

        x = h_den.GetBinCenter(ibin)
        ex = 0.5 * h_den.GetBinWidth(ibin)
        graph.SetPoint(point_idx, x, p)
        graph.SetPointError(point_idx, ex, ex, err, err)
        point_idx += 1

    if point_idx == 0:
        return None, "empty_fallback", num_int, den_int
    return graph, "fallback", num_int, den_int


def make_single_plot(graph, sequence, filter_name, plot_type, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    canvas = ROOT.TCanvas("c", "c", 900, 800)
    canvas.cd()
    ROOT.gPad.SetTopMargin(0.09)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.03)

    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.0)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.SetLineWidth(2)
    graph.SetMinimum(0.0)
    graph.SetMaximum(1.2)
    graph.GetXaxis().SetTitle(X_TITLES[plot_type])
    graph.GetYaxis().SetTitle("Fake rate")
    graph.GetXaxis().SetTitleOffset(1.0)
    graph.GetYaxis().SetTitleOffset(1.2)
    graph.Draw("AP")

    if plot_type == "eta":
        graph.GetXaxis().SetLimits(-4.0, 4.0)
    elif plot_type == "phi":
        graph.GetXaxis().SetLimits(-3.2, 3.2)

    legend = ROOT.TLegend(0.42, 0.76, 0.93, 0.90)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.028)
    legend.AddEntry(graph, f"{sequence} : {filter_name}", "LP")
    legend.Draw()

    CMS_lumi("QCD multijet", canvas, 13, 10, "Phase2 Simulation")

    safe_sequence = sequence.replace("/", "_")
    safe_filter = filter_name.replace("/", "_")
    out_base = os.path.join(output_dir, f"fakeRate_{safe_sequence}_{safe_filter}_{plot_type}")
    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".pdf")


def main():
    parser = argparse.ArgumentParser(
        description="Plot fake rate vs pt/eta/phi for each filter from QCD fake-rate histograms."
    )
    parser.add_argument("--input-file", required=True, help="Input ROOT file from fake-rate analyzer")
    parser.add_argument(
        "--output-dir",
        default="FakeRate_Plots",
        help="Directory to store output plots",
    )
    parser.add_argument(
        "--sequence",
        default="",
        help="Optional regex to select only matching sequences",
    )
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ModTDRStyle()

    if not os.path.isfile(args.input_file):
        raise FileNotFoundError(f"Input file not found: {args.input_file}")

    in_file = ROOT.TFile.Open(args.input_file)
    if not in_file or in_file.IsZombie():
        raise RuntimeError(f"Unable to open input ROOT file: {args.input_file}")

    hist_names = collect_histogram_names(in_file)
    structure = discover_sequences_and_filters(hist_names)
    if not structure:
        raise RuntimeError("No matching fake-rate histograms found in input file.")

    seq_regex = re.compile(args.sequence) if args.sequence else None

    total_plots = 0
    for sequence, by_plot in structure.items():
        if seq_regex and not seq_regex.search(sequence):
            continue

        for plot_type in PLOT_TYPES:
            den = f"{sequence}_den_ele_{plot_type}"
            if den not in hist_names:
                logger.warning("Missing denominator: %s", den)
                continue

            filters = sorted(by_plot.get(plot_type, []))
            for filter_name in filters:
                num = f"{sequence}_num_ele_{plot_type}_{filter_name}"
                if num not in hist_names:
                    continue

                graph, mode, num_int, den_int = build_efficiency_graph(in_file, num, den)
                if not graph:
                    logger.warning(
                        "Skipping: mode=%s num=%s (int=%.3f) den=%s (int=%.3f)",
                        mode, num, num_int, den, den_int
                    )
                    continue

                out_dir = os.path.join(args.output_dir, sequence, plot_type)
                make_single_plot(graph, sequence, filter_name, plot_type, out_dir)
                total_plots += 1

    in_file.Close()
    logger.info("Done. Created %d plots in %s", total_plots, args.output_dir)


if __name__ == "__main__":
    main()
