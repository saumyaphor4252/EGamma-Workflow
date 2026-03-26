#!/usr/bin/env python3

"""
Extract per-filter sigma2vv and sigma2ww distributions from Phase-2 HLT events.

Example:
  python3 extractSigma2vvSigma2ww_PerFilter.py \
      --input-file input.root \
      --output sigma_per_filter.root \
      --max-events 1000
"""

import argparse
import logging
import os
import math
from array import array
from typing import Any, Dict, List

from DataFormats.FWLite import Events, Handle
import ROOT


MAX_DR = 0.1
PT_MIN_GEV = 30.0
PT_MAX_GEV = 3000.0
PT_BINS = 60
ETA_BINS = 60
ETA_MIN = -3.0
ETA_MAX = 3.0

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def get_filters(cms_path: str) -> List[str]:
    """Extract module names containing 'Filter' from a cms.Sequence string."""
    filters = []
    seq_str = str(cms_path).replace("cms.Sequence(", "").replace(")", "")
    for module in seq_str.split("+"):
        module = module.strip()
        if "Filter" in module:
            module = (
                module.replace("process.", "")
                .replace(" ", "")
                .replace("(", "")
                .replace(")", "")
            )
            filters.append(module)
    return filters


def get_filter_index(trig_evt: Any, filter_name: str) -> int:
    """Return filter index in trigger::TriggerEvent or invalid sizeFilters()."""
    for index in range(trig_evt.sizeFilters()):
        label = str(trig_evt.filterLabel(index))
        clean_label = label.split("[")[0].strip('"')
        if clean_label == filter_name:
            return index
    return trig_evt.sizeFilters()


def match_trig_objs(eta: float, phi: float, trig_objs: List[Any], max_dr: float = MAX_DR) -> List[Any]:
    """Match trigger objects to eta/phi by deltaR."""
    max_dr2 = max_dr * max_dr
    return [
        obj
        for obj in trig_objs
        if ROOT.reco.deltaR2(eta, phi, obj.eta(), obj.phi()) < max_dr2
    ]


def make_filter_histograms(output_file: ROOT.TFile, filter_names: List[str]) -> Dict[str, ROOT.TH1D]:
    """Create per-filter sigma histograms and sigma-vs-pt 2D maps."""
    histograms: Dict[str, ROOT.TH1D] = {}
    log_min = math.log10(PT_MIN_GEV)
    log_max = math.log10(PT_MAX_GEV)
    step = (log_max - log_min) / PT_BINS
    pt_edges = array("d", [10 ** (log_min + i * step) for i in range(PT_BINS + 1)])

    for filt in filter_names:
        h_vv = ROOT.TH1D(
            f"h_sigma2vv_{filt}",
            f"sigma2vv for {filt};sigma2vv;Entries",
            200,
            -0.02,
            0.20,
        )
        h_ww = ROOT.TH1D(
            f"h_sigma2ww_{filt}",
            f"sigma2ww for {filt};sigma2ww;Entries",
            200,
            -0.02,
            0.20,
        )
        h_vv.SetDirectory(output_file)
        h_ww.SetDirectory(output_file)
        histograms[h_vv.GetName()] = h_vv
        histograms[h_ww.GetName()] = h_ww

        # Mean sigma maps vs (eta, pT), in the same style as colz temperature maps
        p2_vv = ROOT.TProfile2D(
            f"p2_sigma2vv_{filt}",
            f"sigma2vv vs eta and pT for {filt};#eta;p_{{T}} [GeV];Mean #sigma^{{2}}_{{vv}}",
            ETA_BINS,
            ETA_MIN,
            ETA_MAX,
            PT_BINS,
            pt_edges,
        )
        p2_ww = ROOT.TProfile2D(
            f"p2_sigma2ww_{filt}",
            f"sigma2ww vs eta and pT for {filt};#eta;p_{{T}} [GeV];Mean #sigma^{{2}}_{{ww}}",
            ETA_BINS,
            ETA_MIN,
            ETA_MAX,
            PT_BINS,
            pt_edges,
        )
        p2_vv.SetDirectory(output_file)
        p2_ww.SetDirectory(output_file)
        histograms[p2_vv.GetName()] = p2_vv
        histograms[p2_ww.GetName()] = p2_ww

        # Requested COLZ maps: sigma2uu / sigma2ww (y) vs pT (x)
        h2_uu_pt = ROOT.TH2D(
            f"h2_sigma2uu_vs_pt_{filt}",
            f"sigma2uu vs pT for {filt};p_{{T}} [GeV];#sigma^{{2}}_{{uu}}",
            PT_BINS,
            pt_edges,
            80,
            0.0,
            10.0,
        )
        h2_vv_pt = ROOT.TH2D(
            f"h2_sigma2vv_vs_pt_{filt}",
            f"sigma2vv vs pT for {filt};p_{{T}} [GeV];#sigma^{{2}}_{{vv}}",
            PT_BINS,
            pt_edges,
            80,
            0.0,
            3.0,
        )
        h2_ww_pt = ROOT.TH2D(
            f"h2_sigma2ww_vs_pt_{filt}",
            f"sigma2ww vs pT for {filt};p_{{T}} [GeV];#sigma^{{2}}_{{ww}}",
            PT_BINS,
            pt_edges,
            80,
            0.0,
            200.0,
        )
        h2_uu_pt.SetDirectory(output_file)
        h2_vv_pt.SetDirectory(output_file)
        h2_ww_pt.SetDirectory(output_file)
        histograms[h2_uu_pt.GetName()] = h2_uu_pt
        histograms[h2_vv_pt.GetName()] = h2_vv_pt
        histograms[h2_ww_pt.GetName()] = h2_ww_pt
    return histograms


def get_sigma_vars_for_filter(eg_obj: Any, filter_name: str) -> tuple:
    """
    Read sigma2vv/sigma2ww from EgammaObject var-map for a given filter family.
    Returns (sigma2vv, sigma2ww).
    """
    if "L1Seeded" in filter_name:
        uu_name = "hltEgammaHGCALIDVarsL1Seeded_sigma2uu"
        vv_name = "hltEgammaHGCALIDVarsL1Seeded_sigma2vv"
        ww_name = "hltEgammaHGCALIDVarsL1Seeded_sigma2ww"
    else:
        uu_name = "hltEgammaHGCALIDVarsUnseeded_sigma2uu"
        vv_name = "hltEgammaHGCALIDVarsUnseeded_sigma2vv"
        ww_name = "hltEgammaHGCALIDVarsUnseeded_sigma2ww"

    sigma2uu = eg_obj.var(uu_name, 0)
    sigma2vv = eg_obj.var(vv_name, 0)
    sigma2ww = eg_obj.var(ww_name, 0)
    return sigma2uu, sigma2vv, sigma2ww


def process_events(
    events: Events,
    output_file: ROOT.TFile,
    all_filters: List[str],
    max_events: int,
    max_dr: float,
) -> None:
    """Main event loop: collect sigma2vv/sigma2ww per filter."""
    eg_unseeded_handle = Handle("std::vector<trigger::EgammaObject>")
    eg_l1seeded_handle = Handle("std::vector<trigger::EgammaObject>")
    trig_evt_handle = Handle("trigger::TriggerEvent")

    histograms = make_filter_histograms(output_file, all_filters)

    total_events = events.size() if max_events < 0 else min(events.size(), max_events)
    logger.info("Processing %d events", total_events)

    for i, event in enumerate(events):
        if max_events >= 0 and i >= max_events:
            break

        if (i + 1) % 1000 == 0:
            logger.info("Processed %d / %d events", i + 1, total_events)

        try:
            event.getByLabel("hltEgammaHLTExtra:Unseeded", eg_unseeded_handle)
            event.getByLabel("hltEgammaHLTExtra:L1Seeded", eg_l1seeded_handle)
            event.getByLabel("hltTriggerSummaryAOD::HLTX", trig_evt_handle)
        except Exception as exc:
            logger.debug("Skipping event due to getByLabel failure: %s", exc)
            continue

        if not trig_evt_handle.isValid():
            continue

        trig_evt = trig_evt_handle.product()
        eg_unseeded = eg_unseeded_handle.product() if eg_unseeded_handle.isValid() else []
        eg_l1seeded = eg_l1seeded_handle.product() if eg_l1seeded_handle.isValid() else []

        for filt in all_filters:
            filter_index = get_filter_index(trig_evt, filt)
            if filter_index >= trig_evt.sizeFilters():
                continue

            trigger_keys = trig_evt.filterKeys(filter_index)
            passed_trig_objs = []
            for key in trigger_keys:
                try:
                    passed_trig_objs.append(trig_evt.getObjects()[key])
                except Exception:
                    continue

            eg_collection = eg_l1seeded if "L1Seeded" in filt else eg_unseeded
            if not eg_collection:
                continue

            for eg_obj in eg_collection:
                matched = match_trig_objs(eg_obj.eta(), eg_obj.phi(), passed_trig_objs, max_dr=max_dr)
                if not matched:
                    continue

                sigma2uu, sigma2vv, sigma2ww = get_sigma_vars_for_filter(eg_obj, filt)
                pt = eg_obj.et()
                eta = eg_obj.eta()
                histograms[f"h_sigma2vv_{filt}"].Fill(sigma2vv)
                histograms[f"h_sigma2ww_{filt}"].Fill(sigma2ww)
                if PT_MIN_GEV <= pt <= PT_MAX_GEV and ETA_MIN <= eta <= ETA_MAX:
                    histograms[f"p2_sigma2vv_{filt}"].Fill(eta, pt, sigma2vv)
                    histograms[f"p2_sigma2ww_{filt}"].Fill(eta, pt, sigma2ww)
                    histograms[f"h2_sigma2uu_vs_pt_{filt}"].Fill(pt, sigma2uu)
                    histograms[f"h2_sigma2vv_vs_pt_{filt}"].Fill(pt, sigma2vv)
                    histograms[f"h2_sigma2ww_vs_pt_{filt}"].Fill(pt, sigma2ww)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract sigma2vv/sigma2ww per HLT filter from trigger objects."
    )
    parser.add_argument("--input-file", required=True, help="Input ROOT file")
    parser.add_argument("-o", "--output", required=True, help="Output ROOT file")
    parser.add_argument(
        "-n",
        "--max-events",
        type=int,
        default=-1,
        help="Maximum events to process (-1 means all)",
    )
    parser.add_argument(
        "--max-dr",
        type=float,
        default=MAX_DR,
        help="DeltaR for matching EG objects to filter trigger objects",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    if not os.path.isfile(args.input_file):
        raise FileNotFoundError(f"Input file does not exist: {args.input_file}")

    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.gSystem.Load("libDataFormatsFWLite.so")
    ROOT.FWLiteEnabler.enable()

    sequences = {
        "HLTEle26WP70Unseeded": (
            "cms.Sequence( hltEGL1SeedsForSingleEleIsolatedFilter + "
            "HLTDoFullUnpackingEgammaEcalSequence + "
            "HLTPFClusteringForEgammaUnseededSequence + "
            "HLTHgcalTiclPFClusteringForEgammaUnseededSequence + "
            "hltEgammaCandidatesUnseeded + hltEgammaCandidatesWrapperUnseeded + "
            "hltEG26EtUnseededFilter + hltEgammaClusterShapeUnseeded + "
            "hltEle26WP70ClusterShapeUnseededFilter + hltEgammaHGCALIDVarsUnseeded + "
            "hltEle26WP70ClusterShapeSigmavvUnseededFilter + "
            "hltEle26WP70ClusterShapeSigmawwUnseededFilter + "
            "hltEle26WP70HgcalHEUnseededFilter + HLTEGammaDoLocalHcalSequence + "
            "HLTFastJetForEgammaSequence + hltEgammaHoverEUnseeded + "
            "hltEle26WP70HEUnseededFilter + hltEgammaEcalPFClusterIsoUnseeded + "
            "hltEle26WP70EcalIsoUnseededFilter + hltEgammaHGCalLayerClusterIsoUnseeded + "
            "hltEle26WP70HgcalIsoUnseededFilter + HLTPFHcalClusteringForEgammaSequence + "
            "hltEgammaHcalPFClusterIsoUnseeded + hltEle26WP70HcalIsoUnseededFilter + "
            "HLTElePixelMatchUnseededSequence + hltEle26WP70PixelMatchUnseededFilter + "
            "hltEle26WP70PMS2UnseededFilter + HLTGsfElectronUnseededSequence + "
            "hltEle26WP70GsfOneOEMinusOneOPUnseededFilter + "
            "hltEle26WP70GsfDetaUnseededFilter + hltEle26WP70GsfDphiUnseededFilter + "
            "hltEle26WP70BestGsfNLayerITUnseededFilter + "
            "hltEle26WP70BestGsfChi2UnseededFilter + hltEgammaEleL1TrkIsoUnseeded + "
            "hltEle26WP70GsfTrackIsoFromL1TracksUnseededFilter + HLTTrackingSequence + "
            "hltEgammaEleGsfTrackIsoUnseeded + hltEle26WP70GsfTrackIsoUnseededFilter )"
        ),
        "HLTEle26WP70L1Seeded": (
            "cms.Sequence( hltEGL1SeedsForSingleEleIsolatedFilter + "
            "HLTDoFullUnpackingEgammaEcalL1SeededSequence + "
            "HLTPFClusteringForEgammaL1SeededSequence + "
            "HLTHgcalTiclPFClusteringForEgammaL1SeededSequence + "
            "hltEgammaCandidatesL1Seeded + hltEgammaCandidatesWrapperL1Seeded + "
            "hltEG26EtL1SeededFilter + hltEgammaClusterShapeL1Seeded + "
            "hltEle26WP70ClusterShapeL1SeededFilter + hltEgammaHGCALIDVarsL1Seeded + "
            "hltEle26WP70ClusterShapeSigmavvL1SeededFilter + "
            "hltEle26WP70ClusterShapeSigmawwL1SeededFilter + "
            "hltEle26WP70HgcalHEL1SeededFilter + HLTEGammaDoLocalHcalSequence + "
            "HLTFastJetForEgammaSequence + hltEgammaHoverEL1Seeded + "
            "hltEle26WP70HEL1SeededFilter + hltEgammaEcalPFClusterIsoL1Seeded + "
            "hltEle26WP70EcalIsoL1SeededFilter + hltEgammaHGCalLayerClusterIsoL1Seeded + "
            "hltEle26WP70HgcalIsoL1SeededFilter + HLTPFHcalClusteringForEgammaSequence + "
            "hltEgammaHcalPFClusterIsoL1Seeded + hltEle26WP70HcalIsoL1SeededFilter + "
            "HLTElePixelMatchL1SeededSequence + hltEle26WP70PixelMatchL1SeededFilter + "
            "hltEle26WP70PMS2L1SeededFilter + HLTGsfElectronL1SeededSequence + "
            "hltEle26WP70GsfOneOEMinusOneOPL1SeededFilter + "
            "hltEle26WP70GsfDetaL1SeededFilter + hltEle26WP70GsfDphiL1SeededFilter + "
            "hltEle26WP70BestGsfNLayerITL1SeededFilter + "
            "hltEle26WP70BestGsfChi2L1SeededFilter + hltEgammaEleL1TrkIsoL1Seeded + "
            "hltEle26WP70GsfTrackIsoFromL1TracksL1SeededFilter + HLTTrackingSequence + "
            "hltEgammaEleGsfTrackIsoL1Seeded + hltEle26WP70GsfTrackIsoL1SeededFilter )"
        ),
    }

    all_filters: List[str] = []
    for _, seq in sequences.items():
        all_filters.extend(get_filters(seq))
    all_filters = sorted(set(all_filters))
    logger.info("Total filters to analyze: %d", len(all_filters))

    events = Events([args.input_file])
    out_file = ROOT.TFile(args.output, "RECREATE")
    process_events(events, out_file, all_filters, args.max_events, args.max_dr)
    out_file.Write()
    out_file.Close()

    logger.info("Wrote per-filter sigma2vv/sigma2ww histograms to %s", args.output)


if __name__ == "__main__":
    main()

