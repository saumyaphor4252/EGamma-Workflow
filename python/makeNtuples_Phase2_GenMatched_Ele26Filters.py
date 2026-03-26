#!/usr/bin/env python3
"""
Build an egHLTTree ntuple containing
  - sigma variables, kinematics, gen matching
  - One int branch (0/1) per HLT filter for Ele26WP70 Unseeded and L1Seeded paths

Collections (HLTX): hltEgammaHLTExtra Unseeded and L1Seeded.

Filter branches are named passUns_<FilterModuleName> and passL1_<FilterModuleName>
(same module names as in lastFilterEfficiency_Upgrade.py). For each row, only the
branches for the active eg_collection are filled; the other path's branches are 0.

Also: pass_sigmavv (sigma_vv filter only), eg_collection (0=Unseeded, 1=L1Seeded).

Run inside CMSSW with FWLite (same as makeNtuples_Phase2_GenMatched.py):

  cmsenv
  python3 makeNtuples_Phase2_GenMatched_Ele26Filters.py -i input.root -o out.root -n -1

Plotting sigma_vv vs trigger efficiency in (pT, eta, phi) bins (example):

  root -l out.root
  egHLTTree->Draw("eg_sigma2vv:pass_sigmavv", "eg_collection==0", "prof")
  // or bin in pt: use TProfile2D / TEfficiency with pass_sigmavv vs eg_et in eta slices

See also: plot_sigmavv_efficiency_fromNtuple.py.
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import sys
from array import array
from typing import Any, Dict, List, Optional, Tuple

import ROOT
from DataFormats.FWLite import Events, Handle

# --- Constants (aligned with lastFilterEfficiency_Upgrade.py) ---
EB_ETA_MAX = 1.44
EE_ETA_MIN = 1.56
MIN_GEN_PT_DEFAULT = 25.0
MAX_DR = 0.1

# Ele26 WP70 sequences (must match lastFilterEfficiency_Upgrade.py)
SEQUENCES: Dict[str, str] = {
    "HLTEle26WP70Unseeded": (
        "cms.Sequence( hltEGL1SeedsForSingleEleIsolatedFilter + HLTDoFullUnpackingEgammaEcalSequence + HLTPFClusteringForEgammaUnseededSequence + HLTHgcalTiclPFClusteringForEgammaUnseededSequence + hltEgammaCandidatesUnseeded + hltEgammaCandidatesWrapperUnseeded +hltEG26EtUnseededFilter + hltEgammaClusterShapeUnseeded + hltEle26WP70ClusterShapeUnseededFilter + hltEgammaHGCALIDVarsUnseeded +hltEle26WP70ClusterShapeSigmavvUnseededFilter + hltEle26WP70ClusterShapeSigmawwUnseededFilter + hltEle26WP70HgcalHEUnseededFilter +HLTEGammaDoLocalHcalSequence + HLTFastJetForEgammaSequence + hltEgammaHoverEUnseeded + hltEle26WP70HEUnseededFilter +hltEgammaEcalPFClusterIsoUnseeded + hltEle26WP70EcalIsoUnseededFilter + hltEgammaHGCalLayerClusterIsoUnseeded +hltEle26WP70HgcalIsoUnseededFilter + HLTPFHcalClusteringForEgammaSequence + hltEgammaHcalPFClusterIsoUnseeded +hltEle26WP70HcalIsoUnseededFilter + HLTElePixelMatchUnseededSequence + hltEle26WP70PixelMatchUnseededFilter +hltEle26WP70PMS2UnseededFilter + HLTGsfElectronUnseededSequence + hltEle26WP70GsfOneOEMinusOneOPUnseededFilter +hltEle26WP70GsfDetaUnseededFilter + hltEle26WP70GsfDphiUnseededFilter + hltEle26WP70BestGsfNLayerITUnseededFilter +hltEle26WP70BestGsfChi2UnseededFilter + hltEgammaEleL1TrkIsoUnseeded + hltEle26WP70GsfTrackIsoFromL1TracksUnseededFilter +HLTTrackingSequence + hltEgammaEleGsfTrackIsoUnseeded + hltEle26WP70GsfTrackIsoUnseededFilter )"
    ),
    "HLTEle26WP70L1Seeded": (
        "cms.Sequence( hltEGL1SeedsForSingleEleIsolatedFilter + HLTDoFullUnpackingEgammaEcalL1SeededSequence + HLTPFClusteringForEgammaL1SeededSequence + HLTHgcalTiclPFClusteringForEgammaL1SeededSequence + hltEgammaCandidatesL1Seeded + hltEgammaCandidatesWrapperL1Seeded + hltEG26EtL1SeededFilter + hltEgammaClusterShapeL1Seeded+hltEle26WP70ClusterShapeL1SeededFilter + hltEgammaHGCALIDVarsL1Seeded+hltEle26WP70ClusterShapeSigmavvL1SeededFilter + hltEle26WP70ClusterShapeSigmawwL1SeededFilter + hltEle26WP70HgcalHEL1SeededFilter + HLTEGammaDoLocalHcalSequence + HLTFastJetForEgammaSequence + hltEgammaHoverEL1Seeded+hltEle26WP70HEL1SeededFilter + hltEgammaEcalPFClusterIsoL1Seeded + hltEle26WP70EcalIsoL1SeededFilter + hltEgammaHGCalLayerClusterIsoL1Seeded + hltEle26WP70HgcalIsoL1SeededFilter + HLTPFHcalClusteringForEgammaSequence + hltEgammaHcalPFClusterIsoL1Seeded + hltEle26WP70HcalIsoL1SeededFilter + HLTElePixelMatchL1SeededSequence + hltEle26WP70PixelMatchL1SeededFilter + hltEle26WP70PMS2L1SeededFilter + HLTGsfElectronL1SeededSequence + hltEle26WP70GsfOneOEMinusOneOPL1SeededFilter + hltEle26WP70GsfDetaL1SeededFilter + hltEle26WP70GsfDphiL1SeededFilter + hltEle26WP70BestGsfNLayerITL1SeededFilter + hltEle26WP70BestGsfChi2L1SeededFilter + hltEgammaEleL1TrkIsoL1Seeded + hltEle26WP70GsfTrackIsoFromL1TracksL1SeededFilter + HLTTrackingSequence + hltEgammaEleGsfTrackIsoL1Seeded + hltEle26WP70GsfTrackIsoL1SeededFilter )"
    ),
}

FILTER_SIGMAVV_UNSEEDED = "hltEle26WP70ClusterShapeSigmavvUnseededFilter"
FILTER_SIGMAVV_L1SEEDED = "hltEle26WP70ClusterShapeSigmavvL1SeededFilter"


def branch_pass_uns(fn: str) -> str:
    """ROOT branch name for an Unseeded-path filter (0/1)."""
    return "passUns_" + fn


def branch_pass_l1(fn: str) -> str:
    """ROOT branch name for an L1Seeded-path filter (0/1)."""
    return "passL1_" + fn


def getFilters(cms_path: str) -> List[str]:
    """Extract filter module names from a cms.Sequence string."""
    filts: List[str] = []
    seq_str = cms_path.replace("cms.Sequence(", "").replace(")", "")
    for module in seq_str.split("+"):
        module = module.strip()
        if "Filter" in module:
            module = module.replace("process.", "").replace(" ", "").replace("(", "").replace(")", "")
            filts.append(module)
    return filts


def getFilterIndex(trigEvt: Any, filterName: str) -> int:
    """Index of filter in TriggerEvent, or sizeFilters() if missing."""
    for index in range(0, trigEvt.sizeFilters()):
        label = str(trigEvt.filterLabel(index))
        clean_label = label.split("[")[0].strip('"')
        if filterName == clean_label:
            return index
    return trigEvt.sizeFilters()


def match_trig_objs(
    eta: float, phi: float, trig_objs: List[Any], max_dr: float = MAX_DR
) -> List[Any]:
    max_dr2 = max_dr * max_dr
    return [obj for obj in trig_objs if ROOT.reco.deltaR2(eta, phi, obj.eta(), obj.phi()) < max_dr2]


def get_genparts(genparts: Any, pid: int = 11, antipart: bool = True, status: int = 1) -> List[Any]:
    selected: List[Any] = []
    if genparts is None:
        return selected
    for part in genparts:
        pdg_id = part.pdgId()
        if pdg_id == pid or (antipart and abs(pdg_id) == abs(pid)):
            if part.isHardProcess() and status == 1:
                selected.append(part)
    return selected


def match_to_gen(
    eta: float,
    phi: float,
    genparts: Any,
    pid: int = 11,
    antipart: bool = True,
    max_dr: float = MAX_DR,
    status: int = 1,
) -> Tuple[Optional[Any], float, float]:
    best_match = None
    best_dr2 = max_dr * max_dr
    best_pt = -1.0
    for part in get_genparts(genparts, pid, antipart, status):
        dr2 = ROOT.reco.deltaR2(eta, phi, part.eta(), part.phi())
        if dr2 < best_dr2:
            best_match = part
            best_dr2 = dr2
            best_pt = part.pt()
    return best_match, best_dr2, best_pt


def eg_passes_filter(
    trigger_objects: Any,
    eg_eta: float,
    eg_phi: float,
    filter_name: str,
) -> bool:
    """True if at least one trigger object at this filter matches the eg object in dR."""
    idx = getFilterIndex(trigger_objects, filter_name)
    if idx >= trigger_objects.sizeFilters():
        return False
    matched_objs: List[Any] = []
    for key in trigger_objects.filterKeys(idx):
        try:
            matched_objs.append(trigger_objects.getObjects()[key])
        except Exception:
            continue
    matched_objs = match_trig_objs(eg_eta, eg_phi, matched_objs)
    return len(matched_objs) > 0


def vn(name: str, l1: bool) -> str:
    """Map Unseeded HLT var name to L1Seeded."""
    if not l1:
        return name
    return name.replace("Unseeded", "L1Seeded")


def fill_eg_vars(obj: Any, l1: bool, out: Dict[str, "array"]) -> None:
    """Fill scalar arrays from one EgammaObject (same layout as makeNtuples_Phase2_GenMatched)."""

    def V(s: str) -> str:
        return vn(s, l1)

    out["eg_sigmaIEtaIEta"][0] = obj.var(V("hltEgammaClusterShapeUnseeded_sigmaIEtaIEta5x5"), 0)
    out["eg_ecalPFIsol_default"][0] = obj.var(V("hltEgammaEcalPFClusterIsoUnseeded"), 0)
    out["eg_hcalPFIsol_default"][0] = obj.var(V("hltEgammaHcalPFClusterIsoUnseeded"), 0)
    out["eg_hgcalPFIsol_default"][0] = obj.var(V("hltEgammaHGCalPFClusterIsoUnseeded"), 0)
    out["eg_trkIsolV0_default"][0] = obj.var(V("hltEgammaEleGsfTrackIsoUnseeded"), 0)
    out["eg_trkIsolV6_default"][0] = obj.var(V("hltEgammaEleGsfTrackIsoV6Unseeded"), 0)
    out["eg_trkIsolV72_default"][0] = obj.var(V("hltEgammaEleGsfTrackIsoV72Unseeded"), 0)
    out["eg_trkChi2_default"][0] = obj.var(V("hltEgammaGsfTrackVarsUnseeded_Chi2"), 0)
    out["eg_invESeedInvP"][0] = obj.var(V("hltEgammaGsfTrackVarsUnseeded_OneOESeedMinusOneOP"), 0)
    out["eg_invEInvP"][0] = obj.var(V("hltEgammaGsfTrackVarsUnseeded_OneOESuperMinusOneOP"), 0)
    out["eg_trkDEta"][0] = obj.var(V("hltEgammaGsfTrackVarsUnseeded_Deta"), 0)

    out["eg_sigma2uu"][0] = obj.var(V("hltEgammaHGCALIDVarsUnseeded_sigma2uu"), 0)
    out["eg_sigma2vv"][0] = obj.var(V("hltEgammaHGCALIDVarsUnseeded_sigma2vv"), 0)
    out["eg_sigma2ww"][0] = obj.var(V("hltEgammaHGCALIDVarsUnseeded_sigma2ww"), 0)
    out["eg_sigma2xx"][0] = obj.var(V("hltEgammaHGCALIDVarsUnseeded_sigma2xx"), 0)
    out["eg_sigma2xy"][0] = obj.var(V("hltEgammaHGCALIDVarsUnseeded_sigma2xy"), 0)
    out["eg_sigma2yy"][0] = obj.var(V("hltEgammaHGCALIDVarsUnseeded_sigma2yy"), 0)
    out["eg_sigma2yz"][0] = obj.var(V("hltEgammaHGCALIDVarsUnseeded_sigma2yz"), 0)
    out["eg_sigma2zx"][0] = obj.var(V("hltEgammaHGCALIDVarsUnseeded_sigma2zx"), 0)
    out["eg_sigma2zz"][0] = obj.var(V("hltEgammaHGCALIDVarsUnseeded_sigma2zz"), 0)

    out["eg_pms2_default"][0] = obj.var(V("hltEgammaPixelMatchVarsUnseeded_s2"), 0)
    out["eg_hcalHForHoverE"][0] = obj.var(V("hltEgammaHGCALIDVarsUnseeded_hForHOverE"), 0)
    out["eg_l1TrkIsoCMSSW"][0] = obj.var(V("hltEgammaHoverEUnseeded"), 0)
    out["eg_bestTrkChi2"][0] = obj.var(V("hltEgammaEleL1TrkIsoUnseeded"), 0)
    out["eg_bestTrkDEta"][0] = obj.var(V("hltEgammaBestGsfTrackVarsUnseeded_Chi2"), 0)
    out["eg_bestTrkDEtaSeed"][0] = obj.var(V("hltEgammaBestGsfTrackVarsUnseeded_Deta"), 0)
    out["eg_bestTrkDPhi"][0] = obj.var(V("hltEgammaBestGsfTrackVarsUnseeded_DetaSeed"), 0)
    out["eg_bestTrkESeedInvP"][0] = obj.var(V("hltEgammaBestGsfTrackVarsUnseeded_NLayerIT"), 0)
    out["eg_bestTrkInvEInvP"][0] = obj.var(V("hltEgammaBestGsfTrackVarsUnseeded_OneOESeedMinusOneOP"), 0)
    out["eg_hgcaliso_layerclus"][0] = obj.var(V("hltEgammaBestGsfTrackVarsUnseeded_ValidHits"), 0)
    out["eg_hgcaliso_layerclusem"][0] = obj.var(V("hltEgammaHGCalLayerClusterIsoUnseeded"), 0)
    out["eg_hgcaliso_layerclushad"][0] = obj.var(V("hltEgammaHGCalLayerClusterIsoUnseeded_em"), 0)

    s2vv = out["eg_sigma2vv"][0]
    if s2vv >= 0:
        out["sigmavv"][0] = math.sqrt(s2vv)
    else:
        out["sigmavv"][0] = -1.0


def resolve_input_files(args: argparse.Namespace) -> List[str]:
    files: List[str] = []
    if args.input_file:
        files = [args.input_file]
    elif args.input_dir:
        files = sorted(glob.glob(os.path.join(args.input_dir, "*.root")))
    elif args.input_pattern:
        files = sorted(glob.glob(args.input_pattern))
    return files


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Gen-matched Ele26 ntuple with sigma vars and per-filter pass flags."
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-i", "--input-file", help="Single input AOD/MINIAOD with HLT")
    input_group.add_argument("-d", "--input-dir", help="Directory of ROOT files")
    input_group.add_argument("-p", "--input-pattern", help='Glob, e.g. "data*.root"')
    parser.add_argument("-o", "--output", required=True, help="Output ROOT file")
    parser.add_argument("-n", "--max-events", type=int, default=-1, help="-1 = all")
    parser.add_argument("--max-dr", type=float, default=MAX_DR, help="Gen match dR")
    parser.add_argument("--min-gen-pt", type=float, default=MIN_GEN_PT_DEFAULT)
    parser.add_argument(
        "--skip-transition-region",
        action="store_true",
        default=True,
        help="Skip 1.44 < |eta| < 1.56 (default on)",
    )
    parser.add_argument(
        "--no-skip-transition-region",
        dest="skip_transition_region",
        action="store_false",
    )
    parser.add_argument(
        "--require-gen-match",
        action="store_true",
        help="Only store objects matched to gen e (hard process)",
    )
    parser.add_argument("--gen-pid", type=int, default=11)
    parser.add_argument(
        "--gen-label",
        default="genParticles",
        help="Gen collection label (module name)",
    )
    parser.add_argument(
        "--gen-instance",
        default="",
        help="Gen collection instance (often empty)",
    )
    parser.add_argument(
        "--gen-process",
        default="",
        help="Process name for genParticles (try '' or HLT to match your file)",
    )
    parser.add_argument(
        "--hlt-process",
        default="HLTX",
        help="Process name for HLT products (default HLTX)",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    input_files = resolve_input_files(args)
    if not input_files:
        print("No input files found.")
        sys.exit(1)

    # FWLite
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.gSystem.Load("libDataFormatsFWLite.so")
    ROOT.FWLiteEnabler.enable()

    filters_unseeded = getFilters(SEQUENCES["HLTEle26WP70Unseeded"])
    filters_l1seeded = getFilters(SEQUENCES["HLTEle26WP70L1Seeded"])

    # Output file
    out_file = ROOT.TFile(args.output, "RECREATE")
    tree = ROOT.TTree("egHLTTree", "EGamma + Ele26 filter flags")

    # Scalars per row
    br: Dict[str, Any] = {}
    for name, typ, init in [
        ("run", "I", 0),
        ("lumi", "I", 0),
        ("event", "I", 0),
        ("eg_collection", "I", 0),  # 0 = Unseeded, 1 = L1Seeded
        ("eg_et", "F", 0.0),
        ("eg_energy", "F", 0.0),
        ("eg_eta", "F", 0.0),
        ("eg_phi", "F", 0.0),
        ("eg_sigmaIEtaIEta", "F", 0.0),
        ("eg_ecalPFIsol_default", "F", 0.0),
        ("eg_hcalPFIsol_default", "F", 0.0),
        ("eg_hgcalPFIsol_default", "F", 0.0),
        ("eg_trkIsolV0_default", "F", 0.0),
        ("eg_trkIsolV6_default", "F", 0.0),
        ("eg_trkIsolV72_default", "F", 0.0),
        ("eg_trkChi2_default", "F", 0.0),
        ("eg_invESeedInvP", "F", 0.0),
        ("eg_invEInvP", "F", 0.0),
        ("eg_trkDEta", "F", 0.0),
        ("eg_sigma2uu", "F", 0.0),
        ("eg_sigma2vv", "F", 0.0),
        ("eg_sigma2ww", "F", 0.0),
        ("eg_sigma2xx", "F", 0.0),
        ("eg_sigma2xy", "F", 0.0),
        ("eg_sigma2yy", "F", 0.0),
        ("eg_sigma2yz", "F", 0.0),
        ("eg_sigma2zx", "F", 0.0),
        ("eg_sigma2zz", "F", 0.0),
        ("eg_pms2_default", "F", 0.0),
        ("eg_hcalHForHoverE", "F", 0.0),
        ("eg_l1TrkIsoCMSSW", "F", 0.0),
        ("eg_bestTrkChi2", "F", 0.0),
        ("eg_bestTrkDEta", "F", 0.0),
        ("eg_bestTrkDEtaSeed", "F", 0.0),
        ("eg_bestTrkDPhi", "F", 0.0),
        ("eg_bestTrkESeedInvP", "F", 0.0),
        ("eg_bestTrkInvEInvP", "F", 0.0),
        ("eg_hgcaliso_layerclus", "F", 0.0),
        ("eg_hgcaliso_layerclusem", "F", 0.0),
        ("eg_hgcaliso_layerclushad", "F", 0.0),
        ("sigmavv", "F", 0.0),
        ("gen_matched", "I", 0),
        ("gen_pt", "F", -1.0),
        ("gen_eta", "F", -999.0),
        ("gen_phi", "F", -999.0),
        ("gen_energy", "F", -1.0),
        ("gen_pdgId", "I", 0),
        ("gen_dr", "F", -1.0),
        ("pass_sigmavv", "I", 0),
    ]:
        arr_type = "i" if typ == "I" else "f"
        br[name] = array(arr_type, [init])
        tree.Branch(name, br[name], f"{name}/{typ}")

    # One int branch per filter (like other flags); prefix distinguishes Unseeded vs L1Seeded
    for fn in filters_unseeded:
        bn = branch_pass_uns(fn)
        br[bn] = array("i", [0])
        tree.Branch(bn, br[bn], f"{bn}/I")
    for fn in filters_l1seeded:
        bn = branch_pass_l1(fn)
        br[bn] = array("i", [0])
        tree.Branch(bn, br[bn], f"{bn}/I")

    # Handles
    eg_handle = Handle("std::vector<trigger::EgammaObject>")
    trig_handle = Handle("trigger::TriggerEvent")
    gen_handle = Handle("vector<reco::GenParticle>")

    labels_collections: List[Tuple[Tuple[str, str, str], bool, List[str], str]] = [
        (
            ("hltEgammaHLTExtra", "Unseeded", args.hlt_process),
            False,
            filters_unseeded,
            FILTER_SIGMAVV_UNSEEDED,
        ),
        (
            ("hltEgammaHLTExtra", "L1Seeded", args.hlt_process),
            True,
            filters_l1seeded,
            FILTER_SIGMAVV_L1SEEDED,
        ),
    ]

    events = Events(input_files)

    print(f"Filters Unseeded ({len(filters_unseeded)}): first={filters_unseeded[:3]}...")
    print(f"Filters L1Seeded ({len(filters_l1seeded)}): first={filters_l1seeded[:3]}...")
    lim = "all" if args.max_events < 0 else str(args.max_events)
    print(f"Processing up to {lim} events from {len(input_files)} file(s).")

    for ev_i, event in enumerate(events):
        if args.max_events >= 0 and ev_i >= args.max_events:
            break

        event.getByLabel("hltTriggerSummaryAOD::" + args.hlt_process, trig_handle)
        if not trig_handle.isValid():
            if args.verbose:
                print(f"Event {ev_i}: invalid TriggerEvent")
            continue
        trig_ev = trig_handle.product()

        event.getByLabel((args.gen_label, args.gen_instance, args.gen_process), gen_handle)
        genparts = gen_handle.product() if gen_handle.isValid() else None

        for coll_tag, is_l1, filt_list, sigmavv_filter_name in labels_collections:
            event.getByLabel(coll_tag, eg_handle)
            if not eg_handle.isValid():
                if args.verbose:
                    print(f"  missing E/gamma collection {coll_tag}")
                continue

            egs = eg_handle.product()
            for j in range(egs.size()):
                eg = egs[j]
                eta = eg.eta()
                phi = eg.phi()
                pt = eg.pt()

                if args.skip_transition_region:
                    if abs(eta) > EB_ETA_MAX and abs(eta) < EE_ETA_MIN:
                        continue

                gm = 0
                gpt = -1.0
                geta = -999.0
                gphi = -999.0
                gen_e = -1.0
                gpdg = 0
                gdr = -1.0

                if genparts is not None:
                    best, _, _ = match_to_gen(
                        eta, phi, genparts, pid=args.gen_pid, max_dr=args.max_dr
                    )
                    if best is not None:
                        gdr = ROOT.reco.deltaR(eta, phi, best.eta(), best.phi())
                        gpt = best.pt()
                        geta = best.eta()
                        gphi = best.phi()
                        gen_e = best.energy()
                        gpdg = best.pdgId()
                        gm = 1
                        if gpt < args.min_gen_pt and args.require_gen_match:
                            continue
                    elif args.require_gen_match:
                        continue
                elif args.require_gen_match:
                    continue

                br["run"][0] = event.eventAuxiliary().run()
                br["lumi"][0] = event.eventAuxiliary().luminosityBlock()
                br["event"][0] = event.eventAuxiliary().event()
                br["eg_collection"][0] = 1 if is_l1 else 0
                br["eg_et"][0] = pt
                br["eg_energy"][0] = eg.energy()
                br["eg_eta"][0] = eta
                br["eg_phi"][0] = phi

                fill_eg_vars(eg, is_l1, br)

                br["gen_matched"][0] = gm
                br["gen_pt"][0] = gpt
                br["gen_eta"][0] = geta
                br["gen_phi"][0] = gphi
                br["gen_energy"][0] = gen_e
                br["gen_pdgId"][0] = gpdg
                br["gen_dr"][0] = gdr

                # Reset all filter branches; then fill the path that matches this row
                for fn in filters_unseeded:
                    br[branch_pass_uns(fn)][0] = 0
                for fn in filters_l1seeded:
                    br[branch_pass_l1(fn)][0] = 0
                for fn in filt_list:
                    ok = 1 if eg_passes_filter(trig_ev, eta, phi, fn) else 0
                    if is_l1:
                        br[branch_pass_l1(fn)][0] = ok
                    else:
                        br[branch_pass_uns(fn)][0] = ok

                br["pass_sigmavv"][0] = (
                    1 if eg_passes_filter(trig_ev, eta, phi, sigmavv_filter_name) else 0
                )

                tree.Fill()

    out_file.Write()
    out_file.Close()
    print(f"Done. Wrote {args.output}")


if __name__ == "__main__":
    main()
