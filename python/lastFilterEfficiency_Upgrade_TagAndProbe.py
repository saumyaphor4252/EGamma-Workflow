#!/usr/bin/env python

# Run: python3 lastFilterEfficiency_Upgrade_TagAndProbe.py --input-file inputFile.root -o outFile.root
# Condor: one input file → one output file in current dir (so cp *.root finds it). 

from array import array
import re
import argparse
import sys
import math
import logging
import numpy as np
from typing import List, Dict, Any
from DataFormats.FWLite import Events, Handle
import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, edm
import json
import FWCore.ParameterSet.Config as cms
from CondCore.CondDB.CondDB_cfi import *
from PhysicsTools.PythonAnalysis import *
import time
import os

# Constants
EB_ETA_MAX = 1.44
EE_ETA_MIN = 1.56
MIN_GEN_PT = 25.
MIN_TAG_PT = 35.
MIN_PROBE_PT = 10.
MAX_DR = 0.1
MAX_DR_TAG_PROBE = 0.3
Z_MASS_LOW = 60.
Z_MASS_HIGH = 120.

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class HistogramManager:
    def __init__(self, pt_bins: array, pt_bins_TurnOn: array, eta_bins: array, phi_bins: array):
        self.pt_bins = pt_bins
        self.pt_bins_TurnOn = pt_bins_TurnOn
        self.eta_bins = eta_bins
        self.phi_bins = phi_bins
        self.histograms = {}
        
    def create_histograms(self, sequences: Dict[str, List[str]]) -> None:
        """Create all necessary histograms for the analysis.
        
        Args:
            sequences: Dictionary mapping sequence names to their filter lists
        """
        # Create denominator histograms for each sequence
        for sequence_name in sequences.keys():
            self.histograms[f'{sequence_name}_den_ele_pt_EB'] = ROOT.TH1D(f"{sequence_name}_den_ele_pt_EB", "pT", len(self.pt_bins)-1, self.pt_bins)
            self.histograms[f'{sequence_name}_den_ele_pt_EE'] = ROOT.TH1D(f"{sequence_name}_den_ele_pt_EE", "pT", len(self.pt_bins)-1, self.pt_bins)
            self.histograms[f'{sequence_name}_den_ele_pt'] = ROOT.TH1D(f"{sequence_name}_den_ele_pt", "pT", len(self.pt_bins)-1, self.pt_bins)
            
            # Create eta and phi denominator histograms
            self.histograms[f'{sequence_name}_den_ele_eta'] = ROOT.TH1D(f"{sequence_name}_den_ele_eta", "eta", len(self.eta_bins)-1, self.eta_bins)
            self.histograms[f'{sequence_name}_den_ele_phi'] = ROOT.TH1D(f"{sequence_name}_den_ele_phi", "phi", len(self.phi_bins)-1, self.phi_bins)
            
            # Create numerator histograms for each filter in this sequence
            for filter_name in sequences[sequence_name]:
                self.histograms[f'{sequence_name}_num_ele_pt_EB_{filter_name}'] = ROOT.TH1D(f"{sequence_name}_num_ele_pt_EB_{filter_name}", "pT", len(self.pt_bins)-1, self.pt_bins)
                self.histograms[f'{sequence_name}_num_ele_pt_EE_{filter_name}'] = ROOT.TH1D(f"{sequence_name}_num_ele_pt_EE_{filter_name}", "pT", len(self.pt_bins)-1, self.pt_bins)
                self.histograms[f'{sequence_name}_num_ele_pt_{filter_name}'] = ROOT.TH1D(f"{sequence_name}_num_ele_pt_{filter_name}", "pT", len(self.pt_bins)-1, self.pt_bins)
                        
                # Create eta and phi numerator histograms
                self.histograms[f'{sequence_name}_num_ele_eta_{filter_name}'] = ROOT.TH1D(f"{sequence_name}_num_ele_eta_{filter_name}", "eta", len(self.eta_bins)-1, self.eta_bins)
                self.histograms[f'{sequence_name}_num_ele_phi_{filter_name}'] = ROOT.TH1D(f"{sequence_name}_num_ele_phi_{filter_name}", "phi", len(self.phi_bins)-1, self.phi_bins)

    def write_histograms(self, output_file: TFile) -> None:
        """Write all histograms to the output file."""
        for hist in self.histograms.values():
            hist.Write()

def match_trig_objs(eta: float, phi: float, trig_objs: List[Any], max_dr: float = MAX_DR) -> List[Any]:    
    """Match trigger objects to a given eta,phi position.
    
    Args:
        eta: Eta coordinate of offline electron
        phi: Phi coordinate of offline electron
        trig_objs: List of trigger objects to match against
        max_dr: Maximum delta R for matching
        
    Returns:
        List of matched trigger objects
    """
    max_dr2 = max_dr * max_dr
    matched_objs = [obj for obj in trig_objs if ROOT.reco.deltaR2(eta, phi, obj.eta(), obj.phi()) < max_dr2]
    return matched_objs

def getFilterIndex(trigEvt, filterName):
    """Get the index of a specific filter in the trigger event.
    
    Args:
        trigEvt: Trigger event object (trigger::TriggerEvent)
        filterName: Name of the filter to search for
        
    Returns:
        Index of the filter if found, otherwise sizeFilters() (invalid index)
        
    Note:
        Handles AOD filter label format by cleaning quotes and bracket notation.
        Example: "hltEle26WP70Filter"[25] -> hltEle26WP70Filter
    """
    for index in range(0, trigEvt.sizeFilters()):
        # Convert C++ string_view to Python string for comparison
        label = str(trigEvt.filterLabel(index))
        
        # Clean the label - remove quotes and bracket notation like [25]
        # Format: "hltDiEG12EtL1SeededFilter"[25] -> hltDiEG12EtL1SeededFilter
        clean_label = label.split('[')[0]  # Remove [38] part first
        clean_label = clean_label.strip('"')  # Remove surrounding quotes
        
        if filterName == clean_label:
                return index
    return trigEvt.sizeFilters()

def get_filter_objects(trig_evt: Any, filter_name: str) -> List[Any]:
    """Return trigger objects that passed a specific filter."""
    passed_objs = []
    filter_index = getFilterIndex(trig_evt, filter_name)
    if filter_index >= trig_evt.sizeFilters():
        return passed_objs

    trigger_keys = trig_evt.filterKeys(filter_index)
    for key in trigger_keys:
        try:
            passed_objs.append(trig_evt.getObjects()[key])
        except Exception as e:
            logger.debug(f"Failed to access trigger object key {key} for {filter_name}: {e}")
    return passed_objs

def get_genparts(genparts: Any, pid: int = 11, antipart: bool = True, status: int = 1) -> List[Any]:
    """Get a list of gen particles matching the given criteria.
    
    Args:
        genparts: GenParticle collection
        pid: PDG ID to match
        antipart: Whether to include antiparticles
        status: Particle status to match
        
    Returns:
        List of matching gen particles
    """
    selected = []
    if genparts is None:
        return []
        
    for part in genparts:
        pdg_id = part.pdgId()
        if pdg_id == pid or (antipart and abs(pdg_id) == abs(pid)):
            if part.isHardProcess():
                if status == 1:
                    selected.append(part)
    return selected

def match_to_gen(eta: float, phi: float, genparts: Any, pid: int = 11, 
                antipart: bool = True, max_dr: float = MAX_DR, status: int = 1) -> tuple:
    """Match an eta,phi to gen level particle.
    
    Args:
        eta: Eta coordinate
        phi: Phi coordinate
        genparts: GenParticle collection
        pid: PDG ID to match
        antipart: Whether to include antiparticles
        max_dr: Maximum delta R for matching
        status: Particle status to match
        
    Returns:
        Tuple of (best_match, best_dr2, best_pt)
    """
    best_match = None
    best_dr2 = max_dr * max_dr
    best_pt = -1
    
    selected_parts = get_genparts(genparts, pid, antipart, status)
    for part in selected_parts:
        dr2 = ROOT.reco.deltaR2(eta, phi, part.eta(), part.phi())
        if dr2 < best_dr2:
            best_match = part
            best_dr2 = dr2
            best_pt = part.pt()
            
    return best_match, best_dr2, best_pt

def in_accepted_eta(eta: float) -> bool:
    """Apply EB/EE acceptance with transition veto."""
    abs_eta = abs(eta)
    return abs_eta <= EB_ETA_MAX or abs_eta >= EE_ETA_MIN

def build_tag_probe_pairs(eles: Any, genobjs: Any, tag_objs: List[Any]) -> List[tuple]:
    """Build tag-and-probe pairs using HLT egamma objects and gen-matched DY electrons."""
    selected = []
    for eg in eles:
        gen_match, _, gen_pt = match_to_gen(eg.eta(), eg.phi(), genobjs, pid=11)
        if not gen_match:
            continue
        if not in_accepted_eta(eg.eta()):
            continue
        if gen_pt < MIN_PROBE_PT:
            continue
        selected.append((eg, gen_pt))

    pairs = []
    tag_used = set()
    for i, (tag_candidate, _) in enumerate(selected):
        if tag_candidate.pt() < MIN_TAG_PT:
            continue
        if len(match_trig_objs(tag_candidate.eta(), tag_candidate.phi(), tag_objs)) == 0:
            continue

        for j, (probe_candidate, probe_gen_pt) in enumerate(selected):
            if i == j:
                continue
            if probe_gen_pt < MIN_GEN_PT:
                continue
            if ROOT.reco.deltaR2(tag_candidate.eta(), tag_candidate.phi(), probe_candidate.eta(), probe_candidate.phi()) < MAX_DR_TAG_PROBE * MAX_DR_TAG_PROBE:
                continue

            mass = (tag_candidate.p4() + probe_candidate.p4()).M()
            if mass < Z_MASS_LOW or mass > Z_MASS_HIGH:
                continue

            pair_key = (i, j)
            if pair_key in tag_used:
                continue
            tag_used.add(pair_key)
            pairs.append((tag_candidate, probe_candidate))
    return pairs

def getFilters(cms_path: str) -> List[str]:
    """Extract filter names from a CMS path sequence.
    
    Args:
        cms_path: CMS path string or cms.Sequence object
        
    Returns:
        List of filter names
    """
    filts = []
    
    # Handle both string and cms.Sequence inputs
    if isinstance(cms_path, str):
        # Remove cms.Sequence wrapper if present
        seq_str = cms_path.replace("cms.Sequence(", "").replace(")", "")
        # Split by + and clean each module
        for module in seq_str.split("+"):
            module = module.strip()
            if "Filter" in module:
                # Remove all whitespace, process. prefix, and parentheses
                module = module.replace("process.", "").replace(" ", "").replace("(", "").replace(")", "")
                filts.append(module)
    else:
        # Handle cms.Sequence object
        # Convert sequence to string and clean it
        seq_str = str(cms_path)
        # Remove the cms.Sequence wrapper
        seq_str = seq_str.replace("cms.Sequence(", "").replace(")", "")
        # Split by + and clean each module
        for module in seq_str.split("+"):
            module = module.strip()
            if "Filter" in module:
                # Remove all whitespace, process. prefix, and parentheses
                module = module.replace("process.", "").replace(" ", "").replace("(", "").replace(")", "")
                filts.append(module)
    
    return filts

def process_events(events: Events, hist_manager: HistogramManager, sequences: Dict[str, List[str]],
                   tag_filter: str, max_events: int = -1) -> None:
    """Process events and fill histograms for multiple sequences.
    
    Args:
        events: Events to process
        hist_manager: Histogram manager instance
        sequences: Dictionary mapping sequence names to their filter lists
        max_events: Maximum number of events to process (-1 for all events)
    """
    # Initialize handles
    #ele_handle, ele_label = Handle("vector<reco::Electron>"), ("hltEgammaGsfElectronsUnseeded", "", "HLTX")
    ele_handle, ele_label = Handle("std::vector<trigger::EgammaObject>"), "hltEgammaHLTExtra:Unseeded"
    hlt_handle, hlt_label = Handle("edm::TriggerResults"), ("TriggerResults", "", "HLTX")
    hltevt_handle, hltevt_label = Handle("trigger::TriggerEvent"), "hltTriggerSummaryAOD::HLTX"
    #triggerObjects_handle, triggerObjects_label = Handle("vector<pat::TriggerObjectStandAlone>"), "slimmedPatTrigger"
    triggerObjects_handle, triggerObjects_label = Handle("trigger::TriggerEvent"), "hltTriggerSummaryAOD::HLTX"
    gen_handle, gen_label = Handle("vector<reco::GenParticle>"), ("genParticles", "", "HLT")
    
    percent_step = 1
    start_time = time.time()
    total_entries = events.size() if max_events == -1 else min(events.size(), max_events)
    
    for event_nr, event in enumerate(events):
        if max_events != -1 and event_nr >= max_events:
            break
            
        # Progress reporting
        current_percent = (event_nr + 1) / total_entries * 100
        if current_percent % percent_step == 0:
            elapsed_time = time.time()-start_time
            est_finish = "n/a"
            if event_nr!=0 or elapsed_time==0:
                remaining = float(total_entries-event_nr)/event_nr*elapsed_time 
                est_finish = time.ctime(remaining+start_time+elapsed_time)
                print("{} / {} time: {:.1f}s, est finish {}".format(event_nr+1,total_entries,elapsed_time,est_finish))
        
        # Get event data
        try:
            event.getByLabel(ele_label, ele_handle)
            event.getByLabel(hlt_label, hlt_handle)
            event.getByLabel(triggerObjects_label, triggerObjects_handle)
            event.getByLabel(gen_label, gen_handle)
        except Exception as e:
            logger.error(f"Error getting event data: {str(e)}")
            continue
            
        eles = ele_handle.product()
        hlts = hlt_handle.product()
        trigger_objects = triggerObjects_handle.product()
        genobjs = gen_handle.product()
        
        tag_passed_objects = get_filter_objects(trigger_objects, tag_filter)
        if len(tag_passed_objects) == 0:
            continue

        tag_probe_pairs = build_tag_probe_pairs(eles, genobjs, tag_passed_objects)
        if len(tag_probe_pairs) == 0:
            continue
        
        # Process each sequence
        for sequence_name, filter_names in sequences.items():
            for _, probe in tag_probe_pairs:
                if not in_accepted_eta(probe.eta()):
                    continue

                if abs(probe.eta()) <= EB_ETA_MAX:
                    hist_manager.histograms[f'{sequence_name}_den_ele_pt_EB'].Fill(probe.pt())
                if abs(probe.eta()) >= EE_ETA_MIN:
                    hist_manager.histograms[f'{sequence_name}_den_ele_pt_EE'].Fill(probe.pt())
                hist_manager.histograms[f'{sequence_name}_den_ele_pt'].Fill(probe.pt())
                hist_manager.histograms[f'{sequence_name}_den_ele_eta'].Fill(probe.eta())
                hist_manager.histograms[f'{sequence_name}_den_ele_phi'].Fill(probe.phi())

                for filter_name in filter_names:
                    passed_objs = get_filter_objects(trigger_objects, filter_name)
                    if len(passed_objs) == 0:
                        continue
                    if len(match_trig_objs(probe.eta(), probe.phi(), passed_objs)) == 0:
                        continue

                    if abs(probe.eta()) <= EB_ETA_MAX:
                        hist_manager.histograms[f'{sequence_name}_num_ele_pt_EB_{filter_name}'].Fill(probe.pt())
                    if abs(probe.eta()) >= EE_ETA_MIN:
                        hist_manager.histograms[f'{sequence_name}_num_ele_pt_EE_{filter_name}'].Fill(probe.pt())
                    hist_manager.histograms[f'{sequence_name}_num_ele_pt_{filter_name}'].Fill(probe.pt())
                    hist_manager.histograms[f'{sequence_name}_num_ele_eta_{filter_name}'].Fill(probe.eta())
                    hist_manager.histograms[f'{sequence_name}_num_ele_phi_{filter_name}'].Fill(probe.phi())

def process_single_file(input_file: str, output_file: str, sequences: Dict[str, List[str]], 
                       pt_bins: array, pt_bins_TurnOn: array, eta_bins: array, phi_bins: array, 
                       tag_filter: str, max_events: int = -1) -> bool:
    """Process a single input file and create output histograms.
    
    Args:
        input_file: Path to input ROOT file
        output_file: Path to output ROOT file
        sequences: Dictionary mapping sequence names to their filter lists
        pt_bins, pt_bins_TurnOn, eta_bins, phi_bins: Histogram binning arrays
        max_events: Maximum number of events to process (-1 for all events)
        
    Returns:
        bool: True if processing was successful, False otherwise
    """
    try:
        print(f"🔄 Processing: {os.path.basename(input_file)}")
        
        # Create histogram manager for this file
        hist_manager = HistogramManager(pt_bins, pt_bins_TurnOn, eta_bins, phi_bins)
        hist_manager.create_histograms(sequences)
        
        # Process events from this single file
        events = Events([input_file])
        process_events(events, hist_manager, sequences, tag_filter, max_events)
        
        # Write output for this file
        output_root_file = TFile(output_file, 'recreate')
        hist_manager.write_histograms(output_root_file)
        output_root_file.Close()
        
        print(f"✅ Completed: {os.path.basename(input_file)} -> {os.path.basename(output_file)}")
        return True
        
    except Exception as e:
        logger.error(f"❌ Error processing {input_file}: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description='E/gamma HLT analyzer with DY tag-and-probe (single input/output)')
    parser.add_argument('--input-file', required=True, help='Single input ROOT file')
    parser.add_argument("-o", "--output", required=True, help="Output ROOT file (written in current directory)")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("-n", "--max-events", type=int, default=-1, help="Maximum number of events to process (-1 for all events)")
    parser.add_argument("--tag-filter", default="hltEle32WPTightGsfTrackIsoL1SeededFilter",
                        help="Filter name used for tag leg preselection")
    
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    # Initialize ROOT
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.gSystem.Load("libDataFormatsFWLite.so")
    ROOT.FWLiteEnabler.enable()
    
    input_file = args.input_file
    if not os.path.isfile(input_file):
        print(f"❌ Input file not found: {input_file}")
        sys.exit(1)
    
    # Validate input file
    try:
        file_temp = ROOT.TFile(input_file)
        if file_temp.IsZombie():
            logger.error("Input file is invalid (zombie)")
            sys.exit(1)
        file_temp.Close()
    except Exception as e:
        logger.error(f"Error opening file {input_file}: {str(e)}")
        sys.exit(1)
    
    print(f"🎯 Input: {input_file} → Output: {args.output}")
    print(f"🏷️ Tag filter: {args.tag_filter}")
    
    # Setup histogram binning arrays
    # Custom pt binning: 0-50, 50-100, then steps of 100 till 3000, then steps of 250 till 4000
    pt_bins_list = [0, 50, 100]  # Start with 0-50, 50-100
    #pt_bins_list_TurnOn = [0, 20, 30, 40, 50, 60, 80, 100, 125, 150, 175, 200]
    pt_bins_list_TurnOn = [0, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 125, 150, 175, 200]
    
    # Add steps of 100 from 100 to 3000
    for i in range(200, 3001, 100):
        pt_bins_list.append(i)
    
    # Add steps of 250 from 3000 to 4000
    for i in range(3250, 4001, 250):
        pt_bins_list.append(i)
    
    pt_bins = array('d', pt_bins_list)
    pt_bins_TurnOn = array('d', pt_bins_list_TurnOn)
    eta_bins = np.arange(-4, 4.01, 0.25, dtype='d')
    phi_bins = array('d', [-3.32,-2.97,-2.62,-2.27,-1.92,-1.57,-1.22,-0.87,-0.52,-0.18,0.18,0.52,0.87,1.22,1.57,1.92,2.27,2.62,2.97,3.32])
    
    # Define sequences
    sequences = {
       'HLTEle26WP70Unseeded': "cms.Sequence( hltEGL1SeedsForSingleEleIsolatedFilter + HLTDoFullUnpackingEgammaEcalSequence + HLTPFClusteringForEgammaUnseededSequence + HLTHgcalTiclPFClusteringForEgammaUnseededSequence + hltEgammaCandidatesUnseeded + hltEgammaCandidatesWrapperUnseeded +hltEG26EtUnseededFilter + hltEgammaClusterShapeUnseeded + hltEle26WP70ClusterShapeUnseededFilter + hltEgammaHGCALIDVarsUnseeded +hltEle26WP70ClusterShapeSigmavvUnseededFilter + hltEle26WP70ClusterShapeSigmawwUnseededFilter + hltEle26WP70HgcalHEUnseededFilter +HLTEGammaDoLocalHcalSequence + HLTFastJetForEgammaSequence + hltEgammaHoverEUnseeded + hltEle26WP70HEUnseededFilter +hltEgammaEcalPFClusterIsoUnseeded + hltEle26WP70EcalIsoUnseededFilter + hltEgammaHGCalLayerClusterIsoUnseeded +hltEle26WP70HgcalIsoUnseededFilter + HLTPFHcalClusteringForEgammaSequence + hltEgammaHcalPFClusterIsoUnseeded +hltEle26WP70HcalIsoUnseededFilter + HLTElePixelMatchUnseededSequence + hltEle26WP70PixelMatchUnseededFilter +hltEle26WP70PMS2UnseededFilter + HLTGsfElectronUnseededSequence + hltEle26WP70GsfOneOEMinusOneOPUnseededFilter +hltEle26WP70GsfDetaUnseededFilter + hltEle26WP70GsfDphiUnseededFilter + hltEle26WP70BestGsfNLayerITUnseededFilter +hltEle26WP70BestGsfChi2UnseededFilter + hltEgammaEleL1TrkIsoUnseeded + hltEle26WP70GsfTrackIsoFromL1TracksUnseededFilter +HLTTrackingSequence + hltEgammaEleGsfTrackIsoUnseeded + hltEle26WP70GsfTrackIsoUnseededFilter )",
       'HLTEle26WP70L1Seeded': "cms.Sequence( hltEGL1SeedsForSingleEleIsolatedFilter + HLTDoFullUnpackingEgammaEcalL1SeededSequence + HLTPFClusteringForEgammaL1SeededSequence + HLTHgcalTiclPFClusteringForEgammaL1SeededSequence + hltEgammaCandidatesL1Seeded + hltEgammaCandidatesWrapperL1Seeded + hltEG26EtL1SeededFilter + hltEgammaClusterShapeL1Seeded+hltEle26WP70ClusterShapeL1SeededFilter + hltEgammaHGCALIDVarsL1Seeded+hltEle26WP70ClusterShapeSigmavvL1SeededFilter + hltEle26WP70ClusterShapeSigmawwL1SeededFilter + hltEle26WP70HgcalHEL1SeededFilter + HLTEGammaDoLocalHcalSequence + HLTFastJetForEgammaSequence + hltEgammaHoverEL1Seeded+hltEle26WP70HEL1SeededFilter + hltEgammaEcalPFClusterIsoL1Seeded + hltEle26WP70EcalIsoL1SeededFilter + hltEgammaHGCalLayerClusterIsoL1Seeded + hltEle26WP70HgcalIsoL1SeededFilter + HLTPFHcalClusteringForEgammaSequence + hltEgammaHcalPFClusterIsoL1Seeded + hltEle26WP70HcalIsoL1SeededFilter + HLTElePixelMatchL1SeededSequence + hltEle26WP70PixelMatchL1SeededFilter + hltEle26WP70PMS2L1SeededFilter + HLTGsfElectronL1SeededSequence + hltEle26WP70GsfOneOEMinusOneOPL1SeededFilter + hltEle26WP70GsfDetaL1SeededFilter + hltEle26WP70GsfDphiL1SeededFilter + hltEle26WP70BestGsfNLayerITL1SeededFilter + hltEle26WP70BestGsfChi2L1SeededFilter + hltEgammaEleL1TrkIsoL1Seeded + hltEle26WP70GsfTrackIsoFromL1TracksL1SeededFilter + HLTTrackingSequence + hltEgammaEleGsfTrackIsoL1Seeded + hltEle26WP70GsfTrackIsoL1SeededFilter )",
       'HLTEle32WPTightUnseeded': "cms.Sequence( hltEGL1SeedsForSingleEleIsolatedFilter + HLTDoFullUnpackingEgammaEcalSequence + HLTPFClusteringForEgammaUnseededSequence + HLTHgcalTiclPFClusteringForEgammaUnseededSequence + hltEgammaCandidatesUnseeded + hltEgammaCandidatesWrapperUnseeded + hltEG32EtUnseededFilter+hltEgammaClusterShapeUnseeded + hltEle32WPTightClusterShapeUnseededFilter + hltEgammaHGCALIDVarsUnseeded + hltEle32WPTightClusterShapeSigmavvUnseededFilter + hltEle32WPTightClusterShapeSigmawwUnseededFilter + hltEle32WPTightHgcalHEUnseededFilter + HLTEGammaDoLocalHcalSequence + HLTFastJetForEgammaSequence + hltEgammaHoverEUnseeded + hltEle32WPTightHEUnseededFilter + hltEgammaEcalPFClusterIsoUnseeded + hltEle32WPTightEcalIsoUnseededFilter + hltEgammaHGCalLayerClusterIsoUnseeded + hltEle32WPTightHgcalIsoUnseededFilter + HLTPFHcalClusteringForEgammaSequence + hltEgammaHcalPFClusterIsoUnseeded + hltEle32WPTightHcalIsoUnseededFilter + HLTElePixelMatchUnseededSequence + hltEle32WPTightPixelMatchUnseededFilter + hltEle32WPTightPMS2UnseededFilter + HLTGsfElectronUnseededSequence + hltEle32WPTightGsfOneOEMinusOneOPUnseededFilter + hltEle32WPTightGsfDetaUnseededFilter + hltEle32WPTightGsfDphiUnseededFilter + hltEle32WPTightBestGsfNLayerITUnseededFilter + hltEle32WPTightBestGsfChi2UnseededFilter + hltEgammaEleL1TrkIsoUnseeded + hltEle32WPTightGsfTrackIsoFromL1TracksUnseededFilter + HLTTrackingSequence + hltEgammaEleGsfTrackIsoUnseeded + hltEle32WPTightGsfTrackIsoUnseededFilter )",
       'HLTEle32WPTightL1Seeded': "cms.Sequence( hltEGL1SeedsForSingleEleIsolatedFilter + HLTDoFullUnpackingEgammaEcalL1SeededSequence +HLTPFClusteringForEgammaL1SeededSequence +HLTHgcalTiclPFClusteringForEgammaL1SeededSequence +hltEgammaCandidatesL1Seeded +hltEgammaCandidatesWrapperL1Seeded +hltEG32EtL1SeededFilter +hltEgammaClusterShapeL1Seeded +hltEle32WPTightClusterShapeL1SeededFilter +hltEgammaHGCALIDVarsL1Seeded +hltEle32WPTightClusterShapeSigmavvL1SeededFilter +hltEle32WPTightClusterShapeSigmawwL1SeededFilter +hltEle32WPTightHgcalHEL1SeededFilter +HLTEGammaDoLocalHcalSequence +HLTFastJetForEgammaSequence +hltEgammaHoverEL1Seeded +hltEle32WPTightHEL1SeededFilter +hltEgammaEcalPFClusterIsoL1Seeded +hltEle32WPTightEcalIsoL1SeededFilter +hltEgammaHGCalLayerClusterIsoL1Seeded +hltEle32WPTightHgcalIsoL1SeededFilter +HLTPFHcalClusteringForEgammaSequence +hltEgammaHcalPFClusterIsoL1Seeded +hltEle32WPTightHcalIsoL1SeededFilter +HLTElePixelMatchL1SeededSequence +hltEle32WPTightPixelMatchL1SeededFilter +hltEle32WPTightPMS2L1SeededFilter +HLTGsfElectronL1SeededSequence +hltEle32WPTightGsfOneOEMinusOneOPL1SeededFilter +hltEle32WPTightGsfDetaL1SeededFilter +hltEle32WPTightGsfDphiL1SeededFilter +hltEle32WPTightBestGsfNLayerITL1SeededFilter +hltEle32WPTightBestGsfChi2L1SeededFilter +hltEgammaEleL1TrkIsoL1Seeded +hltEle32WPTightGsfTrackIsoFromL1TracksL1SeededFilter +HLTTrackingSequence +hltEgammaEleGsfTrackIsoL1Seeded +hltEle32WPTightGsfTrackIsoL1SeededFilter )",
       'HLTDoubleEle25CaloIdLPMS2L1Seeded': "cms.Sequence( hltEGL1SeedsForDoubleEleNonIsolatedFilter + HLTDoFullUnpackingEgammaEcalL1SeededSequence + HLTPFClusteringForEgammaL1SeededSequence + HLTHgcalTiclPFClusteringForEgammaL1SeededSequence + hltEgammaCandidatesL1Seeded + hltEgammaCandidatesWrapperL1Seeded + hltDiEG25EtL1SeededFilter + hltEgammaClusterShapeL1Seeded + hltDiEG25CaloIdLClusterShapeL1SeededFilter + hltEgammaHGCALIDVarsL1Seeded + hltDiEG25CaloIdLClusterShapeSigmavvL1SeededFilter + hltDiEG25CaloIdLHgcalHEL1SeededFilter + HLTEGammaDoLocalHcalSequence + hltEgammaHoverEL1Seeded + hltDiEG25CaloIdLHEL1SeededFilter + HLTElePixelMatchL1SeededSequence + hltDiEle25CaloIdLPixelMatchL1SeededFilter + hltDiEle25CaloIdLPMS2L1SeededFilter )",
       'HLTDoubleEle25CaloIdLPMS2Unseeded': "cms.Sequence( HLTL1Sequence + hltEGL1SeedsForDoubleEleNonIsolatedFilter + HLTDoFullUnpackingEgammaEcalSequence + HLTPFClusteringForEgammaUnseededSequence + HLTHgcalTiclPFClusteringForEgammaUnseededSequence + hltEgammaCandidatesUnseeded + hltEgammaCandidatesWrapperUnseeded + hltDiEG25EtUnseededFilter + hltEgammaClusterShapeUnseeded + hltDiEG25CaloIdLClusterShapeUnseededFilter + hltEgammaHGCALIDVarsUnseeded + hltDiEG25CaloIdLClusterShapeSigmavvUnseededFilter + hltDiEG25CaloIdLHgcalHEUnseededFilter + HLTEGammaDoLocalHcalSequence + HLTFastJetForEgammaSequence + hltEgammaHoverEUnseeded + hltDiEG25CaloIdLHEUnseededFilter + HLTElePixelMatchUnseededSequence + hltDiEle25CaloIdLPixelMatchUnseededFilter + hltDiEle25CaloIdLPMS2UnseededFilter )",
       'HLTDoubleEle2312IsoL1Seeded': "cms.Sequence( hltEGL1SeedsForDoubleEleIsolatedFilter +HLTDoFullUnpackingEgammaEcalL1SeededSequence +HLTPFClusteringForEgammaL1SeededSequence +HLTHgcalTiclPFClusteringForEgammaL1SeededSequence +hltEgammaCandidatesL1Seeded +hltEgammaCandidatesWrapperL1Seeded +hltEG23EtL1SeededFilter +hltDiEG12EtL1SeededFilter +hltEgammaClusterShapeL1Seeded +hltDiEG2312IsoClusterShapeL1SeededFilter +hltEgammaHGCALIDVarsL1Seeded +hltDiEG2312IsoClusterShapeSigmavvL1SeededFilter +hltDiEG2312IsoClusterShapeSigmawwL1SeededFilter +hltDiEG2312IsoHgcalHEL1SeededFilter +HLTEGammaDoLocalHcalSequence +HLTFastJetForEgammaSequence +hltEgammaHoverEL1Seeded +hltDiEG2312IsoHEL1SeededFilter +hltEgammaEcalPFClusterIsoL1Seeded +hltDiEG2312IsoEcalIsoL1SeededFilter +hltEgammaHGCalLayerClusterIsoL1Seeded +hltDiEG2312IsoHgcalIsoL1SeededFilter +HLTPFHcalClusteringForEgammaSequence +hltEgammaHcalPFClusterIsoL1Seeded +hltDiEG2312IsoHcalIsoL1SeededFilter +HLTElePixelMatchL1SeededSequence +hltDiEle2312IsoPixelMatchL1SeededFilter +hltDiEle2312IsoPMS2L1SeededFilter +HLTGsfElectronL1SeededSequence +hltDiEle2312IsoGsfOneOEMinusOneOPL1SeededFilter +hltDiEle2312IsoGsfDetaL1SeededFilter +hltDiEle2312IsoGsfDphiL1SeededFilter +hltDiEle2312IsoBestGsfNLayerITL1SeededFilter +hltDiEle2312IsoBestGsfChi2L1SeededFilter +hltEgammaEleL1TrkIsoL1Seeded +hltDiEle2312IsoGsfTrackIsoFromL1TracksL1SeededFilter +HLTTrackingSequence +hltEgammaEleGsfTrackIsoL1Seeded +hltDiEle2312IsoGsfTrackIsoL1SeededFilter )"
    }
    
    # Get filter names for each sequence
    sequence_filters = {name: getFilters(seq) for name, seq in sequences.items()}
    print(f"🔍 Sequence Filters Extracted:")
    for seq_name, filters in sequence_filters.items():
        print(f"  {seq_name}: {len(filters)} filters")
    
    # Process single file and write directly to args.output (Condor-friendly: one file in cwd)
    print(f"\n🚀 Processing...")
    success = process_single_file(
        input_file,
        args.output,
        sequence_filters,
        pt_bins,
        pt_bins_TurnOn,
        eta_bins,
        phi_bins,
        args.tag_filter,
        args.max_events,
    )
    if not success:
        sys.exit(1)
    print(f"\n🎉 Done. Output: {args.output}")

if __name__ == "__main__":
    main()
