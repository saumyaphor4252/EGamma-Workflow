#!/usr/bin/env python

# Run: python3 lastFilterEfficiency_Upgrade_QCDFakeRate.py --input-file inputFile.root -o outFile.root
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
MIN_ELE_PT = 10.
MAX_DR = 0.1
LOOSE_ID_NAMES = [
    "cutBasedElectronID-Winter22-122X-V1-loose",
]

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
    """Apply EB/EE acceptance while removing transition region."""
    abs_eta = abs(eta)
    return abs_eta <= EB_ETA_MAX or abs_eta >= EE_ETA_MIN

def passes_loose_id(ele: Any) -> bool:
    """Best-effort loose-ID for pat::Electron across collections/campaigns."""
    tried_any = False
    for id_name in LOOSE_ID_NAMES:
        try:
            value = ele.electronID(id_name)
            tried_any = True
            if value > 0.5:
                return True
        except Exception:
            continue

    # If no ID map is available in this collection, keep electron so
    # script remains usable for samples without VID maps embedded.
    if not tried_any:
        return True
    return False

def report_loose_id_availability(all_eles: List[Any]) -> None:
    """Print one-time check of loose-ID availability in this input file."""
    if len(all_eles) == 0:
        logger.warning("Loose-ID check: no electrons available in first inspected event.")
        return

    ele0 = all_eles[0]
    found_ids = []
    missing_ids = []
    for id_name in LOOSE_ID_NAMES:
        try:
            value = ele0.electronID(id_name)
            found_ids.append((id_name, value))
        except Exception:
            missing_ids.append(id_name)

    if found_ids:
        logger.info("Loose-ID check: available IDs on first electron:")
        for id_name, value in found_ids:
            logger.info("  %s = %.3f", id_name, value)
    if missing_ids:
        logger.warning("Loose-ID check: missing IDs on first electron: %s", ", ".join(missing_ids))

def get_filter_objects_trigger_event(trigger_evt: Any, filter_name: str) -> List[Any]:
    """Collect trigger::TriggerObject objects passing a given filter label."""
    passed = []
    filter_index = getFilterIndex(trigger_evt, filter_name)
    if filter_index >= trigger_evt.sizeFilters():
        return passed

    trigger_keys = trigger_evt.filterKeys(filter_index)
    for key in trigger_keys:
        try:
            passed.append(trigger_evt.getObjects()[key])
        except Exception as e:
            logger.debug(f"Error accessing trigger object key {key} for {filter_name}: {e}")
    return passed

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

def process_events(events: Events, hist_manager: HistogramManager, sequences: Dict[str, List[str]], max_events: int = -1) -> None:
    """Process events and fill histograms for multiple sequences.
    
    Args:
        events: Events to process
        hist_manager: Histogram manager instance
        sequences: Dictionary mapping sequence names to their filter lists
        max_events: Maximum number of events to process (-1 for all events)
    """
    # Initialize handles
    ele_handle, ele_label = Handle("vector<pat::Electron>"), "slimmedElectrons"
    ele_HGC_handle, ele_HGC_label = Handle("vector<pat::Electron>"), "slimmedElectronsHGC"
    ele_lowpT_handle, ele_lowpT_label = Handle("vector<pat::Electron>"), "slimmedLowPtElectrons"
    hlt_handle, hlt_label = Handle("edm::TriggerResults"), ("TriggerResults", "", "HLTX")
    triggerObjects_handle, triggerObjects_label = Handle("trigger::TriggerEvent"), "hltTriggerSummaryAOD::HLTX"
    gen_handle, gen_label = Handle("vector<reco::GenParticle>"), "prunedGenParticles"
    
    percent_step = 1
    start_time = time.time()
    total_entries = events.size() if max_events == -1 else min(events.size(), max_events)
    printed_loose_id_report = False
    
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
            event.getByLabel(ele_HGC_label, ele_HGC_handle)
            event.getByLabel(ele_lowpT_label, ele_lowpT_handle)
            event.getByLabel(hlt_label, hlt_handle)
            event.getByLabel(triggerObjects_label, triggerObjects_handle)
            event.getByLabel(gen_label, gen_handle)
        except Exception as e:
            logger.error(f"Error getting event data: {str(e)}")
            continue
            
        eles = list(ele_handle.product()) if ele_handle.isValid() else []
        eles_hgc = list(ele_HGC_handle.product()) if ele_HGC_handle.isValid() else []
        eles_lowpt = list(ele_lowpT_handle.product()) if ele_lowpT_handle.isValid() else []
        all_eles = eles + eles_hgc + eles_lowpt
        if not printed_loose_id_report:
            report_loose_id_availability(all_eles)
            printed_loose_id_report = True
        hlts = hlt_handle.product()
        trigger_objects = triggerObjects_handle.product()
        genobjs = gen_handle.product()
        
        # Process each sequence
        for sequence_name, filter_names in sequences.items():
            # Pre-build trigger objects per filter once per event for this sequence
            filter_to_trigger_objects = {}
            for filter_name in filter_names:
                filter_to_trigger_objects[filter_name] = get_filter_objects_trigger_event(trigger_objects, filter_name)

            # Process fake-electron denominator and numerator
            for eg in all_eles:
                if eg.pt() < MIN_ELE_PT:
                    continue
                if not in_accepted_eta(eg.eta()):
                    continue
                if not passes_loose_id(eg):
                    continue

                # Fake definition: not gen-matched to electron
                gen_match_ele, _, _ = match_to_gen(eg.eta(), eg.phi(), genobjs, pid=11, max_dr=MAX_DR)
                if gen_match_ele:
                    continue

                if abs(eg.eta()) <= EB_ETA_MAX:
                    hist_manager.histograms[f'{sequence_name}_den_ele_pt_EB'].Fill(eg.pt())
                if abs(eg.eta()) >= EE_ETA_MIN:
                    hist_manager.histograms[f'{sequence_name}_den_ele_pt_EE'].Fill(eg.pt())
                hist_manager.histograms[f'{sequence_name}_den_ele_pt'].Fill(eg.pt())
                hist_manager.histograms[f'{sequence_name}_den_ele_eta'].Fill(eg.eta())
                hist_manager.histograms[f'{sequence_name}_den_ele_phi'].Fill(eg.phi())

                for filter_name in filter_names:
                    matched_objs = match_trig_objs(
                        eg.eta(),
                        eg.phi(),
                        filter_to_trigger_objects[filter_name]
                    )
                    if len(matched_objs) == 0:
                        continue

                    if abs(eg.eta()) <= EB_ETA_MAX:
                        hist_manager.histograms[f'{sequence_name}_num_ele_pt_EB_{filter_name}'].Fill(eg.pt())
                    if abs(eg.eta()) >= EE_ETA_MIN:
                        hist_manager.histograms[f'{sequence_name}_num_ele_pt_EE_{filter_name}'].Fill(eg.pt())
                    hist_manager.histograms[f'{sequence_name}_num_ele_pt_{filter_name}'].Fill(eg.pt())
                    hist_manager.histograms[f'{sequence_name}_num_ele_eta_{filter_name}'].Fill(eg.eta())
                    hist_manager.histograms[f'{sequence_name}_num_ele_phi_{filter_name}'].Fill(eg.phi())

def process_single_file(input_file: str, output_file: str, sequences: Dict[str, List[str]], 
                       pt_bins: array, pt_bins_TurnOn: array, eta_bins: array, phi_bins: array, 
                       max_events: int = -1) -> bool:
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
        process_events(events, hist_manager, sequences, max_events)
        
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
    parser = argparse.ArgumentParser(description='QCD fake-electron filter efficiency analyzer (single input/output)')
    parser.add_argument('--input-file', required=True, help='Single input ROOT file')
    parser.add_argument("-o", "--output", required=True, help="Output ROOT file (written in current directory)")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("-n", "--max-events", type=int, default=-1, help="Maximum number of events to process (-1 for all events)")
    
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
        args.max_events,
    )
    if not success:
        sys.exit(1)
    print(f"\n🎉 Done. Output: {args.output}")

if __name__ == "__main__":
    main()
