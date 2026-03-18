#!/usr/bin/env python3

import ROOT
from DataFormats.FWLite import Events, Handle
import argparse
import glob
import os
import math

def calculate_invariant_mass(e1, e2):
    """Calculate invariant mass of two electrons using their 4-momenta"""
    # Get 4-momenta
    eta1, phi1, e1_val = e1.eta(), e1.phi(), e1.energy()
    eta2, phi2, e2_val = e2.eta(), e2.phi(), e2.energy()
    
    # Calculate px, py, pz for each electron
    pt1 = e1.et() / math.cosh(eta1) if abs(eta1) < 10 else e1.et()
    pt2 = e2.et() / math.cosh(eta2) if abs(eta2) < 10 else e2.et()
    
    px1 = pt1 * math.cos(phi1)
    py1 = pt1 * math.sin(phi1)
    pz1 = pt1 * math.sinh(eta1)
    
    px2 = pt2 * math.cos(phi2)
    py2 = pt2 * math.sin(phi2)
    pz2 = pt2 * math.sinh(eta2)
    
    # Calculate invariant mass: M^2 = (E1 + E2)^2 - (p1 + p2)^2
    e_tot = e1_val + e2_val
    px_tot = px1 + px2
    py_tot = py1 + py2
    pz_tot = pz1 + pz2
    
    m2 = e_tot * e_tot - (px_tot * px_tot + py_tot * py_tot + pz_tot * pz_tot)
    return math.sqrt(max(0.0, m2))

def get_charge(obj):
    """Try to get charge from electron object"""
    try:
        # Try to get charge from GSF track if available
        if hasattr(obj, 'gsfTrack') and obj.gsfTrack().isNonnull():
            return obj.gsfTrack().charge()
        # Try to get from best track
        if hasattr(obj, 'bestTrack') and obj.bestTrack().isNonnull():
            return obj.bestTrack().charge()
    except:
        pass
    return 0  # Return 0 if charge cannot be determined

def is_tag_electron(obj, tag_et_min=30.0):
    """Tag selection: tight criteria for well-identified electrons"""
    # Basic kinematic cuts
    if obj.et() < tag_et_min:
        return False
    
    # Additional tag criteria can be added here
    # For example: tight ID cuts, isolation requirements, etc.
    # For now, we use ET cut as the main tag requirement
    
    return True

def is_probe_electron(obj, probe_et_min=10.0):
    """Probe selection: looser criteria"""
    # Basic kinematic cuts
    if obj.et() < probe_et_min:
        return False
    
    # Probe can have looser selection - just basic ET requirement
    return True

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Create debug ntuples for E/Gamma trigger objects with Tag-and-Probe selection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 makeNtuples_Phase2_TagProbe.py input.root -o output.root
  python3 makeNtuples_Phase2_TagProbe.py --output output.root input.root
  python3 makeNtuples_Phase2_TagProbe.py --max-events 100 -o output.root input.root
  python3 makeNtuples_Phase2_TagProbe.py --input-dir /path/to/files/ -o output.root
  python3 makeNtuples_Phase2_TagProbe.py --input-pattern "*.root" -o output.root
  python3 makeNtuples_Phase2_TagProbe.py --help
        """
    )
    
    # Add new input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("input_file", nargs='?',
                           help="Input ROOT file containing E/Gamma trigger objects")
    input_group.add_argument("--input-dir", "-d",
                           help="Directory containing ROOT files to process")
    input_group.add_argument("--input-pattern", "-p",
                           help="Pattern to match ROOT files (e.g., '*.root')")
    
    parser.add_argument("--output", "-o", 
                       required=True,
                       help="Output ROOT file for debug ntuples")
    parser.add_argument("--max-events", "-n", 
                       type=int, default=10,
                       help="Maximum number of events to process (default: 10)")
    parser.add_argument("--verbose", "-v", 
                       action="store_true",
                       help="Enable verbose output")
    parser.add_argument("--tag-et-min", 
                       type=float, default=30.0,
                       help="Minimum ET for tag electron (default: 30.0 GeV)")
    parser.add_argument("--probe-et-min", 
                       type=float, default=10.0,
                       help="Minimum ET for probe electron (default: 10.0 GeV)")
    parser.add_argument("--mass-min", 
                       type=float, default=60.0,
                       help="Minimum invariant mass for tag-probe pair (default: 60.0 GeV)")
    parser.add_argument("--mass-max", 
                       type=float, default=120.0,
                       help="Maximum invariant mass for tag-probe pair (default: 120.0 GeV)")
    
    args = parser.parse_args()
    
    # Determine input files
    input_files = []
    if args.input_file:
        input_files = [args.input_file]
    elif args.input_dir:
        # Handle EOS paths better
        if "eos" in args.input_dir:
            # Use shell expansion for EOS paths
            import subprocess
            try:
                result = subprocess.run(f"ls {args.input_dir}/*.root", 
                                     capture_output=True, text=True, shell=True)
                if result.returncode == 0:
                    input_files = result.stdout.strip().split('\n')
                else:
                    # Fallback to glob
                    input_files = glob.glob(os.path.join(args.input_dir, "*.root"))
            except:
                input_files = glob.glob(os.path.join(args.input_dir, "*.root"))
        else:
            input_files = glob.glob(os.path.join(args.input_dir, "*.root"))
    elif args.input_pattern:
        # Handle EOS paths in patterns
        if "eos" in args.input_pattern:
            import subprocess
            try:
                result = subprocess.run(f"ls {args.input_pattern}", 
                                     capture_output=True, text=True, shell=True)
                if result.returncode == 0:
                    input_files = result.stdout.strip().split('\n')
                else:
                    input_files = glob.glob(args.input_pattern)
            except:
                input_files = glob.glob(args.input_pattern)
        else:
            input_files = glob.glob(args.input_pattern)
    
    if not input_files:
        print("❌ No input files found!")
        print("💡 Debugging info:")
        if args.input_dir:
            print(f"   Directory: {args.input_dir}")
            print(f"   Directory exists: {os.path.exists(args.input_dir)}")
            if os.path.exists(args.input_dir):
                print(f"   Directory contents:")
                try:
                    for item in os.listdir(args.input_dir):
                        print(f"     {item}")
                except Exception as e:
                    print(f"     Error listing directory: {e}")
        elif args.input_pattern:
            print(f"   Pattern: {args.input_pattern}")
        return
    
    # Sort files for consistent processing order
    input_files.sort()
    
    output_file = args.output
    max_events = args.max_events
    verbose = args.verbose
    tag_et_min = args.tag_et_min
    probe_et_min = args.probe_et_min
    mass_min = args.mass_min
    mass_max = args.mass_max
    
    print(f"📁 Found {len(input_files)} input files:")
    for f in input_files:
        print(f"   {f}")
    print(f"💾 Output will be saved to: {output_file}")
    print(f"⚡ Processing up to {max_events} events")
    print(f"🏷️  Tag selection: ET > {tag_et_min} GeV")
    print(f"🔍 Probe selection: ET > {probe_et_min} GeV")
    print(f"📊 Mass window: {mass_min} - {mass_max} GeV")
    if verbose:
        print("🔍 Verbose mode enabled")
    
    # Create handles for the collections
    egtrigobjs_handle = Handle("std::vector<trigger::EgammaObject>")
    egtrigobjs_unseeded_handle = Handle("std::vector<trigger::EgammaObject>")
    
    # Create output file and tree
    out_file = ROOT.TFile(output_file, "RECREATE")
    tree = ROOT.TTree("egHLTTree", "EGamma Tree - Tag and Probe")
    
    # Create variables for the tree
    run = ROOT.std.vector("int")()
    lumi = ROOT.std.vector("int")()
    event = ROOT.std.vector("int")()
    
    # E/Gamma object variables (for probe electrons only)
    eg_et = ROOT.std.vector("float")()
    eg_energy = ROOT.std.vector("float")()
    eg_eta = ROOT.std.vector("float")()
    eg_phi = ROOT.std.vector("float")()
    eg_rawEnergy = ROOT.std.vector("float")()
    eg_nrClus = ROOT.std.vector("int")()
    eg_phiWidth = ROOT.std.vector("float")()
    eg_seedId = ROOT.std.vector("unsigned int")()
    eg_seedDet = ROOT.std.vector("int")()
    eg_sigmaIEtaIEta = ROOT.std.vector("float")()
    
    # HLT Isolation variables
    eg_ecalPFIsol_default = ROOT.std.vector("float")()
    eg_hcalPFIsol_default = ROOT.std.vector("float")()
    eg_hgcalPFIsol_default = ROOT.std.vector("float")()
    eg_trkIsolV0_default = ROOT.std.vector("float")()
    eg_trkIsolV6_default = ROOT.std.vector("float")()
    eg_trkIsolV72_default = ROOT.std.vector("float")()
    
    # Track variables
    eg_trkChi2_default = ROOT.std.vector("float")()
    eg_trkMissHits = ROOT.std.vector("int")()
    eg_trkValidHits = ROOT.std.vector("int")()
    eg_invESeedInvP = ROOT.std.vector("float")()
    eg_invEInvP = ROOT.std.vector("float")()
    eg_trkDEta = ROOT.std.vector("float")()
    eg_trkDEtaSeed = ROOT.std.vector("float")()
    eg_trkDPhi = ROOT.std.vector("float")()
    eg_trkNrLayerIT = ROOT.std.vector("int")()
    
    # HGCAL ID variables
    eg_rVar = ROOT.std.vector("float")()
    eg_sigma2uu = ROOT.std.vector("float")()
    eg_sigma2vv = ROOT.std.vector("float")()
    eg_sigma2ww = ROOT.std.vector("float")()
    eg_sigma2xx = ROOT.std.vector("float")()
    eg_sigma2xy = ROOT.std.vector("float")()
    eg_sigma2yy = ROOT.std.vector("float")()
    eg_sigma2yz = ROOT.std.vector("float")()
    eg_sigma2zx = ROOT.std.vector("float")()
    eg_sigma2zz = ROOT.std.vector("float")()
    
    # Pixel match and other variables
    eg_pms2_default = ROOT.std.vector("float")()
    eg_hcalHForHoverE = ROOT.std.vector("float")()
    eg_l1TrkIsoCMSSW = ROOT.std.vector("float")()
    eg_bestTrkChi2 = ROOT.std.vector("float")()
    eg_bestTrkDEta = ROOT.std.vector("float")()
    eg_bestTrkDEtaSeed = ROOT.std.vector("float")()
    eg_bestTrkDPhi = ROOT.std.vector("float")()
    eg_bestTrkMissHits = ROOT.std.vector("int")()
    eg_bestTrkNrLayerIT = ROOT.std.vector("int")()
    eg_bestTrkESeedInvP = ROOT.std.vector("float")()
    eg_bestTrkInvEInvP = ROOT.std.vector("float")()
    eg_bestTrkValitHits = ROOT.std.vector("int")()
    
    # HGCAL isolation variables
    eg_hgcaliso_layerclus = ROOT.std.vector("float")()
    eg_hgcaliso_layerclusem = ROOT.std.vector("float")()
    eg_hgcaliso_layerclushad = ROOT.std.vector("float")()
    
    # Legacy variables (keeping for compatibility)
    eg_ecalPFIsol = ROOT.std.vector("float")()
    eg_hcalPFIsol = ROOT.std.vector("float")()
    eg_trkIsolV6 = ROOT.std.vector("float")()
    
    # Tag-Probe specific variables
    tag_et = ROOT.std.vector("float")()
    tag_eta = ROOT.std.vector("float")()
    tag_phi = ROOT.std.vector("float")()
    tag_charge = ROOT.std.vector("int")()
    probe_charge = ROOT.std.vector("int")()
    inv_mass = ROOT.std.vector("float")()
    
    # Collection info
    collection_name = ROOT.std.vector("string")()
    nr_objects = ROOT.std.vector("int")()
    
    # Create branches
    tree.Branch("run", run)
    tree.Branch("lumi", lumi)
    tree.Branch("event", event)
    tree.Branch("eg_et", eg_et)
    tree.Branch("eg_energy", eg_energy)
    tree.Branch("eg_eta", eg_eta)
    tree.Branch("eg_phi", eg_phi)
    tree.Branch("eg_rawEnergy", eg_rawEnergy)
    tree.Branch("eg_nrClus", eg_nrClus)
    tree.Branch("eg_phiWidth", eg_phiWidth)
    tree.Branch("eg_seedId", eg_seedId)
    tree.Branch("eg_seedDet", eg_seedDet)
    tree.Branch("eg_sigmaIEtaIEta", eg_sigmaIEtaIEta)
    
    # HLT Isolation branches
    tree.Branch("eg_ecalPFIsol_default", eg_ecalPFIsol_default)
    tree.Branch("eg_hcalPFIsol_default", eg_hcalPFIsol_default)
    tree.Branch("eg_hgcalPFIsol_default", eg_hgcalPFIsol_default)
    tree.Branch("eg_trkIsolV0_default", eg_trkIsolV0_default)
    tree.Branch("eg_trkIsolV6_default", eg_trkIsolV6_default)
    tree.Branch("eg_trkIsolV72_default", eg_trkIsolV72_default)
    
    # Track variable branches
    tree.Branch("eg_trkChi2_default", eg_trkChi2_default)
    tree.Branch("eg_trkMissHits", eg_trkMissHits)
    tree.Branch("eg_trkValidHits", eg_trkValidHits)
    tree.Branch("eg_invESeedInvP", eg_invESeedInvP)
    tree.Branch("eg_invEInvP", eg_invEInvP)
    tree.Branch("eg_trkDEta", eg_trkDEta)
    tree.Branch("eg_trkDEtaSeed", eg_trkDEtaSeed)
    tree.Branch("eg_trkDPhi", eg_trkDPhi)
    tree.Branch("eg_trkNrLayerIT", eg_trkNrLayerIT)
    
    # HGCAL ID variable branches
    tree.Branch("eg_rVar", eg_rVar)
    tree.Branch("eg_sigma2uu", eg_sigma2uu)
    tree.Branch("eg_sigma2vv", eg_sigma2vv)
    tree.Branch("eg_sigma2ww", eg_sigma2ww)
    tree.Branch("eg_sigma2xx", eg_sigma2xx)
    tree.Branch("eg_sigma2xy", eg_sigma2xy)
    tree.Branch("eg_sigma2yy", eg_sigma2yy)
    tree.Branch("eg_sigma2yz", eg_sigma2yz)
    tree.Branch("eg_sigma2zx", eg_sigma2zx)
    tree.Branch("eg_sigma2zz", eg_sigma2zz)
    
    # Pixel match and other variable branches
    tree.Branch("eg_pms2_default", eg_pms2_default)
    tree.Branch("eg_hcalHForHoverE", eg_hcalHForHoverE)
    tree.Branch("eg_l1TrkIsoCMSSW", eg_l1TrkIsoCMSSW)
    tree.Branch("eg_bestTrkChi2", eg_bestTrkChi2)
    tree.Branch("eg_bestTrkDEta", eg_bestTrkDEta)
    tree.Branch("eg_bestTrkDEtaSeed", eg_bestTrkDEtaSeed)
    tree.Branch("eg_bestTrkDPhi", eg_bestTrkDPhi)
    tree.Branch("eg_bestTrkMissHits", eg_bestTrkMissHits)
    tree.Branch("eg_bestTrkNrLayerIT", eg_bestTrkNrLayerIT)
    tree.Branch("eg_bestTrkESeedInvP", eg_bestTrkESeedInvP)
    tree.Branch("eg_bestTrkInvEInvP", eg_bestTrkInvEInvP)
    tree.Branch("eg_bestTrkValitHits", eg_bestTrkValitHits)
    
    # HGCAL isolation branches
    tree.Branch("eg_hgcaliso_layerclus", eg_hgcaliso_layerclus)
    tree.Branch("eg_hgcaliso_layerclusem", eg_hgcaliso_layerclusem)
    tree.Branch("eg_hgcaliso_layerclushad", eg_hgcaliso_layerclushad)
    
    # Legacy branches (keeping for compatibility)
    tree.Branch("eg_ecalPFIsol", eg_ecalPFIsol)
    tree.Branch("eg_hcalPFIsol", eg_hcalPFIsol)
    tree.Branch("eg_trkIsolV6", eg_trkIsolV6)
    
    # Tag-Probe branches
    tree.Branch("tag_et", tag_et)
    tree.Branch("tag_eta", tag_eta)
    tree.Branch("tag_phi", tag_phi)
    tree.Branch("tag_charge", tag_charge)
    tree.Branch("probe_charge", probe_charge)
    tree.Branch("inv_mass", inv_mass)
    
    tree.Branch("collection_name", collection_name)
    tree.Branch("nr_objects", nr_objects)
    
    # Open events
    events = Events(input_files)
    
    print("Processing events...")
    for i, evt in enumerate(events):
        if i >= max_events:  # Process up to max_events
            break
            
        print(f"\n=== Event {i} ===")
        print(f"Run: {evt.eventAuxiliary().run()}, Lumi: {evt.eventAuxiliary().luminosityBlock()}, Event: {evt.eventAuxiliary().event()}")
                
        # Clear vectors for this event
        run.clear()
        lumi.clear()
        event.clear()
        eg_et.clear()
        eg_energy.clear()
        eg_eta.clear()
        eg_phi.clear()
        eg_rawEnergy.clear()
        eg_nrClus.clear()
        eg_phiWidth.clear()
        eg_seedId.clear()
        eg_seedDet.clear()
        eg_sigmaIEtaIEta.clear()
        
        # Clear HLT Isolation vectors
        eg_ecalPFIsol_default.clear()
        eg_hcalPFIsol_default.clear()
        eg_hgcalPFIsol_default.clear()
        eg_trkIsolV0_default.clear()
        eg_trkIsolV6_default.clear()
        eg_trkIsolV72_default.clear()
        
        # Clear Track variable vectors
        eg_trkChi2_default.clear()
        eg_trkMissHits.clear()
        eg_trkValidHits.clear()
        eg_invESeedInvP.clear()
        eg_invEInvP.clear()
        eg_trkDEta.clear()
        eg_trkDEtaSeed.clear()
        eg_trkDPhi.clear()
        eg_trkNrLayerIT.clear()
        
        # Clear HGCAL ID variable vectors
        eg_rVar.clear()
        eg_sigma2uu.clear()
        eg_sigma2vv.clear()
        eg_sigma2ww.clear()
        eg_sigma2xx.clear()
        eg_sigma2xy.clear()
        eg_sigma2yy.clear()
        eg_sigma2yz.clear()
        eg_sigma2zx.clear()
        eg_sigma2zz.clear()
        
        # Clear Pixel match and other variable vectors
        eg_pms2_default.clear()
        eg_hcalHForHoverE.clear()
        eg_l1TrkIsoCMSSW.clear()
        eg_bestTrkChi2.clear()
        eg_bestTrkDEta.clear()
        eg_bestTrkDEtaSeed.clear()
        eg_bestTrkDPhi.clear()
        eg_bestTrkMissHits.clear()
        eg_bestTrkNrLayerIT.clear()
        eg_bestTrkESeedInvP.clear()
        eg_bestTrkInvEInvP.clear()
        eg_bestTrkValitHits.clear()
        
        # Clear HGCAL isolation vectors
        eg_hgcaliso_layerclus.clear()
        eg_hgcaliso_layerclusem.clear()
        eg_hgcaliso_layerclushad.clear()
        
        # Clear Legacy vectors
        eg_ecalPFIsol.clear()
        eg_hcalPFIsol.clear()
        eg_trkIsolV6.clear()
        
        # Clear Tag-Probe vectors
        tag_et.clear()
        tag_eta.clear()
        tag_phi.clear()
        tag_charge.clear()
        probe_charge.clear()
        inv_mass.clear()
        
        collection_name.clear()
        nr_objects.clear()
        
        # Add event info
        run.push_back(evt.eventAuxiliary().run())
        lumi.push_back(evt.eventAuxiliary().luminosityBlock())
        event.push_back(evt.eventAuxiliary().event())
        
        # Use the working collections found by check_collections.py
        collections_to_try = [
            ("egtrigobjs_unseeded", "hltEgammaHLTExtra:Unseeded"),
            ("egtrigobjs_unseeded", "hltEgammaHLTExtra:Unseeded", "HLTX")
        ]
        
        for collection_info in collections_to_try:
            if len(collection_info) == 2:
                col_name, label = collection_info
                process_name = ""
            else:
                col_name, label, process_name = collection_info
                
            print(f"\nTrying collection: {col_name} with label: {label}" + (f" process: {process_name}" if process_name else ""))
            
            try:
                if col_name == "egtrigobjs":
                    handle = egtrigobjs_handle
                else:
                    handle = egtrigobjs_unseeded_handle
                
                if process_name:
                    evt.getByLabel(label, process_name, handle)
                else:
                    evt.getByLabel(label, handle)
                
                if handle.isValid():
                    egobjs = handle.product()
                    total_objects = egobjs.size()
                    print(f"  ✓ Found {total_objects} objects")
                    
                    if total_objects < 2:
                        if verbose:
                            print(f"    → Need at least 2 objects for tag-probe, skipping event")
                        continue
                    
                    # Identify tag and probe candidates
                    tags = []
                    probes = []
                    
                    for j, obj in enumerate(egobjs):
                        obj_charge = get_charge(obj)
                        if is_tag_electron(obj, tag_et_min):
                            tags.append((j, obj, obj_charge))
                        if is_probe_electron(obj, probe_et_min):
                            probes.append((j, obj, obj_charge))
                    
                    if verbose:
                        print(f"    → Found {len(tags)} tag candidates and {len(probes)} probe candidates")
                    
                    # Match tag-probe pairs
                    probe_count = 0
                    for tag_idx, tag_obj, tag_chg in tags:
                        for probe_idx, probe_obj, probe_chg in probes:
                            # Skip if same object
                            if tag_idx == probe_idx:
                                continue
                            
                            # Require opposite charge (if charge is available)
                            if tag_chg != 0 and probe_chg != 0:
                                if tag_chg * probe_chg > 0:  # Same sign
                                    if verbose:
                                        print(f"    → Skipping pair: same charge (tag={tag_chg}, probe={probe_chg})")
                                    continue
                            
                            # Calculate invariant mass
                            mass = calculate_invariant_mass(tag_obj, probe_obj)
                            
                            # Check mass window
                            if mass < mass_min or mass > mass_max:
                                if verbose:
                                    print(f"    → Skipping pair: mass={mass:.2f} GeV (outside {mass_min}-{mass_max} GeV window)")
                                continue
                            
                            # Found a valid tag-probe pair! Store probe electron
                            if verbose:
                                print(f"    → Valid tag-probe pair found: mass={mass:.2f} GeV, tag ET={tag_obj.et():.2f} GeV, probe ET={probe_obj.et():.2f} GeV")
                            
                            probe_count += 1
                            
                            # Store tag information
                            tag_et.push_back(tag_obj.et())
                            tag_eta.push_back(tag_obj.eta())
                            tag_phi.push_back(tag_obj.phi())
                            tag_charge.push_back(tag_chg)
                            probe_charge.push_back(probe_chg)
                            inv_mass.push_back(mass)
                            
                            # Store probe electron properties
                            eg_et.push_back(probe_obj.et())
                            eg_energy.push_back(probe_obj.energy())
                            eg_eta.push_back(probe_obj.eta())
                            eg_phi.push_back(probe_obj.phi())
                            
                            # SuperCluster properties
                            eg_nrClus.push_back(probe_obj.superCluster().clusters().size())
                            eg_rawEnergy.push_back(probe_obj.superCluster().rawEnergy())
                            eg_phiWidth.push_back(probe_obj.superCluster().phiWidth())
                            eg_seedId.push_back(probe_obj.superCluster().seed().seed().rawId())
                            eg_seedDet.push_back(probe_obj.superCluster().seed().seed().det())                                                                           
                                                                   
                            # HLT variables - try different naming conventions
                            eg_sigmaIEtaIEta.push_back(probe_obj.var("hltEgammaClusterShapeUnseeded_sigmaIEtaIEta5x5",0))
                            eg_ecalPFIsol_default.push_back(probe_obj.var("hltEgammaEcalPFClusterIsoUnseeded",0))
                            eg_hcalPFIsol_default.push_back(probe_obj.var("hltEgammaHcalPFClusterIsoUnseeded",0))
                            eg_hgcalPFIsol_default.push_back(probe_obj.var("hltEgammaHGCalPFClusterIsoUnseeded",0))
                            eg_trkIsolV0_default.push_back(probe_obj.var("hltEgammaEleGsfTrackIsoUnseeded",0))
                            eg_trkIsolV6_default.push_back(probe_obj.var("hltEgammaEleGsfTrackIsoV6Unseeded",0))
                            eg_trkIsolV72_default.push_back(probe_obj.var("hltEgammaEleGsfTrackIsoV72Unseeded",0))
                            eg_trkChi2_default.push_back(probe_obj.var("hltEgammaGsfTrackVarsUnseeded_Chi2",0))
                            #eg_trkMissHits.push_back(probe_obj.var("hltEgammaGsfTrackVarsUnseeded_MissingHits",0))
                            #eg_trkValidHits.push_back(probe_obj.var("hltEgammaGsfTrackVarsUnseeded_ValidHits",0))
                            eg_invESeedInvP.push_back(probe_obj.var("hltEgammaGsfTrackVarsUnseeded_OneOESeedMinusOneOP",0))
                            eg_invEInvP.push_back(probe_obj.var("hltEgammaGsfTrackVarsUnseeded_OneOESuperMinusOneOP",0))
                            eg_trkDEta.push_back(probe_obj.var("hltEgammaGsfTrackVarsUnseeded_Deta",0))
                            #eg_trkDEtaSeed.push_back(probe_obj.var("hltEgammaGsfTrackVarsUnseeded_DetaSeed",0))
                            #eg_trkDPhi.push_back(probe_obj.var("hltEgammaGsfTrackVarsUnseeded_Dphi",0))
                            #eg_trkNrLayerIT.push_back(probe_obj.var("hltEgammaGsfTrackVarsUnseeded_NLayerIT",0))
                            #eg_rVar.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_rVar",0))
                            
                            eg_sigma2uu.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2uu",0))
                            eg_sigma2vv.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2vv",0))
                            eg_sigma2ww.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2ww",0))
                            eg_sigma2xx.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2xx",0))
                            eg_sigma2xy.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2xy",0))
                            eg_sigma2yy.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2yy",0))
                            eg_sigma2yz.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2yz",0))
                            eg_sigma2zx.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2zx",0))
                            eg_sigma2zz.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2zz",0))
                            eg_pms2_default.push_back(probe_obj.var("hltEgammaPixelMatchVarsUnseeded_s2",0))
                            eg_hcalHForHoverE.push_back(probe_obj.var("hltEgammaHGCALIDVarsUnseeded_hForHOverE",0))
                            
                            eg_l1TrkIsoCMSSW.push_back(probe_obj.var("hltEgammaHoverEUnseeded",0))
                            eg_bestTrkChi2.push_back(probe_obj.var("hltEgammaEleL1TrkIsoUnseeded",0))
                            eg_bestTrkDEta.push_back(probe_obj.var("hltEgammaBestGsfTrackVarsUnseeded_Chi2",0))
                            eg_bestTrkDEtaSeed.push_back(probe_obj.var("hltEgammaBestGsfTrackVarsUnseeded_Deta",0))
                            eg_bestTrkDPhi.push_back(probe_obj.var("hltEgammaBestGsfTrackVarsUnseeded_DetaSeed",0))
                            #eg_bestTrkMissHits.push_back(probe_obj.var("hltEgammaBestGsfTrackVarsUnseeded_Dphi",0))
                            #eg_bestTrkNrLayerIT.push_back(probe_obj.var("hltEgammaBestGsfTrackVarsUnseeded_MissingHits",0))
                            eg_bestTrkESeedInvP.push_back(probe_obj.var("hltEgammaBestGsfTrackVarsUnseeded_NLayerIT",0))
                            eg_bestTrkInvEInvP.push_back(probe_obj.var("hltEgammaBestGsfTrackVarsUnseeded_OneOESeedMinusOneOP",0))
                            #eg_bestTrkValitHits.push_back(probe_obj.var("hltEgammaBestGsfTrackVarsUnseeded_OneOESuperMinusOneOP",0))
                            eg_hgcaliso_layerclus.push_back(probe_obj.var("hltEgammaBestGsfTrackVarsUnseeded_ValidHits",0))
                            eg_hgcaliso_layerclusem.push_back(probe_obj.var("hltEgammaHGCalLayerClusterIsoUnseeded",0))
                            eg_hgcaliso_layerclushad.push_back(probe_obj.var("hltEgammaHGCalLayerClusterIsoUnseeded_em",0))
                            
                            # Only process first valid pair per event (can be modified to process all pairs)
                            break
                        
                        if probe_count > 0:
                            break  # Found a valid pair, move to next tag
                    
                    # Store collection info
                    collection_name.push_back(col_name)
                    nr_objects.push_back(probe_count)
                    
                    if probe_count > 0:
                        if verbose:
                            print(f"    → Stored {probe_count} probe electron(s) from tag-probe pairs")
                        # Fill tree for this collection
                        tree.Fill()
                        break  # Found a working collection, stop trying others
                    else:
                        if verbose:
                            print(f"    → No valid tag-probe pairs found")
                                                                                    
                else:
                    print(f"  ✗ Invalid handle for {label}")
                    
            except Exception as e:
                print(f"  ✗ Error accessing {label}: {e}")
                import traceback
                if verbose:
                    traceback.print_exc()
    
    # Write and close
    out_file.Write()
    out_file.Close()
    print(f"\n✅ Debug ntuple saved to: {output_file}")
    print("You can open this file in ROOT to inspect the variables:")
    print("Note: Only probe electrons from valid tag-probe pairs are stored.")

if __name__ == "__main__":
    main() 

