#!/usr/bin/env python3

import ROOT
from DataFormats.FWLite import Events, Handle
import argparse
import glob
import os


def get_genparts(genparts, pid=11, antipart=True, status=1):
    selected = []
    if genparts is None:
        return selected

    for part in genparts:
        pdg_id = part.pdgId()
        if pdg_id == pid or (antipart and abs(pdg_id) == abs(pid)):
            if part.isHardProcess():
                if status == 1:
                    selected.append(part)
    return selected


def match_to_gen(eta, phi, genparts, pid=11, antipart=True, max_dr=0.1, status=1):
    best_match = None
    best_dr2 = max_dr * max_dr
    best_pt = -1.0
    selected_parts = get_genparts(genparts, pid, antipart, status)
    for part in selected_parts:
        dr2 = ROOT.reco.deltaR2(eta, phi, part.eta(), part.phi())
        if dr2 < best_dr2:
            best_match = part
            best_dr2 = dr2
            best_pt = part.pt()
    return best_match, best_dr2, best_pt


def main():
    parser = argparse.ArgumentParser(
        description="Create debug ntuples for E/Gamma trigger objects from Phase2 data with Gen matching (pt-filtered)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "input_file",
        nargs="?",
        help="Input ROOT file containing E/Gamma trigger objects",
    )
    input_group.add_argument(
        "--input-dir",
        "-d",
        help="Directory containing ROOT files to process",
    )
    input_group.add_argument(
        "--input-pattern",
        "-p",
        help='Pattern to match ROOT files (e.g., "*.root")',
    )

    parser.add_argument("--output", "-o", required=True, help="Output ROOT file")
    parser.add_argument("--max-events", "-n", type=int, default=10, help="Maximum events")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose")
    parser.add_argument("--gen-label", default="genParticles", help="Gen particle label")
    parser.add_argument("--gen-process", default="", help="Gen process name")
    parser.add_argument("--max-dr", type=float, default=0.1, help="Max dR for gen matching")
    parser.add_argument("--gen-pid", type=int, default=11, help="PDG ID for gen matching")
    # New: pt-range filtering (in GeV). Objects with et outside this range will be skipped.
    parser.add_argument("--pt-min", type=float, default=0.0, help="Minimum pt/et to store (GeV)")
    parser.add_argument("--pt-max", type=float, default=500.0, help="Maximum pt/et to store (GeV)")

    args = parser.parse_args()

    input_files = []
    if args.input_file:
        input_files = [args.input_file]
    elif args.input_dir:
        if "eos" in args.input_dir:
            import subprocess

            try:
                result = subprocess.run(
                    f"ls {args.input_dir}/*.root",
                    capture_output=True,
                    text=True,
                    shell=True,
                )
                if result.returncode == 0:
                    input_files = result.stdout.strip().split("\n")
                else:
                    input_files = glob.glob(os.path.join(args.input_dir, "*.root"))
            except Exception:
                input_files = glob.glob(os.path.join(args.input_dir, "*.root"))
        else:
            input_files = glob.glob(os.path.join(args.input_dir, "*.root"))
    elif args.input_pattern:
        if "eos" in args.input_pattern:
            import subprocess

            try:
                result = subprocess.run(
                    f"ls {args.input_pattern}",
                    capture_output=True,
                    text=True,
                    shell=True,
                )
                if result.returncode == 0:
                    input_files = result.stdout.strip().split("\n")
                else:
                    input_files = glob.glob(args.input_pattern)
            except Exception:
                input_files = glob.glob(args.input_pattern)
        else:
            input_files = glob.glob(args.input_pattern)

    if not input_files:
        print("❌ No input files found!")
        return

    input_files.sort()

    output_file = args.output
    max_events = args.max_events
    verbose = args.verbose
    gen_label = args.gen_label
    gen_process = args.gen_process
    max_dr = args.max_dr
    gen_pid = args.gen_pid
    pt_min = args.pt_min
    pt_max = args.pt_max

    print(f"📁 Found {len(input_files)} input files")
    print(f"💾 Output: {output_file}")
    print(f"⚡ Max events: {max_events}, pt-range: [{pt_min}, {pt_max}] GeV")

    egtrigobjs_handle = Handle("std::vector<trigger::EgammaObject>")
    egtrigobjs_unseeded_handle = Handle("std::vector<trigger::EgammaObject>")
    genparts_handle = Handle("std::vector<reco::GenParticle>")

    out_file = ROOT.TFile(output_file, "RECREATE")
    tree = ROOT.TTree("egHLTTree", "EGamma Tree with Gen Matching (pt-filtered)")

    # Define vectors (same as original)
    run = ROOT.std.vector("int")()
    lumi = ROOT.std.vector("int")()
    event = ROOT.std.vector("int")()
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
    eg_ecalPFIsol_default = ROOT.std.vector("float")()
    eg_hcalPFIsol_default = ROOT.std.vector("float")()
    eg_hgcalPFIsol_default = ROOT.std.vector("float")()
    eg_trkIsolV0_default = ROOT.std.vector("float")()
    eg_trkIsolV6_default = ROOT.std.vector("float")()
    eg_trkIsolV72_default = ROOT.std.vector("float")()
    eg_trkChi2_default = ROOT.std.vector("float")()
    eg_trkMissHits = ROOT.std.vector("int")()
    eg_trkValidHits = ROOT.std.vector("int")()
    eg_invESeedInvP = ROOT.std.vector("float")()
    eg_invEInvP = ROOT.std.vector("float")()
    eg_trkDEta = ROOT.std.vector("float")()
    eg_trkDEtaSeed = ROOT.std.vector("float")()
    eg_trkDPhi = ROOT.std.vector("float")()
    eg_trkNrLayerIT = ROOT.std.vector("int")()
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
    eg_hgcaliso_layerclus = ROOT.std.vector("float")()
    eg_hgcaliso_layerclusem = ROOT.std.vector("float")()
    eg_hgcaliso_layerclushad = ROOT.std.vector("float")()
    eg_ecalPFIsol = ROOT.std.vector("float")()
    eg_hcalPFIsol = ROOT.std.vector("float")()
    eg_trkIsolV6 = ROOT.std.vector("float")()
    gen_matched = ROOT.std.vector("int")()
    gen_pt = ROOT.std.vector("float")()
    gen_eta = ROOT.std.vector("float")()
    gen_phi = ROOT.std.vector("float")()
    gen_energy = ROOT.std.vector("float")()
    gen_pdgId = ROOT.std.vector("int")()
    gen_dr = ROOT.std.vector("float")()
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
    tree.Branch("eg_ecalPFIsol_default", eg_ecalPFIsol_default)
    tree.Branch("eg_hcalPFIsol_default", eg_hcalPFIsol_default)
    tree.Branch("eg_hgcalPFIsol_default", eg_hgcalPFIsol_default)
    tree.Branch("eg_trkIsolV0_default", eg_trkIsolV0_default)
    tree.Branch("eg_trkIsolV6_default", eg_trkIsolV6_default)
    tree.Branch("eg_trkIsolV72_default", eg_trkIsolV72_default)
    tree.Branch("eg_trkChi2_default", eg_trkChi2_default)
    tree.Branch("eg_trkMissHits", eg_trkMissHits)
    tree.Branch("eg_trkValidHits", eg_trkValidHits)
    tree.Branch("eg_invESeedInvP", eg_invESeedInvP)
    tree.Branch("eg_invEInvP", eg_invEInvP)
    tree.Branch("eg_trkDEta", eg_trkDEta)
    tree.Branch("eg_trkDEtaSeed", eg_trkDEtaSeed)
    tree.Branch("eg_trkDPhi", eg_trkDPhi)
    tree.Branch("eg_trkNrLayerIT", eg_trkNrLayerIT)
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
    tree.Branch("eg_hgcaliso_layerclus", eg_hgcaliso_layerclus)
    tree.Branch("eg_hgcaliso_layerclusem", eg_hgcaliso_layerclusem)
    tree.Branch("eg_hgcaliso_layerclushad", eg_hgcaliso_layerclushad)
    tree.Branch("eg_ecalPFIsol", eg_ecalPFIsol)
    tree.Branch("eg_hcalPFIsol", eg_hcalPFIsol)
    tree.Branch("eg_trkIsolV6", eg_trkIsolV6)
    tree.Branch("gen_matched", gen_matched)
    tree.Branch("gen_pt", gen_pt)
    tree.Branch("gen_eta", gen_eta)
    tree.Branch("gen_phi", gen_phi)
    tree.Branch("gen_energy", gen_energy)
    tree.Branch("gen_pdgId", gen_pdgId)
    tree.Branch("gen_dr", gen_dr)
    tree.Branch("collection_name", collection_name)
    tree.Branch("nr_objects", nr_objects)

    events = Events(input_files)

    print("Processing events...")
    for i, evt in enumerate(events):
        if i >= max_events:
            break

        if verbose:
            print(f"\n=== Event {i} ===")

        # Clear vectors
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
        eg_ecalPFIsol_default.clear()
        eg_hcalPFIsol_default.clear()
        eg_hgcalPFIsol_default.clear()
        eg_trkIsolV0_default.clear()
        eg_trkIsolV6_default.clear()
        eg_trkIsolV72_default.clear()
        eg_trkChi2_default.clear()
        eg_trkMissHits.clear()
        eg_trkValidHits.clear()
        eg_invESeedInvP.clear()
        eg_invEInvP.clear()
        eg_trkDEta.clear()
        eg_trkDEtaSeed.clear()
        eg_trkDPhi.clear()
        eg_trkNrLayerIT.clear()
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
        eg_hgcaliso_layerclus.clear()
        eg_hgcaliso_layerclusem.clear()
        eg_hgcaliso_layerclushad.clear()
        eg_ecalPFIsol.clear()
        eg_hcalPFIsol.clear()
        eg_trkIsolV6.clear()
        gen_matched.clear()
        gen_pt.clear()
        gen_eta.clear()
        gen_phi.clear()
        gen_energy.clear()
        gen_pdgId.clear()
        gen_dr.clear()
        collection_name.clear()
        nr_objects.clear()

        # Event info
        run.push_back(evt.eventAuxiliary().run())
        lumi.push_back(evt.eventAuxiliary().luminosityBlock())
        event.push_back(evt.eventAuxiliary().event())

        # Get gen particles
        genparts = None
        try:
            if gen_process:
                evt.getByLabel(gen_label, gen_process, genparts_handle)
            else:
                evt.getByLabel(gen_label, genparts_handle)
            if genparts_handle.isValid():
                genparts = genparts_handle.product()
                if verbose:
                    print(f"  ✓ Found {genparts.size()} gen particles")
            else:
                if verbose:
                    print(f"  ⚠ Gen particle collection '{gen_label}' not found")
        except Exception as e:
            if verbose:
                print(f"  ⚠ Error accessing gen particles: {e}")

        collections_to_try = [
            ("egtrigobjs_unseeded", "hltEgammaHLTExtra:Unseeded"),
            ("egtrigobjs_unseeded", "hltEgammaHLTExtra:Unseeded", "HLTX"),
        ]

        for collection_info in collections_to_try:
            if len(collection_info) == 2:
                col_name, label = collection_info
                process_name = ""
            else:
                col_name, label, process_name = collection_info

            if verbose:
                print(f"\nTrying collection: {col_name} label: {label}")

            try:
                handle = egtrigobjs_handle if col_name == "egtrigobjs" else egtrigobjs_unseeded_handle
                if process_name:
                    evt.getByLabel(label, process_name, handle)
                else:
                    evt.getByLabel(label, handle)

                if handle.isValid():
                    egobjs = handle.product()
                    if verbose:
                        print(f"  ✓ Found {egobjs.size()} objects")

                    if egobjs.size() > 0:
                        # We'll count how many objects we actually store after pt-filtering
                        stored_count = 0
                        for j, obj in enumerate(egobjs):
                            obj_et = obj.et()
                            # Apply pt/et filter (only store objects within [pt_min, pt_max])
                            if obj_et < pt_min or obj_et > pt_max:
                                continue

                            stored_count += 1
                            # Basic properties
                            eg_et.push_back(obj_et)
                            eg_energy.push_back(obj.energy())
                            obj_eta = obj.eta()
                            obj_phi = obj.phi()
                            eg_eta.push_back(obj_eta)
                            eg_phi.push_back(obj_phi)

                            # SuperCluster properties
                            eg_nrClus.push_back(obj.superCluster().clusters().size())
                            eg_rawEnergy.push_back(obj.superCluster().rawEnergy())
                            eg_phiWidth.push_back(obj.superCluster().phiWidth())
                            eg_seedId.push_back(obj.superCluster().seed().seed().rawId())
                            eg_seedDet.push_back(obj.superCluster().seed().seed().det())

                            # HLT variables (use existing keys)
                            eg_sigmaIEtaIEta.push_back(
                                obj.var("hltEgammaClusterShapeUnseeded_sigmaIEtaIEta5x5", 0)
                            )
                            eg_ecalPFIsol_default.push_back(obj.var("hltEgammaEcalPFClusterIsoUnseeded", 0))
                            eg_hcalPFIsol_default.push_back(obj.var("hltEgammaHcalPFClusterIsoUnseeded", 0))
                            eg_hgcalPFIsol_default.push_back(obj.var("hltEgammaHGCalPFClusterIsoUnseeded", 0))
                            eg_trkIsolV0_default.push_back(obj.var("hltEgammaEleGsfTrackIsoUnseeded", 0))
                            eg_trkIsolV6_default.push_back(obj.var("hltEgammaEleGsfTrackIsoV6Unseeded", 0))
                            eg_trkIsolV72_default.push_back(obj.var("hltEgammaEleGsfTrackIsoV72Unseeded", 0))
                            eg_trkChi2_default.push_back(obj.var("hltEgammaGsfTrackVarsUnseeded_Chi2", 0))
                            eg_invESeedInvP.push_back(obj.var("hltEgammaGsfTrackVarsUnseeded_OneOESeedMinusOneOP", 0))
                            eg_invEInvP.push_back(obj.var("hltEgammaGsfTrackVarsUnseeded_OneOESuperMinusOneOP", 0))
                            eg_trkDEta.push_back(obj.var("hltEgammaGsfTrackVarsUnseeded_Deta", 0))
                            eg_sigma2uu.push_back(obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2uu", 0))
                            eg_sigma2vv.push_back(obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2vv", 0))
                            eg_sigma2ww.push_back(obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2ww", 0))
                            eg_sigma2xx.push_back(obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2xx", 0))
                            eg_sigma2xy.push_back(obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2xy", 0))
                            eg_sigma2yy.push_back(obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2yy", 0))
                            eg_sigma2yz.push_back(obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2yz", 0))
                            eg_sigma2zx.push_back(obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2zx", 0))
                            eg_sigma2zz.push_back(obj.var("hltEgammaHGCALIDVarsUnseeded_sigma2zz", 0))
                            eg_pms2_default.push_back(obj.var("hltEgammaPixelMatchVarsUnseeded_s2", 0))
                            eg_hcalHForHoverE.push_back(obj.var("hltEgammaHGCALIDVarsUnseeded_hForHOverE", 0))
                            eg_l1TrkIsoCMSSW.push_back(obj.var("hltEgammaHoverEUnseeded", 0))
                            eg_bestTrkChi2.push_back(obj.var("hltEgammaEleL1TrkIsoUnseeded", 0))
                            eg_bestTrkDEta.push_back(obj.var("hltEgammaBestGsfTrackVarsUnseeded_Chi2", 0))
                            eg_bestTrkDEtaSeed.push_back(obj.var("hltEgammaBestGsfTrackVarsUnseeded_Deta", 0))
                            eg_bestTrkDPhi.push_back(obj.var("hltEgammaBestGsfTrackVarsUnseeded_DetaSeed", 0))
                            eg_bestTrkESeedInvP.push_back(obj.var("hltEgammaBestGsfTrackVarsUnseeded_NLayerIT", 0))
                            eg_bestTrkInvEInvP.push_back(obj.var("hltEgammaBestGsfTrackVarsUnseeded_OneOESeedMinusOneOP", 0))
                            eg_hgcaliso_layerclus.push_back(obj.var("hltEgammaBestGsfTrackVarsUnseeded_ValidHits", 0))
                            eg_hgcaliso_layerclusem.push_back(obj.var("hltEgammaHGCalLayerClusterIsoUnseeded", 0))
                            eg_hgcaliso_layerclushad.push_back(obj.var("hltEgammaHGCalLayerClusterIsoUnseeded_em", 0))

                            # Gen matching
                            if genparts is not None:
                                best_match, best_dr2, best_pt = match_to_gen(
                                    obj_eta, obj_phi, genparts, pid=gen_pid, antipart=True, max_dr=max_dr, status=1
                                )
                                if best_match is not None:
                                    this_dr = ROOT.reco.deltaR(obj_eta, obj_phi, best_match.eta(), best_match.phi())
                                    gen_matched.push_back(1)
                                    gen_pt.push_back(best_pt)
                                    gen_eta.push_back(best_match.eta())
                                    gen_phi.push_back(best_match.phi())
                                    gen_energy.push_back(best_match.energy())
                                    gen_pdgId.push_back(best_match.pdgId())
                                    gen_dr.push_back(this_dr)
                                else:
                                    gen_matched.push_back(0)
                                    gen_pt.push_back(-1.0)
                                    gen_eta.push_back(-999.0)
                                    gen_phi.push_back(-999.0)
                                    gen_energy.push_back(-1.0)
                                    gen_pdgId.push_back(0)
                                    gen_dr.push_back(-1.0)
                            else:
                                gen_matched.push_back(0)
                                gen_pt.push_back(-1.0)
                                gen_eta.push_back(-999.0)
                                gen_phi.push_back(-999.0)
                                gen_energy.push_back(-1.0)
                                gen_pdgId.push_back(0)
                                gen_dr.push_back(-1.0)

                        # After processing objects, store collection info and number stored
                        collection_name.push_back(col_name)
                        nr_objects.push_back(stored_count)
                        # Fill tree for this collection and stop trying other collections
                        tree.Fill()
                        break
                else:
                    if verbose:
                        print(f"  ✗ Invalid handle for {label}")
            except Exception as e:
                print(f"  ✗ Error accessing {label}: {e}")
                if verbose:
                    import traceback

                    traceback.print_exc()

    out_file.Write()
    out_file.Close()
    print(f"\n✅ Debug ntuple saved to: {output_file}")


if __name__ == "__main__":
    main()

