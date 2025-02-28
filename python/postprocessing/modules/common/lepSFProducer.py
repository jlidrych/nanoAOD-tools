from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True


class lepSFProducer(Module):
    def __init__(self, muonSelectionTag, electronSelectionTag):
        if muonSelectionTag == "LooseWP_2016":
            mu_f = ["Mu_Trg.root", "Mu_ID.root", "Mu_Iso.root"]
            mu_h = [
                "IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio",
                "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio",
                "LooseISO_LooseID_pt_eta/pt_abseta_ratio"
            ]
        elif muonSelectionTag == "MediumWP_UL2016_preVPF":
            mu_f = ["Run2016_UL_HIPM_SingleMuonTriggers.root","Run2016_UL_HIPM_ID.root","Run2016_UL_HIPM_ISO.root"]
            mu_h = [
                "NUM_IsoMu24_or_IsoTkMu24_or_Mu50_or_TkMu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_MediumID_DEN_TrackerMuons_abseta_pt",
                "NUM_TightRelIso_DEN_MediumID_abseta_pt"
            ]
        elif muonSelectionTag == "MediumWP_UL2016":
            mu_f = ["Run2016_UL_SingleMuonTriggers.root","Run2016_UL_ID.root","Run2016_UL_ISO.root"]
            mu_h = [
                "NUM_IsoMu24_or_IsoTkMu24_or_Mu50_or_TkMu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_MediumID_DEN_TrackerMuons_abseta_pt",
                "NUM_TightRelIso_DEN_MediumID_abseta_pt"
            ]
        elif muonSelectionTag == "MediumWP_UL2017":
            mu_f = ["Run2017_UL_SingleMuonTriggers.root","Run2017_UL_ID.root","Run2017_UL_ISO.root"]
            mu_h = [
                "NUM_IsoMu27_or_Mu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_MediumID_DEN_TrackerMuons_abseta_pt",
                "NUM_TightRelIso_DEN_MediumID_abseta_pt"
            ]
        elif muonSelectionTag == "MediumWP_UL2018":
            mu_f = ["Run2018_UL_SingleMuonTriggers.root","Run2018_UL_ID.root","Run2018_UL_ISO.root"]
            mu_h = [
                "NUM_IsoMu24_or_Mu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_MediumID_DEN_TrackerMuons_abseta_pt",
                "NUM_TightRelIso_DEN_MediumID_abseta_pt"
            ]
        elif muonSelectionTag == "HighPt_UL2017":
            mu_f = ["Run2017_UL_ID.root","Run2017_UL_ISO.root"]
#            mu_f = ["Run2017_UL_SingleMuonTriggers.root","Run2017_UL_ID.root","Run2017_UL_ISO.root"]
            mu_h = [
#                "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt",
                "NUM_TrkHighPtID_DEN_TrackerMuons_abseta_pt",
                "NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut_abseta_pt"
            ]
        elif muonSelectionTag == "LooseWP_UL2016_preVPF":
            mu_f = ["Run2016_UL_HIPM_SingleMuonTriggers.root","Run2016_UL_HIPM_ID.root","Run2016_UL_HIPM_ISO.root"]
            mu_h = [
                "NUM_IsoMu24_or_IsoTkMu24_or_Mu50_or_TkMu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                "NUM_LooseRelIso_DEN_LooseID_abseta_pt"
            ]
        elif muonSelectionTag == "LooseWP_UL2016":
            mu_f = ["Run2016_UL_SingleMuonTriggers.root","Run2016_UL_ID.root","Run2016_UL_ISO.root"]
            mu_h = [
                "NUM_IsoMu24_or_IsoTkMu24_or_Mu50_or_TkMu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                "NUM_LooseRelIso_DEN_LooseID_abseta_pt"
            ]
        elif muonSelectionTag == "LooseWP_UL2017":
            mu_f = ["Run2017_UL_SingleMuonTriggers.root","Run2017_UL_ID.root","Run2017_UL_ISO.root"]
            mu_h = [
                "NUM_IsoMu27_or_Mu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                "NUM_LooseRelIso_DEN_LooseID_abseta_pt"
            ]
        elif muonSelectionTag == "LooseWP_UL2018":
            mu_f = ["Run2018_UL_SingleMuonTriggers.root","Run2018_UL_ID.root","Run2018_UL_ISO.root"]
            mu_h = [
                "NUM_IsoMu24_or_Mu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                "NUM_LooseRelIso_DEN_LooseID_abseta_pt"
            ]
        if electronSelectionTag == "GPMVA90_2016":
            el_f = ["EGM2D_eleGSF.root", "EGM2D_eleMVA90.root"]
            el_h = ["EGamma_SF2D", "EGamma_SF2D"]
        elif electronSelectionTag == "MediumWP_UL2016_preVPF":
            el_f = ["Run2016_UL_HIPM_SingleMuonTriggers.root","Run2016_UL_HIPM_ID.root","Run2016_UL_HIPM_ISO.root"]
            el_h = [
                "NUM_IsoMu24_or_IsoTkMu24_or_Mu50_or_TkMu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_MediumID_DEN_TrackerMuons_abseta_pt",
                "NUM_TightRelIso_DEN_MediumID_abseta_pt"
            ]
        elif electronSelectionTag == "MediumWP_UL2016":
            el_f = ["Run2016_UL_SingleMuonTriggers.root","Run2016_UL_ID.root","Run2016_UL_ISO.root"]
            el_h = [
                "NUM_IsoMu24_or_IsoTkMu24_or_Mu50_or_TkMu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_MediumID_DEN_TrackerMuons_abseta_pt",
                "NUM_TightRelIso_DEN_MediumID_abseta_pt"
            ]
        elif electronSelectionTag == "MediumWP_UL2017":
            el_f = ["Run2017_UL_SingleMuonTriggers.root","Run2017_UL_ID.root","Run2017_UL_ISO.root"]
            el_h = [
                "NUM_IsoMu27_or_Mu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_MediumID_DEN_TrackerMuons_abseta_pt",
                "NUM_TightRelIso_DEN_MediumID_abseta_pt"
            ]
        elif electronSelectionTag == "MediumWP_UL2018":
            el_f = ["Run2018_UL_SingleMuonTriggers.root","Run2018_UL_ID.root","Run2018_UL_ISO.root"]
            el_h = [
                "NUM_IsoMu24_or_Mu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_MediumID_DEN_TrackerMuons_abseta_pt",
                "NUM_TightRelIso_DEN_MediumID_abseta_pt"
            ]
        elif electronSelectionTag == "LooseWP_UL2016_preVPF":
            el_f = ["Run2016_UL_HIPM_SingleMuonTriggers.root","Run2016_UL_HIPM_ID.root","Run2016_UL_HIPM_ISO.root"]
            el_h = [
                "NUM_IsoMu24_or_IsoTkMu24_or_Mu50_or_TkMu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                "NUM_LooseRelIso_DEN_LooseID_abseta_pt"
            ]
        elif electronSelectionTag == "LooseWP_UL2016":
            el_f = ["Run2016_UL_SingleMuonTriggers.root","Run2016_UL_ID.root","Run2016_UL_ISO.root"]
            el_h = [
                "NUM_IsoMu24_or_IsoTkMu24_or_Mu50_or_TkMu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                "NUM_LooseRelIso_DEN_LooseID_abseta_pt"
            ]
        elif electronSelectionTag == "LooseWP_UL2017":
            el_f = ["Run2017_UL_SingleMuonTriggers.root","Run2017_UL_ID.root","Run2017_UL_ISO.root"]
            el_h = [
                "NUM_IsoMu27_or_Mu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                "NUM_LooseRelIso_DEN_LooseID_abseta_pt"
            ]
        elif electronSelectionTag == "LooseWP_UL2018":
            el_f = ["Run2018_UL_SingleMuonTriggers.root","Run2018_UL_ID.root","Run2018_UL_ISO.root"]
            el_h = [
                "NUM_IsoMu24_or_Mu50_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",
                "NUM_LooseID_DEN_TrackerMuons_abseta_pt",
                "NUM_LooseRelIso_DEN_LooseID_abseta_pt"
            ]
        mu_f = [
            "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/leptonSF/"
            % os.environ['CMSSW_BASE'] + f for f in mu_f
        ]
        el_f = [
            "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/leptonSF/"
            % os.environ['CMSSW_BASE'] + f for f in el_f
        ]

        self.mu_f = ROOT.std.vector(str)(len(mu_f))
        self.mu_h = ROOT.std.vector(str)(len(mu_f))
        for i in range(len(mu_f)):
            self.mu_f[i] = mu_f[i]
            self.mu_h[i] = mu_h[i]
        self.el_f = ROOT.std.vector(str)(len(el_f))
        self.el_h = ROOT.std.vector(str)(len(el_f))
        for i in range(len(el_f)):
            self.el_f[i] = el_f[i]
            self.el_h[i] = el_h[i]

        if "/LeptonEfficiencyCorrector_cc.so" not in ROOT.gSystem.GetLibraries(
        ):
            print("Load C++ Worker")
            ROOT.gROOT.ProcessLine(
                ".L %s/src/PhysicsTools/NanoAODTools/python/postprocessing/helpers/LeptonEfficiencyCorrector.cc+"
                % os.environ['CMSSW_BASE'])

    def beginJob(self):
        self._worker_mu = ROOT.LeptonEfficiencyCorrector(self.mu_f, self.mu_h)
        self._worker_el = ROOT.LeptonEfficiencyCorrector(self.el_f, self.el_h)

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Muon_SF"       , "F", lenVar="nMuon"    )
        self.out.branch("Muon_SFErr"    , "F", lenVar="nMuon"    )
        self.out.branch("Electron_SF"   , "F", lenVar="nElectron")
        self.out.branch("Electron_SFErr", "F", lenVar="nElectron")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")

        sf_mu = [ self._worker_mu.getSF(mu.pdgId,mu.pt,mu.eta) for mu in muons ]
        sferr_mu = [ self._worker_mu.getSFErr(mu.pdgId,mu.pt,mu.eta) for mu in muons ]
        sf_el = [ self._worker_el.getSF(el.pdgId,el.pt,el.eta) for el in electrons ]
        sferr_el = [ self._worker_el.getSFErr(el.pdgId,el.pt,el.eta) for el in electrons ]

        self.out.fillBranch("Muon_SF"       , sf_mu)
        self.out.fillBranch("Muon_SFErr"    , sferr_mu)
        self.out.fillBranch("Electron_SF"   , sf_el)
        self.out.fillBranch("Electron_SFErr", sferr_el)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid
# having them loaded when not needed

lepSF = lambda: lepSFProducer("LooseWP_2016", "GPMVA90_2016")
lepSF_UL2016_preVPF_WPL = lambda : lepSFProducer( "LooseWP_UL2016_preVPF", "LooseWP_UL2016_preVPF")
lepSF_UL2016_WPL = lambda : lepSFProducer( "LooseWP_UL2016", "LooseWP_UL2016")
lepSF_UL2017_WPL = lambda : lepSFProducer( "LooseWP_UL2017", "LooseWP_UL2017")
lepSF_UL2018_WPL = lambda : lepSFProducer( "LooseWP_UL2018", "LooseWP_UL2018")
lepSF_UL2016_preVPF_WPM = lambda : lepSFProducer( "MediumWP_UL2016_preVPF","MediumWP_UL2016_preVPF")
lepSF_UL2016_WPM = lambda : lepSFProducer( "MediumWP_UL2016","MediumWP_UL2016")
lepSF_UL2017_WPM = lambda : lepSFProducer( "MediumWP_UL2017","MediumWP_UL2017")
lepSF_UL2018_WPM = lambda : lepSFProducer( "MediumWP_UL2018","MediumWP_UL2018")
lepSF_UL2017_HighPt = lambda : lepSFProducer("HighPt_UL2017","MediumWP_UL2017")
