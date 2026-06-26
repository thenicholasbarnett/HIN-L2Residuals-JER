#include "RunAsymmetry.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"
#include "Rtypes.h"

#include "JetCorrector.h"
#include "JetSelection.h"
// #include "JSON_handler.h"   // requires: sudo pacman -S nlohmann-json

#include "BranchMapping.h"
#include "EventStructs.h"
#include "JetStruct.h"
#include "ProgressBar.h"
#include "Binning.h"
#include "Dijet.h"
#include "DijetHistograms.h"

#include "2024ppRef.h"

#include <vector>
#include <string>
#include <iostream>

enum class RunMode { MC, ZeroBias, HardProbes };

static constexpr Int_t  kNRefMax = 200;
static constexpr float  kVzCut   = 15.0f;
static constexpr float  kMinPt   = 10.0f;
static constexpr float  kMaxAbsA = 0.7f;

void runAsymmetry(TString input, TString output, TString modeFlag) {

    RunMode mode = RunMode::HardProbes;
    if      (modeFlag == "--mc")          mode = RunMode::MC;
    else if (modeFlag == "--zero-bias")   mode = RunMode::ZeroBias;
    else if (modeFlag == "--hard-probes") mode = RunMode::HardProbes;

    const char* modeStr = (mode == RunMode::MC)       ? "MC"
                        : (mode == RunMode::ZeroBias)  ? "zero bias"
                                                       : "hard probes";
    std::cout << "Mode:   " << modeStr << "\n";
    std::cout << "Input:  " << input   << "\n";
    std::cout << "Output: " << output  << "\n";

    JetCorrector jec(kJECFiles);
    JetSelect js(kVetoMapPath);
    // JSON_handler dcs(kJSONPath);   // requires: sudo pacman -S nlohmann-json

    // ---- TTree branch structs ----
    JetStruct<kNRefMax> jets;
    EventStruct         event;
    FiltersStruct       filters;
    Int_t               hlt_j80 = 0;

    // ---- open input file ----
    TFile* fi = TFile::Open(input, "read");
    if (!fi || fi->IsZombie()) {
        std::cerr << "Cannot open " << input << "\n";
        return;
    }

    // ---- TTree array: jet trees first, then event / skim / trig ----
    const size_t nJetTrees = kConeLabels.size();
    const size_t kEvtIdx   = nJetTrees;
    const size_t kSkimIdx  = nJetTrees + 1;
    const size_t kTrigIdx  = nJetTrees + 2;

    std::vector<TTree*> trees(nJetTrees + 3, nullptr);

    for (size_t c = 0; c < nJetTrees; c++) {
        trees[c] = (TTree*)fi->Get(kJetTreePaths[c]);
        if (!trees[c]) {
            std::cerr << "Missing jet tree " << kJetTreePaths[c] << " in " << input << "\n";
            return;
        }
    }
    trees[kEvtIdx] = (TTree*)fi->Get(kHiTreePath);
    if (!trees[kEvtIdx]) {
        std::cerr << "Missing HiTree in " << input << "\n";
        return;
    }
    if (mode != RunMode::MC) {
        trees[kSkimIdx] = (TTree*)fi->Get(kSkimTreePath);
        if (!trees[kSkimIdx]) {
            std::cerr << "Missing skim tree in " << input
                      << "\n(check kSkimTreePath in cfg/2024ppRef.h)\n";
            return;
        }
    }
    if (mode == RunMode::HardProbes) {
        trees[kTrigIdx] = (TTree*)fi->Get(kTrigTreePath);
        if (!trees[kTrigIdx]) {
            std::cerr << "Missing HLT tree in " << input
                      << "\n(check kTrigTreePath in cfg/2024ppRef.h)\n";
            return;
        }
    }

    // ---- branch setup ----
    const bool isMC = (mode == RunMode::MC);
    SetBranches(trees[kEvtIdx], event.BranchMap(isMC));
    SetBranches(trees[0],       jets.BranchMap(isMC));
    if (mode != RunMode::MC) {
        SetBranches(trees[kSkimIdx], filters.BranchMap());
    }
    if (mode == RunMode::HardProbes) {
        trees[kTrigIdx]->SetBranchStatus("*", 0);
        trees[kTrigIdx]->SetBranchStatus(kHLTJ80Branch, 1);
        trees[kTrigIdx]->SetBranchAddress(kHLTJ80Branch, &hlt_j80);
    }

    // ---- QA histograms ----
    TH1D* hvz_all = new TH1D("hvz_all", "all events;v_{z} (cm);N",         40, -20, 20);
    TH1D* hvz     = new TH1D("hvz",     "after vz+filter;v_{z} (cm);N",    40, -20, 20);
    TH1I* hfilt   = (mode != RunMode::MC)         ? new TH1I("hfilt",     "ppvF;filter;N",         2, 0, 2) : nullptr;
    TH1I* h_j80   = (mode == RunMode::HardProbes) ? new TH1I("h_hlt_j80", "HLT_AK4PFJet80;bit;N", 2, 0, 2) : nullptr;

    // ---- physics histograms ----
    BinningConfig bins;
    std::vector<ConeHistograms> cones(nJetTrees);
    for (size_t c = 0; c < nJetTrees; c++) cones[c].Init(kConeLabels[c], bins);

    // ---- corrected pT buffer ----
    float corrPt[kNRefMax] = {};

    // ---- event loop ----
    const Long64_t nEvents = trees[kEvtIdx]->GetEntries();
    ProgressBar pb(modeStr, (int)nEvents);

    for (Long64_t i = 0; i < nEvents; i++) {
        pb.Update();

        trees[kEvtIdx]->GetEntry(i);
        hvz_all->Fill(event.vz);
        if (TMath::Abs(event.vz) > kVzCut) continue;

        if (mode != RunMode::MC) {
            trees[kSkimIdx]->GetEntry(i);
            hfilt->Fill(filters.ppvF);
            if (filters.ppvF == 0) continue;
        }

        hvz->Fill(event.vz);

        if (mode == RunMode::HardProbes) {
            trees[kTrigIdx]->GetEntry(i);
            h_j80->Fill(hlt_j80);
        }

        for (size_t c = 0; c < nJetTrees; c++) trees[c]->GetEntry(i);
        if (jets.reco.nref < 2) continue;

        // DATA modes: golden JSON — uncomment after: sudo pacman -S nlohmann-json
        // if (mode != RunMode::MC && !dcs.isGood(event.run, event.lumi)) continue;

        const float weight = event.w;

        for (int j = 0; j < jets.reco.nref; j++) {
            jec.SetJetPT(jets.reco.rawpt[j]);
            jec.SetJetEta(jets.reco.eta[j]);
            jec.SetJetPhi(jets.reco.phi[j]);
            corrPt[j] = (float)jec.GetCorrectedPT();

            if (corrPt[j] < kMinPt) continue;
            cones[0].FillInclJet(corrPt[j], jets.reco.eta[j], jets.reco.phi[j], weight);
        }

        SortedJets sorted = FindLeadingJets(corrPt, jets.reco.nref);
        if (sorted.lead == -1 || sorted.sublead == -1) continue;

        if (mode == RunMode::HardProbes) {
            bool noTrig     = (hlt_j80 == 0);
            bool j80Ineffic = (hlt_j80 == 1 && corrPt[sorted.lead] <= kHLTJ80Thresh);
            if (noTrig || j80Ineffic) continue;
        }

        auto& pf = jets.reco.pf;
        if (!js.JetSelection(jets.reco.eta[sorted.lead],    jets.reco.phi[sorted.lead],
                             pf.CHF[sorted.lead], pf.NHF[sorted.lead], pf.CEF[sorted.lead],
                             pf.NEF[sorted.lead], pf.MUF[sorted.lead],
                             pf.CHM[sorted.lead], pf.NHM[sorted.lead], pf.CEM[sorted.lead],
                             pf.NEM[sorted.lead], pf.MUM[sorted.lead])) continue;
        if (!js.JetSelection(jets.reco.eta[sorted.sublead], jets.reco.phi[sorted.sublead],
                             pf.CHF[sorted.sublead], pf.NHF[sorted.sublead], pf.CEF[sorted.sublead],
                             pf.NEF[sorted.sublead], pf.MUF[sorted.sublead],
                             pf.CHM[sorted.sublead], pf.NHM[sorted.sublead], pf.CEM[sorted.sublead],
                             pf.NEM[sorted.sublead], pf.MUM[sorted.sublead])) continue;

        bool hasThird = (sorted.third != -1);
        if (hasThird) {
            hasThird = js.JetSelection(jets.reco.eta[sorted.third], jets.reco.phi[sorted.third],
                                       pf.CHF[sorted.third], pf.NHF[sorted.third], pf.CEF[sorted.third],
                                       pf.NEF[sorted.third], pf.MUF[sorted.third],
                                       pf.CHM[sorted.third], pf.NHM[sorted.third], pf.CEM[sorted.third],
                                       pf.NEM[sorted.third], pf.MUM[sorted.third]);
        }

        DijetResult dijet = MakeDijet(sorted, hasThird,
                                      corrPt, jets.reco.eta, jets.reco.phi,
                                      event.event);
        if (!dijet.valid) continue;
        if (TMath::Abs(dijet.A) > kMaxAbsA) continue;

        for (size_t c = 0; c < cones.size(); c++)
            cones[c].Fill(dijet, corrPt, jets.reco.eta, jets.reco.phi, weight);
    }

    pb.Finish();

    // ---- write output ----
    TFile* fo = new TFile(output, "recreate");
    fo->cd();
    hvz_all->Write();
    hvz->Write();
    if (hfilt) hfilt->Write();
    if (h_j80) h_j80->Write();
    for (auto& cone : cones) cone.Write();
    fo->Close();
    fi->Close();
}
