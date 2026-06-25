// Asymmetry generator for 2024 pp reference run L2 Residual JEC.
//
// Unified DATA/MC entry point.  Compile via CMake then run:
//   ./runAsymmetry <input.root> <output.root> [--mc]
//
// What this produces:
//   - One ConeHistograms per cone size (label from kConeLabels in 2024ppRef.h).
//     Each contains a 4D THnSparse (eta_probe, pT_avg, alpha, A) filled once
//     per valid dijet event, plus 9 control TH1Ds.
//   - A few loose event-level QA histograms (vz, filter, trigger).
//
// JSON_handler.h requires nlohmann/json — install with: sudo pacman -S nlohmann-json
// All other dependencies are now in include/.

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
#include <stdexcept>
#include <cstring>
#include <iostream>

static constexpr Int_t   kNRefMax   = 200;
static constexpr float   kVzCut     = 15.0f;    // |vz| < 15 cm
static constexpr float   kMinPt     = 10.0f;    // minimum corrected pT for inclusive hists
static constexpr float   kMaxAbsA   = 0.7f;     // |A| cut

int main(int argc, char* argv[]) {

    if (argc < 3) {
        std::cerr << "Usage: runAsymmetry <input.root> <output.root> [--mc]\n";
        return 1;
    }

    TString input  = argv[1];
    TString output = argv[2];
    bool    isMC   = false;
    for (int i = 3; i < argc; i++) {
        if (std::string(argv[i]) == "--mc") isMC = true;
    }

    std::cout << "Mode: " << (isMC ? "MC" : "DATA") << "\n";
    std::cout << "Input:  " << input  << "\n";
    std::cout << "Output: " << output << "\n";

    // ---- jet energy corrections ----
    JetCorrector jec(kJECFiles);

    // ---- jet ID + veto map ----
    // Applied to leading and subleading jets; events are rejected if either fails.
    // Applied to third jet to set hasThird (events are NOT rejected if third fails).
    JetSelect js(kVetoMapPath);

    // ---- golden JSON (DATA only) ----
    // Requires nlohmann-json: sudo pacman -S nlohmann-json
    // JSON_handler dcs(kJSONPath);

    // ---- TTree branch structs ----
    JetStruct<kNRefMax> jets;
    EventStruct         event;
    FiltersStruct       filters;
    Int_t hlt_j40 = 0, hlt_j80 = 0;   // trigger decisions (DATA only)

    // ---- open input file and fetch TTrees ----
    TFile* fi = TFile::Open(input, "read");
    if (!fi || fi->IsZombie()) {
        std::cerr << "Cannot open " << input << "\n";
        return 1;
    }

    TTree* tEvt  = (TTree*)fi->Get(kHiTreePath);
    TTree* tJet  = (TTree*)fi->Get(kJetTreePath);
    TTree* tSkim = isMC ? nullptr : (TTree*)fi->Get(kSkimTreePath);
    TTree* tTrig = isMC ? nullptr : (TTree*)fi->Get(kTrigTreePath);

    if (!tEvt || !tJet) {
        std::cerr << "Missing HiTree or jet tree in " << input << "\n";
        return 1;
    }
    if (!isMC && (!tSkim || !tTrig)) {
        std::cerr << "Missing skim or HLT tree in " << input
                  << "\n(check kSkimTreePath and kTrigTreePath in configs/2024ppRef.h)\n";
        return 1;
    }

    // bind branches — SetBranches disables all others first
    SetBranches(tEvt, event.BranchMap(isMC));
    SetBranches(tJet, jets.BranchMap(isMC));
    if (!isMC) {
        SetBranches(tSkim, filters.BranchMap());

        tTrig->SetBranchStatus("*", 0);
        tTrig->SetBranchStatus(kHLTJ40Branch, 1);
        tTrig->SetBranchStatus(kHLTJ80Branch, 1);
        tTrig->SetBranchAddress(kHLTJ40Branch, &hlt_j40);
        tTrig->SetBranchAddress(kHLTJ80Branch, &hlt_j80);
    }

    // ---- QA histograms (event-level, not per-cone) ----
    TH1D* hvz_all = new TH1D("hvz_all", "all events;v_{z} (cm);N",         40, -20, 20);
    TH1D* hvz     = new TH1D("hvz",     "after vz+filter;v_{z} (cm);N",    40, -20, 20);
    TH1I* hfilt   = isMC ? nullptr : new TH1I("hfilt",    "ppvF;filter;N",         2, 0, 2);
    TH1I* h_j40   = isMC ? nullptr : new TH1I("h_hlt_j40","HLT_AK4PFJet40;bit;N", 2, 0, 2);
    TH1I* h_j80   = isMC ? nullptr : new TH1I("h_hlt_j80","HLT_AK4PFJet80;bit;N", 2, 0, 2);

    // ---- physics histograms (one ConeHistograms per cone size) ----
    BinningConfig bins;
    std::vector<ConeHistograms> cones(kConeLabels.size());
    for (size_t c = 0; c < kConeLabels.size(); c++) {
        cones[c].Init(kConeLabels[c], bins);
    }
    // TODO: when adding more cone sizes, extend kConeLabels in configs/2024ppRef.h
    // and change the jet tree path to match (e.g., "ak3PFJetAnalyzer/t")

    // ---- corrected pT buffer ----
    float corrPt[kNRefMax] = {};

    // ---- event loop ----
    const Long64_t nEvents = tEvt->GetEntries();
    ProgressBar pb(isMC ? "MC events:" : "DATA events:", (int)nEvents);

    for (Long64_t i = 0; i < nEvents; i++) {
        pb.Update();

        tEvt->GetEntry(i);
        hvz_all->Fill(event.vz);

        if (TMath::Abs(event.vz) > kVzCut) continue;

        // DATA-only: primary vertex filter
        if (!isMC) {
            tSkim->GetEntry(i);
            hfilt->Fill(filters.ppvF);
            if (filters.ppvF == 0) continue;
        }

        hvz->Fill(event.vz);

        // DATA-only: read trigger decisions now (cut applied below after sorting)
        if (!isMC) {
            tTrig->GetEntry(i);
            h_j40->Fill(hlt_j40);
            h_j80->Fill(hlt_j80);
        }

        // load jets
        tJet->GetEntry(i);
        if (jets.reco.nref < 2) continue;

        // DATA-only: golden JSON — uncomment after: sudo pacman -S nlohmann-json
        // if (!dcs.isGood(event.run, event.lumi)) continue;

        const float weight = event.w;   // 1.0 for DATA, pthat weight for MC

        // apply JEC to every jet and fill inclusive control histograms
        for (int j = 0; j < jets.reco.nref; j++) {
            jec.SetJetPT(jets.reco.rawpt[j]);
            jec.SetJetEta(jets.reco.eta[j]);
            jec.SetJetPhi(jets.reco.phi[j]);
            corrPt[j] = (float)jec.GetCorrectedPT();

            if (corrPt[j] < kMinPt) continue;
            cones[0].FillInclJet(corrPt[j], jets.reco.eta[j], jets.reco.phi[j], weight);
        }

        // rank jets by corrected pT (single-pass O(n))
        SortedJets sorted = FindLeadingJets(corrPt, jets.reco.nref);
        if (sorted.lead == -1 || sorted.sublead == -1) continue;

        // DATA-only: trigger cut — leading jet must be in the fully efficient region of
        // at least one fired trigger.  Mirrors original logic exactly: if a trigger fires
        // but the leading jet pT is below its threshold the event is rejected even if
        // another trigger also fired and is efficient.
        if (!isMC) {
            bool noTrig       = (hlt_j40 == 0 && hlt_j80 == 0);
            bool j40Ineffic   = (hlt_j40 == 1 && corrPt[sorted.lead] <= kHLTJ40Thresh);
            bool j80Ineffic   = (hlt_j80 == 1 && corrPt[sorted.lead] <= kHLTJ80Thresh);
            if (noTrig || j40Ineffic || j80Ineffic) continue;
        }

        // jet ID + veto map on leading and subleading — reject event if either fails
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

        // jet ID on third jet — failure sets hasThird=false and alpha=0, does not reject event.
        // This is stricter than the original (which used the third jet pT unconditionally):
        // a fake/bad third jet won't incorrectly suppress the event in tight alpha slices.
        bool hasThird = (sorted.third != -1);
        if (hasThird) {
            hasThird = js.JetSelection(jets.reco.eta[sorted.third], jets.reco.phi[sorted.third],
                                       pf.CHF[sorted.third], pf.NHF[sorted.third], pf.CEF[sorted.third],
                                       pf.NEF[sorted.third], pf.MUF[sorted.third],
                                       pf.CHM[sorted.third], pf.NHM[sorted.third], pf.CEM[sorted.third],
                                       pf.NEM[sorted.third], pf.MUM[sorted.third]);
        }

        // build the dijet — returns invalid result if dphi/barrel/sublead-pT cuts fail
        DijetResult dijet = MakeDijet(sorted, hasThird,
                                      corrPt, jets.reco.eta, jets.reco.phi,
                                      event.event);
        if (!dijet.valid) continue;

        if (TMath::Abs(dijet.A) > kMaxAbsA) continue;

        // fill physics histograms for all cone sizes
        // (loop is trivial for now — extend for AK3PF, AK6PF etc.)
        for (size_t c = 0; c < cones.size(); c++) {
            cones[c].Fill(dijet, corrPt, jets.reco.eta, jets.reco.phi, weight);
        }
    }

    pb.Finish();

    // ---- write output ----
    TFile* fo = new TFile(output, "recreate");
    fo->cd();

    hvz_all->Write();
    hvz->Write();
    if (!isMC) {
        hfilt->Write();
        h_j40->Write();
        h_j80->Write();
    }
    for (auto& cone : cones) cone.Write();

    fo->Close();
    fi->Close();

    return 0;
}
