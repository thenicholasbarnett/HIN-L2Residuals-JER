#include "RunAsymmetry.h"

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"
#include "Rtypes.h"

#include "JetCorrector.h"
#include "JetSelection.h"
#include "JSON_handler.h"

#include "BranchMapping.h"
#include "EventStructs.h"
#include "JetStruct.h"
#include "Binning.h"
#include "Dijet.h"
#include "DijetHistograms.h"

#include "AnalysisConfig.h"

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

enum class RunMode { MC, ZeroBias, HardProbes };

static constexpr Int_t kNRefMax = 200;
static constexpr float kVzCut = 15.0f;

void runAsymmetry( TString input, TString output, TString modeFlag, Long64_t maxEvents ){

    const AnalysisConfig& cfg = Config();

    RunMode mode = RunMode::HardProbes;
    if( modeFlag == "--monte-carlo" || modeFlag == "-mc" ) mode = RunMode::MC;
    else if( modeFlag == "--zero-bias" || modeFlag == "-zb" ) mode = RunMode::ZeroBias;
    else if( modeFlag == "--hard-probes" || modeFlag == "-hp" ) mode = RunMode::HardProbes;
    else {
        throw std::invalid_argument( Form( "ERROR: invalid runAsymmetry mode '%s'. Expected --monte-carlo/-mc, --zero-bias/-zb, or --hard-probes/-hp.", modeFlag.Data() ) );
    }

    // checking configuration
    const size_t nCones = cfg.coneLabels.size();
    if( cfg.jecFilesPerCone.size() != nCones || cfg.jetTreePaths.size() != nCones ){
        std::cerr << "ERROR: cone labels, JEC files, and jet tree paths must be the same length\n";
        return;
    }

    // ?
    size_t trigConeIdx = nCones;
    for( size_t c = 0; c < nCones; c++ ){
        if( cfg.coneLabels[c] == cfg.trigCone ){ trigConeIdx = c; break; }
    }
    if( trigConeIdx == nCones ){
        std::cerr << "ERROR: trigger cone '" << cfg.trigCone << "' not found in cone labels\n";
        return;
    }

    // JetCorrector — chain per-cone L2Residual files onto the base JEC, but
    // only for data (HP/ZB); MC must never see residual corrections, so this
    // check is the only thing that decides that, not what's in the TOML.
    std::vector<JetCorrector> jecs;
    jecs.reserve( nCones );
    for( size_t c = 0; c < nCones; c++ ){
        std::vector<std::string> chain = cfg.jecFilesPerCone[c];
        if( mode != RunMode::MC && !cfg.residualFilesPerCone.empty() ){
            for( const auto& f : cfg.residualFilesPerCone[c] ){ chain.push_back( f ); }
        }
        jecs.emplace_back( chain );
    }

    // loading jet ID, jet veto map, and golden json
    JetSelect js( cfg.vetoMapPath, cfg.vetoMapHist );
    std::unique_ptr<JSON_handler> dcs;
    if( mode != RunMode::MC ){ dcs = std::make_unique<JSON_handler>( cfg.jsonPath.Data() ); }

    // structures
    std::vector<JetStruct<kNRefMax>> jets( nCones );
    EventStruct event;
    FiltersStruct filters;
    Int_t hlt_j80 = 0;

    // inpt
    TFile* fi = TFile::Open( input, "read" );
    if( !fi || fi->IsZombie() ){
        std::cerr << "Cannot open " << input << "\n";
        return;
    }

    // ttrees: jet collections, event info, filter, trigger
    const size_t kEvtIdx = nCones;
    const size_t kSkimIdx = nCones + 1;
    const size_t kTrigIdx = nCones + 2;
    std::vector<TTree*> trees( nCones + 3, nullptr );

    for( size_t c = 0; c < nCones; c++ ){
        trees[c] = ( TTree* )fi->Get( cfg.jetTreePaths[c] );
        if( !trees[c] ){
            std::cerr << "Missing jet tree " << cfg.jetTreePaths[c] << " in " << input << "\n";
            return;
        }
    }
    trees[kEvtIdx] = ( TTree* )fi->Get( cfg.hiTreePath );
    if( !trees[kEvtIdx] ){ std::cerr << "Missing HiTree in " << input << "\n"; return; }
    if( mode != RunMode::MC ){
        trees[kSkimIdx] = ( TTree* )fi->Get( cfg.skimTreePath );
        if( !trees[kSkimIdx] ){
            std::cerr << "Missing skim tree in " << input << "\n(check cfg/2024ppRef.toml)\n";
            return;
        }
    }
    if( mode == RunMode::HardProbes ){
        trees[kTrigIdx] = ( TTree* )fi->Get( cfg.trigTreePath );
        if( !trees[kTrigIdx] ){
            std::cerr << "Missing HLT tree in " << input << "\n(check cfg/2024ppRef.toml)\n";
            return;
        }
    }

    // branch mapping
    const bool isMC = ( mode == RunMode::MC );
    SetBranches( trees[kEvtIdx], event.BranchMap( isMC ) );
    for( size_t c = 0; c < nCones; c++ ){ SetBranches( trees[c], jets[c].BranchMap( isMC ) ); }
    if( mode != RunMode::MC ){ SetBranches( trees[kSkimIdx], filters.BranchMap( cfg.filterBranch ) ); }
    if( mode == RunMode::HardProbes ){ SetBranches( trees[kTrigIdx], {{cfg.hltJ80Branch, &hlt_j80}} ); }

    // event histograms
    TH1D* hvz_all = new TH1D( "hvz_all", "all events;v_{z} (cm);N", 40, -20, 20 );
    TH1D* hvz = new TH1D( "hvz", "after vz+filter;v_{z} (cm);N", 40, -20, 20 );
    TH1I* hfilt = ( mode != RunMode::MC ) ? new TH1I( "hfilt", "ppvF;filter;N", 2, 0, 2 ) : nullptr;
    TH1I* h_j80 = ( mode == RunMode::HardProbes ) ? new TH1I( "h_hlt_j80", "HLT_AK4PFJet80;bit;N", 2, 0, 2 ) : nullptr;

    // histograms
    BinningConfig bins;
    std::vector<ConeHistograms> cones( nCones );
    for( size_t c = 0; c < nCones; c++ ){ cones[c].Init( cfg.coneLabels[c], bins ); }

    // corrected pT buffer: corrPt[cone][jet]
    std::vector<std::vector<float>> corrPt( nCones, std::vector<float>( kNRefMax, 0.0f ) );

    // sorted leading-jet indices
    std::vector<SortedJets> sorted( nCones );

    // event loop
    const Long64_t nEvents = trees[kEvtIdx]->GetEntries();
    const Long64_t nLoop = ( maxEvents < 0 ) ? nEvents : std::min( maxEvents, nEvents );

    for( Long64_t i = 0; i < nLoop; i++ ){

        // vz
        trees[kEvtIdx]->GetEntry( i );
        hvz_all->Fill( event.vz );
        if( TMath::Abs( event.vz ) > kVzCut ){ continue; }

        // ppvF filter
        if( mode != RunMode::MC ){
            trees[kSkimIdx]->GetEntry( i );
            hfilt->Fill( filters.ppvF );
            if( filters.ppvF == 0 ){ continue; }
        }
        hvz->Fill( event.vz );

        // trig
        if( mode == RunMode::HardProbes ){
            trees[kTrigIdx]->GetEntry( i );
            h_j80->Fill( hlt_j80 );
        }

        // applying golden JSON
        if( dcs && !dcs->isGood( event.run, event.lumi ) ){ continue; }

        // jet trees (nref)
        for( size_t c = 0; c < nCones; c++ ){ trees[c]->GetEntry( i ); }
        if( jets[trigConeIdx].reco.nref < 2 ){ continue; }

        const float weight = event.w;

        // JEC: fill corrPt[c][j] for every cone and jet
        for( size_t c = 0; c < nCones; c++ ){
            for( int j = 0; j < jets[c].reco.nref; j++ ){
                jecs[c].SetJetPT( jets[c].reco.rawpt[j] );
                jecs[c].SetJetEta( jets[c].reco.eta[j] );
                jecs[c].SetJetPhi( jets[c].reco.phi[j] );
                corrPt[c][j] = ( float )jecs[c].GetCorrectedPT();
            }
        }

        // sort
        for( size_t c = 0; c < nCones; c++ ){ sorted[c] = FindLeadingJets( corrPt[c].data(), jets[c].reco.nref ); }

        // trigger efficiency cut
        if( mode == RunMode::HardProbes ){
            if( hlt_j80 == 0 ){ continue; }
            if( hlt_j80 == 1 && corrPt[trigConeIdx][sorted[trigConeIdx].lead] <= cfg.hltJ80Thresh ){ continue; }
        }

        // per-cone: jet ID → dijet → |A| → fill
        for( size_t c = 0; c < nCones; c++ ){

            // incl jets: all corrected jets in this cone passing cfg.minJetPt
            for( int j = 0; j < jets[c].reco.nref; j++ ){
                if( corrPt[c][j] >= cfg.minJetPt ){ cones[c].FillInclJet( corrPt[c][j], jets[c].reco.eta[j], jets[c].reco.phi[j], weight ); }
            }

            if( sorted[c].sublead == -1 ){ continue; }

            // jet ID
            auto& pf = jets[c].reco.pf;
            if( !js.JetSelection( jets[c].reco.eta[sorted[c].lead], jets[c].reco.phi[sorted[c].lead],
                                pf.CHF[sorted[c].lead], pf.NHF[sorted[c].lead], pf.CEF[sorted[c].lead],
                                pf.NEF[sorted[c].lead], pf.MUF[sorted[c].lead],
                                pf.CHM[sorted[c].lead], pf.NHM[sorted[c].lead], pf.CEM[sorted[c].lead],
                                pf.NEM[sorted[c].lead], pf.MUM[sorted[c].lead] ) ){ continue; }

            if( !js.JetSelection( jets[c].reco.eta[sorted[c].sublead], jets[c].reco.phi[sorted[c].sublead],
                                pf.CHF[sorted[c].sublead], pf.NHF[sorted[c].sublead], pf.CEF[sorted[c].sublead],
                                pf.NEF[sorted[c].sublead], pf.MUF[sorted[c].sublead],
                                pf.CHM[sorted[c].sublead], pf.NHM[sorted[c].sublead], pf.CEM[sorted[c].sublead],
                                pf.NEM[sorted[c].sublead], pf.MUM[sorted[c].sublead] ) ){ continue; }

            bool hasThird = ( sorted[c].third != -1 );
            if( hasThird ){
                hasThird = js.JetSelection( jets[c].reco.eta[sorted[c].third], jets[c].reco.phi[sorted[c].third],
                                pf.CHF[sorted[c].third], pf.NHF[sorted[c].third], pf.CEF[sorted[c].third],
                                pf.NEF[sorted[c].third], pf.MUF[sorted[c].third],
                                pf.CHM[sorted[c].third], pf.NHM[sorted[c].third], pf.CEM[sorted[c].third],
                                pf.NEM[sorted[c].third], pf.MUM[sorted[c].third] );
            }

            // dijet logic (see include/Dijet.h)
            DijetResult dijet = MakeDijet( sorted[c], hasThird, corrPt[c].data(), jets[c].reco.eta, jets[c].reco.phi,
                                          event.event, cfg.minJetPt, cfg.dphiCut );
            if( !dijet.valid ){ continue; }

            // dijet-level |A| acceptance cut
            if( TMath::Abs( dijet.A ) > cfg.maxAbsA ){ continue; }

            // fill histograms
            cones[c].Fill( dijet, corrPt[c].data(), jets[c].reco.eta, jets[c].reco.phi, weight );
        }
    }

    // output
    { Ssiz_t sl = output.Last( '/' ); if( sl != kNPOS ){ gSystem->mkdir( TString( output( 0, sl ) ), kTRUE ); } }
    TFile* fo = new TFile( output, "recreate" );
    fo->cd();
    hvz_all->Write();
    hvz->Write();
    if( hfilt ){ hfilt->Write(); }
    if( h_j80 ){ h_j80->Write(); }
    for( size_t c = 0; c < nCones; c++ ){
        TDirectory* dir = fo->mkdir( cfg.coneLabels[c].Data() );
        cones[c].Write( dir );
        fo->cd();
    }
    fo->Close();
    fi->Close();
}
