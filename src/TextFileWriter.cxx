#include "TextFileWriter.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

#include "Binning.h"
#include "Naming.h"
#include "AnalysisConfig.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

// Minimum valid pT slices needed to attempt a 3-parameter fit.
static constexpr int kMinSlices = 3;

static constexpr double kPtLo = 40.0;
static constexpr double kPtHi = 1000.0;
static constexpr int    kNPar = 3;

// Fit function: 1/(p0 + p1*log10(0.01*x) + p2/(x/10))
// Matches the formula declared in the JEC header line below.
static Double_t FitFunc( Double_t* x, Double_t* par ){
    double denom = par[0] + par[1] * TMath::Log10( 0.01 * x[0] ) + par[2] / ( x[0] / 10.0 );
    return ( denom != 0.0 ) ? 1.0 / denom : 1.0;
}

// Level label and JEC formula declaration — must match FitFunc exactly.
static const char* const kJECHeader =
    "{1 JetEta 1 JetPt 1./([p0]+[p1]*log10(0.01*x)+[p2]/(x/10.0)) Correction L2Residual}";

struct FitResult {
    double p[kNPar] = {1.0, 0.0, 0.0};
    bool   valid    = false;
};

// Arithmetic midpoint of a pT slice — used as the fit x-value for that slice.
static double SliceCenter( const RangeBin& sl ){
    return 0.5 * ( sl.lo + sl.hi );
}

// Fetch an intercept histogram from the current cone TDirectory layout.
static TH1D* FetchIntercept( TFile* f, const TString& cone, const TString& name ){
    TDirectory* coneDir = ( TDirectory* )f->Get( cone );
    return coneDir ? ( TH1D* )coneDir->Get( name ) : nullptr;
}

// Fit corr(pT_avg) at one eta bin, write the TGraphErrors — with the TF1
// embedded so it draws later — into dGraphs, and return the fit parameters.
// The graph and the text-file output are always derived from this single
// fit; never re-fit the same points twice elsewhere.
static FitResult FitPtSlices(
    const std::vector<double>& ptCenters,
    const std::vector<double>& corr,
    const std::vector<double>& corrErr,
    const TString& graphName,
    TDirectory* dGraphs ){
    FitResult r;
    int n = ( int )ptCenters.size();
    if( n < kMinSlices ) return r;

    std::vector<double> ex( n, 0.0 );
    TGraphErrors* gr = new TGraphErrors( n,
        ptCenters.data(), corr.data(), ex.data(), corrErr.data() );
    gr->SetName( graphName );
    gr->SetTitle( ";p_{T,avg} [GeV];Correction factor" );

    TF1* f = new TF1( graphName + "_fit", FitFunc, kPtLo, kPtHi, kNPar );
    // Start at the unity correction (1/(1+0+0)=1) — physically close for all eta.
    // Starting at (1.5,1.5,1.5) gives a negative denominator at low pT.
    f->SetParameter( 0, 1.0 );
    f->SetParameter( 1, 0.0 );
    f->SetParameter( 2, 0.0 );

    // No "N": the fit function is embedded in the graph so runPlotting can draw it later.
    TFitResultPtr res = gr->Fit( f, "QSR" );
    if( res.Get() && res->IsValid() ){
        for( int p = 0; p < kNPar; p++ ) r.p[p] = res->Parameter( p );
        r.valid = true;
    }

    if( dGraphs ){ dGraphs->cd(); gr->Write(); }
    delete f;    // clone embedded in gr's function list is separately owned
    delete gr;
    return r;
}

// Write one CMS L2Residual JEC data line for the [etaLo, etaHi] range.
static void WriteJECLine( std::ofstream& out, double etaLo, double etaHi, const FitResult& fit ){
    out << etaLo << "\t" << etaHi << "\t" << ( kNPar + 2 )
        << "\t" << kPtLo << "\t" << kPtHi;
    if( fit.valid ){
        out << "\t" << fit.p[0] << "\t" << fit.p[1] << "\t" << fit.p[2];
    } else {
        out << "\t" << 1 << "\t" << 0 << "\t" << 0;
    }
    out << "\n";
}

// Mirrored |eta| text file: negative half outermost→innermost, then positive
// half innermost→outermost, both halves reusing the same |eta| fit results.
static bool WriteAbsEtaTextFile( const TString& path, const std::vector<FitResult>& fits ){
    std::ofstream out( path.Data() );
    if( !out.is_open() ) return false;
    out << kJECHeader << "\n";
    const int nEta = ( int )fits.size();
    for( int ieta = nEta - 1; ieta >= 0; ieta-- )
        WriteJECLine( out, -kAbsEtaEdges[ieta + 1], -kAbsEtaEdges[ieta], fits[ieta] );
    for( int ieta = 0; ieta < nEta; ieta++ )
        WriteJECLine( out, kAbsEtaEdges[ieta], kAbsEtaEdges[ieta + 1], fits[ieta] );
    out.close();
    return true;
}

// Independent full-eta text file: kEtaEdges is already ascending -5.191→5.191,
// so this is a single direct pass — no mirroring, each bin has its own fit.
static bool WriteFullEtaTextFile( const TString& path, const std::vector<FitResult>& fits ){
    std::ofstream out( path.Data() );
    if( !out.is_open() ) return false;
    out << kJECHeader << "\n";
    const int nEta = ( int )fits.size();
    for( int ieta = 0; ieta < nEta; ieta++ )
        WriteJECLine( out, kEtaEdges[ieta], kEtaEdges[ieta + 1], fits[ieta] );
    out.close();
    return true;
}

void runTextFile( TString hpResidualsFile, TString zbResidualsFile,
                 TString outputRootFile, TString outputTextPrefix,
                 TString method, bool useNorm ){

    const TString suffix = useNorm ? "_norm" : "";

    std::cout << "HP residuals: " << hpResidualsFile  << "\n"
              << "ZB residuals: " << zbResidualsFile  << "\n"
              << "Output ROOT:  " << outputRootFile   << "\n"
              << "Text prefix:  " << outputTextPrefix << "\n"
              << "Method:       " << method           << "\n"
              << "Normalized:   " << ( useNorm ? "yes (kFSR-norm)" : "no (direct)" ) << "\n";

    const AnalysisConfig& cfg = Config();
    PrintConfigSummary( cfg );

    TFile* fHP = TFile::Open( hpResidualsFile, "read" );
    TFile* fZB = TFile::Open( zbResidualsFile, "read" );
    if( !fHP || fHP->IsZombie() ){ std::cerr << "Cannot open " << hpResidualsFile << "\n"; return; }
    if( !fZB || fZB->IsZombie() ){ std::cerr << "Cannot open " << zbResidualsFile << "\n"; return; }

    { Ssiz_t sl = outputRootFile.Last( '/' ); if( sl != kNPOS ){ gSystem->mkdir( TString( outputRootFile( 0, sl ) ), kTRUE ); } }
    { Ssiz_t sl = outputTextPrefix.Last( '/' ); if( sl != kNPOS ){ gSystem->mkdir( TString( outputTextPrefix( 0, sl ) ), kTRUE ); } }

    TFile* fOut = new TFile( outputRootFile, "recreate" );

    BinningConfig bins;
    const int nPt = ( int )bins.ptavgSlices.size();

    // pT_avg bin edges for the corrfinal grid — bins.ptavgSlices are contiguous.
    std::vector<Double_t> ptEdges( nPt + 1 );
    ptEdges[0] = bins.ptavgSlices[0].lo;
    for( int ip = 0; ip < nPt; ip++ ) ptEdges[ip + 1] = bins.ptavgSlices[ip].hi;

    for( const TString& cone : cfg.coneLabels ){

        TDirectory* hpConeDir = ( TDirectory* )fHP->Get( cone );
        TDirectory* zbConeDir = ( TDirectory* )fZB->Get( cone );
        if( !hpConeDir || !zbConeDir ){
            std::cerr << "Missing " << cone << " directory in HP or ZB residuals file, skipping\n";
            continue;
        }

        fOut->cd();
        TDirectory* coneDirOut = fOut->mkdir( cone.Data() );
        TDirectory* dGraphs = coneDirOut->mkdir( "graphs" );

        // filled below for abseta then fulleta, used to write the two text files afterward
        std::vector<FitResult> fitsAbsEta, fitsFullEta;

        for( int em = 0; em < 2; em++ ){
            const bool fullEta = ( em == 1 );
            const std::vector<Double_t>& etaEdges = fullEta ? kEtaEdges : kAbsEtaEdges;
            const int nEta = ( int )etaEdges.size() - 1;
            const TString etaMode = L2Name::EtaModeKey( fullEta );

            // one intercept histogram per pT slice, chosen from HP or ZB by the trigger threshold
            std::vector<TH1D*> hSlice( nPt, nullptr );
            int nMissingSlices = 0;
            for( int ip = 0; ip < nPt; ip++ ){
                const auto& ptSlice = bins.ptavgSlices[ip];
                TFile* src = ( ptSlice.lo >= cfg.hltJ80Thresh ) ? fHP : fZB;
                TString name = L2Name::ObjectName( cone, "intercept",
                    {etaMode, L2Name::PtKey( ptSlice )}, {method} ) + suffix;
                hSlice[ip] = FetchIntercept( src, cone, name );
                if( !hSlice[ip] ){ nMissingSlices++; }
            }
            if( nMissingSlices > 0 ){
                std::cerr << cone << " " << etaMode << ": " << nMissingSlices << "/" << nPt
                          << " pT slices missing (intercept histogram not found)\n";
            }

            // final corrections grid: x = eta/|eta|, y = pT_avg slices, z = correction
            TString gridName = L2Name::ObjectName( cone, "corrfinal", {etaMode}, {method} ) + suffix;
            TH2D* hGrid = new TH2D( gridName, "", nEta, etaEdges.data(), nPt, ptEdges.data() );
            hGrid->GetXaxis()->SetTitle( fullEta ? "#eta_{probe}" : "|#eta_{probe}|" );
            hGrid->GetYaxis()->SetTitle( "p_{T,avg} [GeV]" );
            hGrid->GetZaxis()->SetTitle( "Correction factor" );
            hGrid->Sumw2();
            for( int ip = 0; ip < nPt; ip++ ){
                if( !hSlice[ip] ) continue;
                for( int ieta = 0; ieta < nEta; ieta++ ){
                    hGrid->SetBinContent( ieta + 1, ip + 1, hSlice[ip]->GetBinContent( ieta + 1 ) );
                    hGrid->SetBinError( ieta + 1, ip + 1, hSlice[ip]->GetBinError( ieta + 1 ) );
                }
            }
            coneDirOut->cd();
            hGrid->Write();
            delete hGrid;

            // pT-dependence fit, one per eta bin
            std::vector<FitResult>& fits = fullEta ? fitsFullEta : fitsAbsEta;
            fits.assign( nEta, FitResult{} );
            int nUnityFallback = 0;
            for( int ieta = 0; ieta < nEta; ieta++ ){
                std::vector<double> ptX, corr, corrErr;
                for( int ip = 0; ip < nPt; ip++ ){
                    if( !hSlice[ip] ) continue;
                    double v = hSlice[ip]->GetBinContent( ieta + 1 );
                    double e = hSlice[ip]->GetBinError( ieta + 1 );
                    if( v == 0.0 && e == 0.0 ) continue;
                    ptX.push_back( SliceCenter( bins.ptavgSlices[ip] ) );
                    corr.push_back( v );
                    corrErr.push_back( e > 0.0 ? e : 1e-4 );
                }
                TString graphName = L2Name::ObjectName( cone, "ptcorr",
                    {etaMode, L2Name::EtaKey( ieta, fullEta )}, {method} ) + suffix;
                fits[ieta] = FitPtSlices( ptX, corr, corrErr, graphName, dGraphs );
                if( !fits[ieta].valid ){ nUnityFallback++; }
            }
            if( nUnityFallback > 0 ){
                std::cerr << cone << " " << etaMode << ": " << nUnityFallback << "/" << nEta
                          << " eta bins fell back to unity (fewer than " << kMinSlices << " pT slices had data)\n";
            }
        }

        TString absEtaTxt = outputTextPrefix + "_" + cone + "_abseta" + suffix + ".txt";
        TString etaTxt    = outputTextPrefix + "_" + cone + "_eta" + suffix + ".txt";
        if( !WriteAbsEtaTextFile( absEtaTxt, fitsAbsEta ) )
            std::cerr << "Cannot open output file " << absEtaTxt << "\n";
        if( !WriteFullEtaTextFile( etaTxt, fitsFullEta ) )
            std::cerr << "Cannot open output file " << etaTxt << "\n";

        std::cout << "Done. " << cone << ": " << 2 * ( int )fitsAbsEta.size() << " eta bins -> " << absEtaTxt
                  << ", " << ( int )fitsFullEta.size() << " eta bins -> " << etaTxt << "\n";
    }

    fOut->Close();
    fHP->Close();
    fZB->Close();
}

void runTextFile( TString residualsFile, TString outputRootFile,
                 TString outputTextPrefix,
                 TString method, bool useNorm ){
    std::cout << "Single dataset mode (no HP/ZB merge) — using " << residualsFile
              << " for every pT_avg slice\n";
    runTextFile( residualsFile, residualsFile, outputRootFile, outputTextPrefix, method, useNorm );
}
