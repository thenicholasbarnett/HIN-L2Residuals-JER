#include "TextFileWriter.h"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
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

static FitResult FitPtSlices(
    const std::vector<double>& ptCenters,
    const std::vector<double>& corr,
    const std::vector<double>& corrErr,
    int iEta ){
    FitResult r;
    int n = ( int )ptCenters.size();
    if( n < kMinSlices ) return r;

    std::vector<double> ex( n, 0.0 );
    TGraphErrors* gr = new TGraphErrors( n,
        ptCenters.data(), corr.data(), ex.data(), corrErr.data() );

    TF1* f = new TF1( Form( "_ptfit_%d", iEta ), FitFunc, kPtLo, kPtHi, kNPar );
    // Start at the unity correction (1/(1+0+0)=1) — physically close for all eta.
    // Starting at (1.5,1.5,1.5) gives a negative denominator at low pT.
    f->SetParameter( 0, 1.0 );
    f->SetParameter( 1, 0.0 );
    f->SetParameter( 2, 0.0 );

    TFitResultPtr res = gr->Fit( f, "NQSR" );
    if( res.Get() && res->IsValid() ){
        for( int p = 0; p < kNPar; p++ ) r.p[p] = res->Parameter( p );
        r.valid = true;
    }

    delete f;
    delete gr;
    return r;
}

// Arithmetic midpoint of a pT slice — used as the fit x-value for that slice.
static double SliceCenter( const RangeBin& sl ){
    return 0.5 * ( sl.lo + sl.hi );
}

void runTextFile( TString residualsFile, TString outputFile,
                 TString method, TString cone ){

    std::cout << "Residuals: " << residualsFile << "\n"
              << "Output:    " << outputFile    << "\n"
              << "Cone:      " << cone          << "\n"
              << "Method:    " << method        << "\n";

    TFile* fIn = TFile::Open( residualsFile, "read" );
    if( !fIn || fIn->IsZombie() ){
        std::cerr << "Cannot open " << residualsFile << "\n";
        return;
    }

    BinningConfig bins;
    const int nPt  = ( int )bins.ptavgSlices.size();
    const int nEta = ( int )kAbsEtaEdges.size() - 1;  // 18

    std::vector<double> ptCenters( nPt );
    for( int ip = 0; ip < nPt; ip++ ){
        ptCenters[ip] = SliceCenter( bins.ptavgSlices[ip] );
    }

    // Load one intercept TH1D per pT slice. Prefer the canonical global order:
    // {cone}_intercept_abseta_{ptSlice}_{method}
    // and keep the older {cone}_intercept_{method}{ptSlice} form readable.
    TDirectory* coneDir = ( TDirectory* )fIn->Get( cone.Data() );
    std::vector<TH1D*> hSlice( nPt, nullptr );
    for( int ip = 0; ip < nPt; ip++ ){
        TString name = L2Name::ObjectName( cone, "intercept",
            {L2Name::EtaModeKey( false ), L2Name::PtKey( bins.ptavgSlices[ip] )}, {method} );
        if( coneDir ) hSlice[ip] = ( TH1D* )coneDir->Get( name );
        if( !hSlice[ip] ) hSlice[ip] = ( TH1D* )fIn->Get( name );
        if( !hSlice[ip] ){
            TString oldName = Form( "%s_intercept_%s%s",
                cone.Data(), method.Data(), bins.ptavgSlices[ip].shortName.Data() );
            hSlice[ip] = ( TH1D* )fIn->Get( oldName );
            if( !hSlice[ip] ){
                std::cerr << "WARNING: " << name << " not found\n";
            }
        }
    }

    // For each |eta| bin, gather valid (pT, correction, error) points and fit.
    std::vector<FitResult> fits( nEta );
    for( int ieta = 0; ieta < nEta; ieta++ ){
        std::vector<double> ptX, corr, corrErr;
        for( int ip = 0; ip < nPt; ip++ ){
            if( !hSlice[ip] ) continue;
            double v = hSlice[ip]->GetBinContent( ieta + 1 );
            double e = hSlice[ip]->GetBinError( ieta + 1 );
            if( v == 0.0 && e == 0.0 ) continue;
            ptX.push_back( ptCenters[ip] );
            corr.push_back( v );
            corrErr.push_back( e > 0.0 ? e : 1e-4 );
        }
        fits[ieta] = FitPtSlices( ptX, corr, corrErr, ieta );
        if( !fits[ieta].valid ){
            std::cerr << "WARNING: fit failed for |eta| bin " << ieta + 1
                      << " [" << kAbsEtaEdges[ieta] << ", " << kAbsEtaEdges[ieta + 1] << "]"
                      << " (" << ( int )ptX.size() << " pT slices available)\n";
        }
    }

    fIn->Close();

    std::ofstream out( outputFile.Data() );
    if( !out.is_open() ){
        std::cerr << "Cannot open output file " << outputFile << "\n";
        return;
    }

    out << kJECHeader << "\n";

    // Write negative eta half: |eta| bins from outermost (index 17) to innermost (index 0).
    // For |eta| bin i: eta range is [-kAbsEtaEdges[i+1], -kAbsEtaEdges[i]].
    for( int ieta = nEta - 1; ieta >= 0; ieta-- ){
        const FitResult& fit = fits[ieta];
        double etaLo = -kAbsEtaEdges[ieta + 1];
        double etaHi = -kAbsEtaEdges[ieta];
        out << etaLo << "\t" << etaHi << "\t" << ( kNPar + 2 )
            << "\t" << kPtLo << "\t" << kPtHi;
        if( fit.valid ){
            out << "\t" << fit.p[0] << "\t" << fit.p[1] << "\t" << fit.p[2];
        } else {
            out << "\t" << 1 << "\t" << 0 << "\t" << 0;
        }
        out << "\n";
    }

    // Write positive eta half: |eta| bins from innermost (index 0) to outermost (index 17).
    for( int ieta = 0; ieta < nEta; ieta++ ){
        const FitResult& fit = fits[ieta];
        double etaLo = kAbsEtaEdges[ieta];
        double etaHi = kAbsEtaEdges[ieta + 1];
        out << etaLo << "\t" << etaHi << "\t" << ( kNPar + 2 )
            << "\t" << kPtLo << "\t" << kPtHi;
        if( fit.valid ){
            out << "\t" << fit.p[0] << "\t" << fit.p[1] << "\t" << fit.p[2];
        } else {
            out << "\t" << 1 << "\t" << 0 << "\t" << 0;
        }
        out << "\n";
    }

    out.close();
    std::cout << "Done. " << 2 * nEta << " eta bins written to " << outputFile << "\n";
}
