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
#include "jetmet_jer/JetResolutionObject.h"

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>

// Minimum valid pT slices needed to attempt a 3-parameter fit.
static constexpr int kMinSlices = 3;

// outputTag default and fixed output location — see TextFileWriter.h.
static const char* const kDefaultTag = "L2Residual";
static const char* const kTextOutputSubdir = "data/jec/preliminary";

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

// --- JER SF text output, via the vendored JME::JetResolutionObject
// (external/jetmet_jer/, from CMSSW CondFormats/JetMETObjects) ---
//
// Unlike the JEC writer above (which fits a continuous correction(pT_avg)
// function per eta bin), the JER SF file is a direct binned grid: one flat
// value per (eta bin, pT_avg slice) cell, no fit. This matches how real
// official CMS JER SF text files are actually shaped -- a JetEta x JetPt
// grid of flat scale factors, not a parametrized formula. Definition line:
// "{1 JetEta 1 JetPt [0] Resolution}" -- 1 bin variable (JetEta), 1
// structural variable (JetPt, using each pT_avg slice's own [lo,hi) edges
// as its range -- required by the record format even though the formula
// itself doesn't depend on it, since Step 2 only gives discretized pT_avg
// slices, not a continuous fit vs pT), formula "[0]" (the record's first
// parameter *is* the flat JER SF value, retrievable either via
// evaluateFormula() or directly). A second parameter, unc -- the *fractional*
// uncertainty on the SF, unc = error/value, not the raw absolute fit error --
// rides along in every record, unused by evaluateFormula() but directly
// retrievable via record.getParametersValues()[1]. Fractional, not absolute,
// so it matches the standard JER-smearing convention directly (s_up/down =
// sJER * (1 +/- unc), e.g. "Practical Application of JER Smearing") without
// a consumer needing to divide by the SF value themselves.
struct JerRecord {
    double ptLo, ptHi, value, unc;
};

// This eta bin's per-pT-slice JER SF values, skipping slices with no data
// (missing/excluded source, or a genuinely unfit bin left at the histogram's
// zero/zero-error default -- same convention as FitPtSlices's point list above).
static std::vector<JerRecord> CollectJerRecords( int ieta,
    const std::vector<RangeBin>& ptSlices, const std::vector<TH1D*>& hSliceJer ){
    std::vector<JerRecord> out;
    for( size_t ip = 0; ip < hSliceJer.size(); ip++ ){
        if( !hSliceJer[ip] ) continue;
        double v = hSliceJer[ip]->GetBinContent( ieta + 1 );
        double e = hSliceJer[ip]->GetBinError( ieta + 1 );
        if( v == 0.0 && e == 0.0 ) continue;
        double unc = ( v != 0.0 ) ? e / v : 0.0;
        out.push_back( { ptSlices[ip].lo, ptSlices[ip].hi, v, unc } );
    }
    return out;
}

// One record line per JerRecord. The "4" token is 2*nVariables + nParameters
// (2 for the JetPt range this record structurally carries, 2 for the actual
// [value, unc] parameters) -- see JetResolutionObject::Record's own
// parsing for why that combined count, not just nParameters, is what the
// format's record header field means.
static void AppendJerLines( std::stringstream& buf, double etaLo, double etaHi,
    const std::vector<JerRecord>& recs ){
    for( const auto& r : recs ){
        buf << etaLo << " " << etaHi << " 4 " << r.ptLo << " " << r.ptHi
            << " " << r.value << " " << r.unc << "\n";
    }
}

// Assembles the JER SF text-format definition+records into a temp file, then
// round-trips it through the real vendored JME::JetResolutionObject: parsed
// via its file constructor, re-emitted via its own saveToFile(). The final
// on-disk bytes are produced by the vendored code, not reimplemented here --
// this function only assembles the input text and mirrors abseta cells onto
// both eta halves, exactly like WriteAbsEtaTextFile does for JEC (both files
// use "JetEta" binning even in "abseta" mode, for the same reason: the data
// is |eta|-binned, but every record still needs one real eta range).
static bool WriteJerSfTextFile( const TString& path, bool fullEta,
    const std::vector<RangeBin>& ptSlices, const std::vector<TH1D*>& hSliceJer ){

    std::stringstream buf;
    buf << "{1 JetEta 1 JetPt [0] Resolution}\n";

    if( fullEta ){
        const int nEta = ( int )kEtaEdges.size() - 1;
        for( int ieta = 0; ieta < nEta; ieta++ )
            AppendJerLines( buf, kEtaEdges[ieta], kEtaEdges[ieta + 1],
                CollectJerRecords( ieta, ptSlices, hSliceJer ) );
    } else {
        const int nEta = ( int )kAbsEtaEdges.size() - 1;
        for( int ieta = nEta - 1; ieta >= 0; ieta-- )
            AppendJerLines( buf, -kAbsEtaEdges[ieta + 1], -kAbsEtaEdges[ieta],
                CollectJerRecords( ieta, ptSlices, hSliceJer ) );
        for( int ieta = 0; ieta < nEta; ieta++ )
            AppendJerLines( buf, kAbsEtaEdges[ieta], kAbsEtaEdges[ieta + 1],
                CollectJerRecords( ieta, ptSlices, hSliceJer ) );
    }

    TString tmpPath = path + ".tmp";
    {
        std::ofstream tmp( tmpPath.Data() );
        if( !tmp.is_open() ) return false;
        tmp << buf.str();
    }

    bool ok = true;
    try {
        JME::JetResolutionObject obj( tmpPath.Data() );
        obj.saveToFile( path.Data() );
    } catch( const std::exception& e ){
        std::cerr << "ERROR building JER SF text file " << path << ": " << e.what() << "\n";
        ok = false;
    }
    gSystem->Unlink( tmpPath );
    return ok;
}

// Which residuals file backs a given pT_avg slice. Merge is the original
// triggered+non-triggered behavior; TriggeredOnly/NonTriggeredOnly back the
// single-dataset overload -- they are NOT the same threshold branch with one
// side nulled out, because NonTriggeredOnly must ignore the threshold
// entirely (an unbiased dataset has nothing to gate), while TriggeredOnly
// must drop (not fall back on) slices below threshold.
enum class SourceMode { Merge, TriggeredOnly, NonTriggeredOnly };

static TFile* SelectSource( SourceMode mode, const RangeBin& ptSlice, double threshold,
                            TFile* fTrig, TFile* fNoTrig ){
    switch( mode ){
        case SourceMode::Merge:            return ( ptSlice.lo >= threshold ) ? fTrig : fNoTrig;
        case SourceMode::TriggeredOnly:    return ( ptSlice.lo >= threshold ) ? fTrig : nullptr;
        case SourceMode::NonTriggeredOnly: return fNoTrig;
    }
    return nullptr;
}

// Shared implementation for all three public entry points below. fTrig/fNoTrig
// are already-open files; exactly one is null in single-dataset mode
// (TriggeredOnly leaves fNoTrig null, NonTriggeredOnly leaves fTrig null) and
// both are non-null in Merge mode.
static void RunTextFileImpl(
    SourceMode mode, TFile* fTrig, TFile* fNoTrig,
    TString outputRootFile, TString outputTag,
    TString method, bool useNorm ){

    if( outputTag.IsNull() ) outputTag = kDefaultTag;
    if( outputTag.Contains( "/" ) ){
        std::cerr << "ERROR: TAG must be a plain name, not a path (it contains '/'): "
                  << outputTag << "\n";
        return;
    }

    const TString suffix = useNorm ? "_norm" : "";

    const AnalysisConfig& cfg = Config();
    PrintConfigSummary( cfg );

    // Correction text files always land here, not wherever the caller points
    // OUTPUT= — see TextFileWriter.h. Gitignored; created if missing.
    const TString textDir = TString( cfg.repoRoot.c_str() ) + "/" + kTextOutputSubdir;
    if( gSystem->mkdir( textDir, kTRUE ) < 0 && gSystem->AccessPathName( textDir ) ){
        std::cerr << "ERROR: cannot create text output directory: " << textDir << "\n";
        return;
    }

    static const char* const kModeLabel[] = {
        "merge (triggered + non-triggered)", "single (triggered-only)", "single (non-triggered-only)" };
    std::cout << "Mode: " << kModeLabel[( int )mode] << "\n"
              << "Triggered residuals:     " << ( fTrig ? fTrig->GetName() : "(none)" ) << "\n"
              << "Non-triggered residuals: " << ( fNoTrig ? fNoTrig->GetName() : "(none)" ) << "\n"
              << "Output ROOT:  " << outputRootFile   << "\n"
              << "Text output dir: " << textDir        << "\n"
              << "Tag:          " << outputTag << "\n"
              << "Method:       " << method           << "\n"
              << "Normalized:   " << ( useNorm ? "yes (kFSR-norm)" : "no (direct)" ) << "\n";

    { Ssiz_t sl = outputRootFile.Last( '/' ); if( sl != kNPOS ){ gSystem->mkdir( TString( outputRootFile( 0, sl ) ), kTRUE ); } }

    TFile* fOut = new TFile( outputRootFile, "recreate" );

    // [step3] eta_mode selects which of the two text files actually get
    // written; "both" (default) writes both, matching prior behavior.
    const bool wantAbsEta  = ( cfg.etaModeOutput != "eta" );
    const bool wantFullEta = ( cfg.etaModeOutput != "abseta" );

    BinningConfig bins;
    const int nPt = ( int )bins.ptavgSlices.size();

    // pT_avg bin edges for the corrfinal grid — bins.ptavgSlices are contiguous.
    std::vector<Double_t> ptEdges( nPt + 1 );
    ptEdges[0] = bins.ptavgSlices[0].lo;
    for( int ip = 0; ip < nPt; ip++ ) ptEdges[ip + 1] = bins.ptavgSlices[ip].hi;

    for( const TString& cone : cfg.coneLabels ){

        TDirectory* trigConeDir = fTrig ? ( TDirectory* )fTrig->Get( cone ) : nullptr;
        TDirectory* notrigConeDir = fNoTrig ? ( TDirectory* )fNoTrig->Get( cone ) : nullptr;
        if( ( fTrig && !trigConeDir ) || ( fNoTrig && !notrigConeDir ) ){
            std::cerr << "Missing " << cone << " directory in a residuals file, skipping\n";
            continue;
        }

        fOut->cd();
        TDirectory* coneDirOut = fOut->mkdir( cone.Data() );
        TDirectory* dGraphs = coneDirOut->mkdir( "graphs" );

        // filled below for abseta then fulleta, used to write the two text files afterward
        std::vector<FitResult> fitsAbsEta, fitsFullEta;
        std::vector<TH1D*> hSliceJerAbsEta, hSliceJerFullEta;

        for( int em = 0; em < 2; em++ ){
            const bool fullEta = ( em == 1 );
            if( fullEta && !wantFullEta ) continue;
            if( !fullEta && !wantAbsEta ) continue;
            const std::vector<Double_t>& etaEdges = fullEta ? kEtaEdges : kAbsEtaEdges;
            const int nEta = ( int )etaEdges.size() - 1;
            const TString etaMode = L2Name::EtaModeKey( fullEta );

            // one intercept histogram per pT slice, chosen per SelectSource() --
            // a null source (TriggeredOnly below threshold) means the slice is
            // intentionally excluded, not merely missing.
            std::vector<TH1D*> hSlice( nPt, nullptr );
            std::vector<TH1D*>& hSliceJer = fullEta ? hSliceJerFullEta : hSliceJerAbsEta;
            hSliceJer.assign( nPt, nullptr );
            int nMissingSlices = 0;
            for( int ip = 0; ip < nPt; ip++ ){
                const auto& ptSlice = bins.ptavgSlices[ip];
                TFile* src = SelectSource( mode, ptSlice, cfg.hltJ80Thresh, fTrig, fNoTrig );
                if( src ){
                    TString name = L2Name::ObjectName( cone, "intercept",
                        {etaMode, L2Name::PtKey( ptSlice )}, {method} ) + suffix;
                    hSlice[ip] = FetchIntercept( src, cone, name );
                    TString nameJer = L2Name::ObjectName( cone, "intercept_jer",
                        {etaMode, L2Name::PtKey( ptSlice )}, {method} ) + suffix;
                    hSliceJer[ip] = FetchIntercept( src, cone, nameJer );
                }
                if( !hSlice[ip] ){ nMissingSlices++; }
            }
            if( nMissingSlices > 0 ){
                std::cerr << cone << " " << etaMode << ": " << nMissingSlices << "/" << nPt
                          << " pT slices missing or excluded (below trigger threshold in"
                          << " triggered-only mode, or intercept histogram not found)\n";
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

        TString absEtaTxt = textDir + "/" + outputTag + "_" + cone + "_abseta" + suffix + ".txt";
        TString etaTxt    = textDir + "/" + outputTag + "_" + cone + "_eta" + suffix + ".txt";
        if( wantAbsEta && !WriteAbsEtaTextFile( absEtaTxt, fitsAbsEta ) )
            std::cerr << "Cannot open output file " << absEtaTxt << "\n";
        if( wantFullEta && !WriteFullEtaTextFile( etaTxt, fitsFullEta ) )
            std::cerr << "Cannot open output file " << etaTxt << "\n";

        std::cout << "Done. " << cone << ": ";
        if( wantAbsEta )  std::cout << 2 * ( int )fitsAbsEta.size() << " eta bins -> " << absEtaTxt;
        if( wantAbsEta && wantFullEta ) std::cout << ", ";
        if( wantFullEta ) std::cout << ( int )fitsFullEta.size() << " eta bins -> " << etaTxt;
        std::cout << "\n";

        TString absEtaJerTxt = textDir + "/" + outputTag + "_" + cone + "_abseta_jer" + suffix + ".txt";
        TString etaJerTxt    = textDir + "/" + outputTag + "_" + cone + "_eta_jer" + suffix + ".txt";
        if( wantAbsEta && !WriteJerSfTextFile( absEtaJerTxt, false, bins.ptavgSlices, hSliceJerAbsEta ) )
            std::cerr << "Cannot write JER SF output file " << absEtaJerTxt << "\n";
        if( wantFullEta && !WriteJerSfTextFile( etaJerTxt, true, bins.ptavgSlices, hSliceJerFullEta ) )
            std::cerr << "Cannot write JER SF output file " << etaJerTxt << "\n";

        std::cout << "Done (JER SF). " << cone << ": ";
        if( wantAbsEta )  std::cout << "-> " << absEtaJerTxt;
        if( wantAbsEta && wantFullEta ) std::cout << ", ";
        if( wantFullEta ) std::cout << "-> " << etaJerTxt;
        std::cout << "\n";
    }

    fOut->Close();
    if( fTrig )   fTrig->Close();
    if( fNoTrig ) fNoTrig->Close();
}

void runTextFile( TString triggeredResidualsFile, TString nonTriggeredResidualsFile,
                 TString outputRootFile, TString outputTag,
                 TString method, bool useNorm ){
    TFile* fTrig = TFile::Open( triggeredResidualsFile, "read" );
    TFile* fNoTrig = TFile::Open( nonTriggeredResidualsFile, "read" );
    if( !fTrig || fTrig->IsZombie() ){ std::cerr << "Cannot open " << triggeredResidualsFile << "\n"; return; }
    if( !fNoTrig || fNoTrig->IsZombie() ){ std::cerr << "Cannot open " << nonTriggeredResidualsFile << "\n"; return; }
    RunTextFileImpl( SourceMode::Merge, fTrig, fNoTrig, outputRootFile, outputTag, method, useNorm );
}

void runTextFile( TString residualsFile, SingleDatasetKind kind, TString outputRootFile,
                 TString outputTag,
                 TString method, bool useNorm ){
    TFile* f = TFile::Open( residualsFile, "read" );
    if( !f || f->IsZombie() ){ std::cerr << "Cannot open " << residualsFile << "\n"; return; }
    if( kind == SingleDatasetKind::Triggered ){
        RunTextFileImpl( SourceMode::TriggeredOnly, f, nullptr, outputRootFile, outputTag, method, useNorm );
    } else {
        RunTextFileImpl( SourceMode::NonTriggeredOnly, nullptr, f, outputRootFile, outputTag, method, useNorm );
    }
}

void runTextFilePtResolution( TString responseFile, TString outputTag ){
    const AnalysisConfig& cfg = Config();

    TFile* fIn = TFile::Open( responseFile, "read" );
    if( !fIn || fIn->IsZombie() ){
        std::cerr << "Cannot open response file: " << responseFile << "\n";
        return;
    }

    if( outputTag.IsNull() ) outputTag = "JER_ptresolution";
    if( outputTag.Contains( "/" ) ){
        std::cerr << "ERROR: outputTag must not contain '/': " << outputTag << "\n";
        fIn->Close();
        return;
    }

    TString textDir = TString( cfg.repoRoot.c_str() ) + "/" + kTextOutputSubdir;
    if( gSystem->mkdir( textDir, kTRUE ) < 0 && gSystem->AccessPathName( textDir ) ){
        std::cerr << "ERROR: cannot create text output directory: " << textDir << "\n";
        fIn->Close();
        return;
    }

    const bool wantAbsEta  = ( cfg.etaModeOutput != "eta" );
    const bool wantFullEta = ( cfg.etaModeOutput != "abseta" );

    const int nEta = ( int )kAbsEtaEdges.size() - 1;

    // "corr" variant, "incl" collection -- these are the names written by
    // ExtractPerAbsEtaVsPtGen in ResponseExtractor.cxx.
    const TString variant = "corr";
    const TString collection = "incl";

    for( const TString& cone : cfg.coneLabels ){
        TDirectory* dPerEta = ( TDirectory* )fIn->Get( cone + "/JER_per_abseta" );
        if( !dPerEta ){
            std::cerr << "No JER_per_abseta/ directory for " << cone
                      << " -- was this file produced by runResponse?\n";
            continue;
        }

        // Collect JER records per |eta| bin, reading the per-|eta| JER TH1Ds.
        // Each TH1D has one bin per pT_gen interval; bin content = JER, error = JER uncertainty.
        std::vector<std::vector<JerRecord>> etaRecords( nEta );
        bool any = false;
        for( int ieta = 0; ieta < nEta; ieta++ ){
            TString etaKey = L2Name::EtaKey( ieta, false );
            TString hName  = L2Name::ObjectName( cone, "JER",
                { variant, "vs_ptgen", etaKey }, { collection } );
            TH1D* h = ( TH1D* )dPerEta->Get( hName );
            if( !h ) continue;
            for( int ip = 1; ip <= h->GetNbinsX(); ip++ ){
                const double v = h->GetBinContent( ip );
                const double e = h->GetBinError( ip );
                if( v == 0.0 && e == 0.0 ) continue;
                const double ptLo = h->GetXaxis()->GetBinLowEdge( ip );
                const double ptHi = h->GetXaxis()->GetBinUpEdge( ip );
                const double unc  = ( v != 0.0 ) ? e / v : 0.0;
                etaRecords[ieta].push_back( { ptLo, ptHi, v, unc } );
                any = true;
            }
        }

        if( !any ){
            std::cerr << "No per-|eta| JER data found for " << cone << "\n";
            continue;
        }

        // Write abseta file (mirror |eta| bins onto both eta halves).
        if( wantAbsEta ){
            TString path = textDir + "/" + outputTag + "_" + cone + "_abseta_ptresolution.txt";
            std::stringstream buf;
            buf << "{1 JetEta 1 JetPt [0] Resolution}\n";
            for( int ieta = nEta - 1; ieta >= 0; ieta-- )
                AppendJerLines( buf, -kAbsEtaEdges[ieta + 1], -kAbsEtaEdges[ieta], etaRecords[ieta] );
            for( int ieta = 0; ieta < nEta; ieta++ )
                AppendJerLines( buf, kAbsEtaEdges[ieta], kAbsEtaEdges[ieta + 1], etaRecords[ieta] );

            TString tmpPath = path + ".tmp";
            { std::ofstream tmp( tmpPath.Data() ); tmp << buf.str(); }
            bool ok = true;
            try {
                JME::JetResolutionObject obj( tmpPath.Data() );
                obj.saveToFile( path.Data() );
            } catch( const std::exception& ex ){
                std::cerr << "ERROR writing pT resolution file " << path << ": " << ex.what() << "\n";
                ok = false;
            }
            gSystem->Unlink( tmpPath );
            if( ok ) std::cout << "Done (pT resolution). " << cone << ": -> " << path << "\n";
        }

        // Write full-eta file (one record per eta bin using the positive-eta JER for each).
        // Without per-bin negative-eta extraction (the folding sums both sides), full-eta
        // output is the same as abseta but written with explicit full-eta bin edges.
        // A future pass that extracts eta_reco > 0 and < 0 separately can replace this.
        if( wantFullEta ){
            TString path = textDir + "/" + outputTag + "_" + cone + "_eta_ptresolution.txt";
            std::stringstream buf;
            buf << "{1 JetEta 1 JetPt [0] Resolution}\n";
            // Use the |eta|-folded values mirrored across eta=0 -- same content as abseta file,
            // just written with full-eta JEC bin edges. Matches the convention for eta text files
            // elsewhere in TextFileWriter.cxx.
            const int nFullEta = ( int )kEtaEdges.size() - 1;
            for( int ieta = 0; ieta < nFullEta; ieta++ ){
                const int iAbs = ( ieta < nFullEta / 2 )
                    ? ( nFullEta / 2 - 1 - ieta ) : ( ieta - nFullEta / 2 );
                AppendJerLines( buf, kEtaEdges[ieta], kEtaEdges[ieta + 1], etaRecords[iAbs] );
            }

            TString tmpPath = path + ".tmp";
            { std::ofstream tmp( tmpPath.Data() ); tmp << buf.str(); }
            bool ok = true;
            try {
                JME::JetResolutionObject obj( tmpPath.Data() );
                obj.saveToFile( path.Data() );
            } catch( const std::exception& ex ){
                std::cerr << "ERROR writing pT resolution file " << path << ": " << ex.what() << "\n";
                ok = false;
            }
            gSystem->Unlink( tmpPath );
            if( ok ) std::cout << "Done (pT resolution). " << cone << ": -> " << path << "\n";
        }
    }

    fIn->Close();
}
