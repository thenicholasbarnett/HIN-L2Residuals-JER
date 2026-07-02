#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"

#include "Binning.h"
#include "Naming.h"
#include "TextFileWriter.h"

static int nPass = 0;
static int nFail = 0;
static constexpr double kEps = 1e-9;

void Check( bool cond, const char* msg ){
    if( cond ){ std::cout << "  PASS  " << msg << std::endl; nPass++; }
    else { std::cout << "  FAIL  " << msg << std::endl; nFail++; }
}

// ---- helpers ----

// Build a synthetic residuals ROOT file for cone "ak4PF", populating only the
// given pT-slice indices (into BinningConfig::ptavgSlices) with intercept
// histograms — both abseta and fulleta — set to corrValue ± corrErr in every
// bin. The ak4PF TDirectory always exists, even if sliceIndices is empty.
// norm appends "_norm" to the histogram names; recreate=false appends to an
// existing file (opened "update") instead of overwriting it, so a norm
// variant can be added to a file that already has the direct variant.
static TString MakeResiduals( const char* path, double corrValue, double corrErr,
                              const std::vector<int>& sliceIndices,
                              bool norm = false, bool recreate = true ){
    TFile* f = new TFile( path, recreate ? "recreate" : "update" );
    BinningConfig bins;
    TDirectory* coneDir = recreate ? f->mkdir( "ak4PF" ) : ( TDirectory* )f->Get( "ak4PF" );
    coneDir->cd();
    const TString suffix = norm ? "_norm" : "";
    for( int ip : sliceIndices ){
        for( int em = 0; em < 2; em++ ){
            const bool fullEta = ( em == 1 );
            const std::vector<Double_t>& edges = fullEta ? kEtaEdges : kAbsEtaEdges;
            const int nEta = ( int )edges.size() - 1;
            TString name = L2Name::ObjectName( "ak4PF", "intercept",
                { L2Name::EtaModeKey( fullEta ), L2Name::PtKey( bins.ptavgSlices[ip] ) }, { "gauss" } ) + suffix;
            TH1D* h = new TH1D( name, "", nEta, edges.data() );
            for( int ieta = 1; ieta <= nEta; ieta++ ){
                h->SetBinContent( ieta, corrValue );
                h->SetBinError( ieta, corrErr );
            }
            h->Write();
            delete h;
        }
    }
    f->Close();
    delete f;
    return TString( path );
}

// The real cfg.hltJ80Thresh (cfg/2024ppRef.toml) is 100.0: pT slices 0,1
// (30-70, 70-100) fall below threshold -> NonTriggered; slices 2-5 (100-175 ... 500-1000)
// are at/above -> Triggered. Build both files with full slice coverage, possibly at
// different values, so tests can check which source ends up in the merge.
static void MakeTrigNoTrig( const char* trigPath, const char* notrigPath,
                      double trigValue, double notrigValue, double err = 0.001 ){
    MakeResiduals( trigPath, trigValue, err, {0, 1, 2, 3, 4, 5} );
    MakeResiduals( notrigPath, notrigValue, err, {0, 1, 2, 3, 4, 5} );
}

// runTextFile now writes correction text files to a fixed location
// (data/jec/preliminary/, relative to the repo root) regardless of PREFIX --
// PREFIX is just a plain filename prefix now, not a path. Tests still use
// distinct prefixes per case so cleanup doesn't collide across cases.
static TString TxtPath( const TString& prefix, const char* mode, bool norm ){
    return TString( "data/jec/preliminary/" ) + prefix + "_ak4PF_" + mode + ( norm ? "_norm" : "" ) + ".txt";
}

static void CleanupFiles( const char* trigPath, const char* notrigPath,
                         const char* rootPath, const TString& prefix ){
    std::remove( trigPath );
    std::remove( notrigPath );
    std::remove( rootPath );
    std::remove( TxtPath( prefix, "abseta", false ).Data() );
    std::remove( TxtPath( prefix, "eta", false ).Data() );
    std::remove( TxtPath( prefix, "abseta", true ).Data() );
    std::remove( TxtPath( prefix, "eta", true ).Data() );
}

// Parse one data line from the JEC text file.
// Returns false if parsing fails or field count is wrong.
struct JECLine {
    double etaLo, etaHi;
    int npar;
    double ptLo, ptHi;
    double p[3];
};

static bool ParseJECLine( const std::string& line, JECLine& out ){
    std::istringstream ss( line );
    if( !( ss >> out.etaLo >> out.etaHi >> out.npar >> out.ptLo >> out.ptHi )) return false;
    for( int i = 0; i < 3; i++ ) if( !( ss >> out.p[i] )) return false;
    return true;
}

// Read all JEC data lines (skip the header) from a text file.
static std::vector<JECLine> ReadJECFile( const char* path, std::string& header ){
    std::ifstream f( path );
    std::vector<JECLine> lines;
    bool firstLine = true;
    std::string line;
    while( std::getline( f, line )){
        if( firstLine ){ header = line; firstLine = false; continue; }
        if( line.empty() ) continue;
        JECLine jl;
        if( ParseJECLine( line, jl )) lines.push_back( jl );
    }
    return lines;
}

// ---- test cases ----

// [1] Both text files exist with correct header and line counts
void TestFileStructure(){
    std::cout << "\n[1] File structure\n";

    const char* trigPath = "/tmp/tw_trig1.root";
    const char* notrigPath = "/tmp/tw_notrig1.root";
    const char* rootPath = "/tmp/tw_out1.root";
    const char* prefix = "test_tw1";
    CleanupFiles( trigPath, notrigPath, rootPath, prefix );

    MakeTrigNoTrig( trigPath, notrigPath, 1.0, 1.0 );
    runTextFile( trigPath, notrigPath, rootPath, prefix, "gauss", false );

    std::string header, headerEta;
    auto absLines = ReadJECFile( TxtPath( prefix, "abseta", false ).Data(), header );
    auto etaLines = ReadJECFile( TxtPath( prefix, "eta", false ).Data(), headerEta );

    Check( header.rfind( "{1 JetEta", 0 ) == 0, "abseta header starts with {1 JetEta" );
    Check( header.find( "L2Residual" ) != std::string::npos, "abseta header contains L2Residual" );
    Check( ( int )absLines.size() == 36, "abseta file has 36 data lines" );
    Check( headerEta == header, "eta file header matches abseta header" );
    Check( ( int )etaLines.size() == 36, "eta file has 36 data lines" );

    if( !absLines.empty() ){
        Check( absLines.front().npar == 5, "Npar = 5" );
        Check( std::fabs( absLines.front().ptLo - 40.0 ) < kEps, "pT_lo = 40" );
        Check( std::fabs( absLines.front().ptHi - 1000.0 ) < kEps, "pT_hi = 1000" );
    }

    CleanupFiles( trigPath, notrigPath, rootPath, prefix );
}

// [2] Merge source selection: Triggered and NonTriggered carry different values; the merged
// corrfinal grid must pull each pT slice from the correct source.
void TestMergeSourceSelection(){
    std::cout << "\n[2] Merge source selection\n";

    const char* trigPath = "/tmp/tw_trig2.root";
    const char* notrigPath = "/tmp/tw_notrig2.root";
    const char* rootPath = "/tmp/tw_out2.root";
    const char* prefix = "test_tw2";
    CleanupFiles( trigPath, notrigPath, rootPath, prefix );

    MakeTrigNoTrig( trigPath, notrigPath, 5.0, 1.0 );
    runTextFile( trigPath, notrigPath, rootPath, prefix, "gauss", false );

    TFile* fOut = TFile::Open( rootPath, "read" );
    TH2D* hGrid = fOut ? ( TH2D* )fOut->Get( "ak4PF/ak4PF_corrfinal_abseta_gauss" ) : nullptr;
    Check( hGrid != nullptr, "corrfinal_abseta_gauss TH2D exists in output" );
    if( hGrid ){
        // y-bin 1 = ptavg_30_70  (lo=30  < 100 threshold) -> NonTriggered value
        Check( std::fabs( hGrid->GetBinContent( 1, 1 ) - 1.0 ) < 1e-6, "ptavg_30_70 pulled from NonTriggered" );
        // y-bin 3 = ptavg_100_175 (lo=100 >= 100 threshold) -> Triggered value
        Check( std::fabs( hGrid->GetBinContent( 1, 3 ) - 5.0 ) < 1e-6, "ptavg_100_175 pulled from Triggered" );
    }
    if( fOut ) fOut->Close();

    CleanupFiles( trigPath, notrigPath, rootPath, prefix );
}

// [2b] useNorm selects the "_norm"-suffixed intercepts and suffixes the
// output objects/filenames accordingly, leaving the direct variant untouched.
void TestNormSelection(){
    std::cout << "\n[2b] Normalized variant selection\n";

    const char* trigPath = "/tmp/tw_trig2b.root";
    const char* notrigPath = "/tmp/tw_notrig2b.root";
    const char* rootPath = "/tmp/tw_out2b.root";
    const char* prefix = "test_tw2b";
    CleanupFiles( trigPath, notrigPath, rootPath, prefix );

    // direct variant = 1.0, norm variant = 9.0, both present in the same file
    MakeResiduals( trigPath, 1.0, 0.001, {0, 1, 2, 3, 4, 5}, false, true );
    MakeResiduals( trigPath, 9.0, 0.001, {0, 1, 2, 3, 4, 5}, true, false );
    MakeResiduals( notrigPath, 1.0, 0.001, {0, 1, 2, 3, 4, 5}, false, true );
    MakeResiduals( notrigPath, 9.0, 0.001, {0, 1, 2, 3, 4, 5}, true, false );

    runTextFile( trigPath, notrigPath, rootPath, prefix, "gauss", true );

    Check( std::ifstream( TxtPath( prefix, "abseta", true ).Data() ).good(),
        "norm run writes a _norm-suffixed abseta text file" );
    Check( std::ifstream( TxtPath( prefix, "eta", true ).Data() ).good(),
        "norm run writes a _norm-suffixed eta text file" );
    Check( !std::ifstream( TxtPath( prefix, "abseta", false ).Data() ).good(),
        "norm run does not also write the direct-variant filename" );

    TFile* fOut = TFile::Open( rootPath, "read" );
    TH2D* hGrid = fOut ? ( TH2D* )fOut->Get( "ak4PF/ak4PF_corrfinal_abseta_gauss_norm" ) : nullptr;
    Check( hGrid != nullptr, "corrfinal_abseta_gauss_norm TH2D exists in output" );
    if( hGrid ) Check( std::fabs( hGrid->GetBinContent( 1, 1 ) - 9.0 ) < 1e-6,
        "grid picked up the norm value (9.0), not the direct value (1.0)" );
    if( fOut ) fOut->Close();

    CleanupFiles( trigPath, notrigPath, rootPath, prefix );
}

// [2c] Single-file overload: no Triggered/NonTriggered split, every pT slice reads from the
// one residuals file regardless of cfg.hltJ80Thresh.
void TestSingleFileMode(){
    std::cout << "\n[2c] Single-file mode (no Triggered/NonTriggered split)\n";

    const char* resPath = "/tmp/tw_single.root";
    const char* rootPath = "/tmp/tw_out_single.root";
    const char* prefix = "test_tw_single";
    CleanupFiles( resPath, resPath, rootPath, prefix );

    MakeResiduals( resPath, 3.0, 0.001, {0, 1, 2, 3, 4, 5} );
    runTextFile( resPath, rootPath, prefix, "gauss", false );

    Check( std::ifstream( TxtPath( prefix, "abseta", false ).Data() ).good(),
        "single-file run writes the abseta text file" );

    TFile* fOut = TFile::Open( rootPath, "read" );
    TH2D* hGrid = fOut ? ( TH2D* )fOut->Get( "ak4PF/ak4PF_corrfinal_abseta_gauss" ) : nullptr;
    Check( hGrid != nullptr, "corrfinal_abseta_gauss TH2D exists in output" );
    if( hGrid ){
        // y-bin 1 = ptavg_30_70 (below threshold, would be NonTriggered in two-file mode);
        // y-bin 3 = ptavg_100_175 (at/above threshold, would be Triggered) — both must
        // come from the single file regardless, so both equal 3.0.
        Check( std::fabs( hGrid->GetBinContent( 1, 1 ) - 3.0 ) < 1e-6,
            "low-pT slice reads from the single file" );
        Check( std::fabs( hGrid->GetBinContent( 1, 3 ) - 3.0 ) < 1e-6,
            "high-pT slice also reads from the single file" );
    }
    if( fOut ) fOut->Close();

    CleanupFiles( resPath, resPath, rootPath, prefix );
}

// [3] Eta ordering in the abseta (mirrored) text file
void TestEtaOrdering(){
    std::cout << "\n[3] Eta ordering (abseta, mirrored)\n";

    const char* trigPath = "/tmp/tw_trig3.root";
    const char* notrigPath = "/tmp/tw_notrig3.root";
    const char* rootPath = "/tmp/tw_out3.root";
    const char* prefix = "test_tw3";
    CleanupFiles( trigPath, notrigPath, rootPath, prefix );

    MakeTrigNoTrig( trigPath, notrigPath, 1.0, 1.0 );
    runTextFile( trigPath, notrigPath, rootPath, prefix, "gauss", false );

    std::string header;
    auto lines = ReadJECFile( TxtPath( prefix, "abseta", false ).Data(), header );
    if( ( int )lines.size() != 36 ){
        std::cout << "  SKIP  (line count wrong)\n";
    } else {
        Check( lines[0].etaLo < -5.0, "first line etaLo < -5" );
        Check( lines[0].etaHi < 0.0, "first line etaHi < 0" );
        Check( lines[17].etaLo < 0.0, "line 18 etaLo < 0 (inner negative)" );
        Check( std::fabs( lines[17].etaHi ) < kEps, "line 18 etaHi = 0 (boundary)" );
        Check( std::fabs( lines[18].etaLo ) < kEps, "line 19 etaLo = 0 (boundary)" );
        Check( lines[18].etaHi > 0.0, "line 19 etaHi > 0 (inner positive)" );
        Check( lines[35].etaLo > 0.0, "last line etaLo > 0" );
        Check( lines[35].etaHi > 5.0, "last line etaHi > 5" );

        bool negOK = true;
        for( int i = 0; i < 18; i++ )
            if( !( lines[i].etaLo < lines[i].etaHi && lines[i].etaHi <= 0.0 )) negOK = false;
        Check( negOK, "all negative-eta lines have etaLo < etaHi <= 0" );

        bool posOK = true;
        for( int i = 18; i < 36; i++ )
            if( !( lines[i].etaLo >= 0.0 && lines[i].etaLo < lines[i].etaHi )) posOK = false;
        Check( posOK, "all positive-eta lines have 0 <= etaLo < etaHi" );

        bool mirrorOK = true;
        for( int i = 0; i < 18; i++ ){
            double absLoNeg = -lines[17 - i].etaHi;
            double absHiNeg = -lines[17 - i].etaLo;
            double absloPOS = lines[18 + i].etaLo;
            double abshiPOS = lines[18 + i].etaHi;
            if( std::fabs( absLoNeg - absloPOS ) > kEps || std::fabs( absHiNeg - abshiPOS ) > kEps )
                mirrorOK = false;
        }
        Check( mirrorOK, "negative and positive eta bin edges are symmetric" );

        bool paramsMatch = true;
        for( int i = 0; i < 18; i++ ){
            const JECLine& neg = lines[17 - i];
            const JECLine& pos = lines[18 + i];
            for( int p = 0; p < 3; p++ )
                if( std::fabs( neg.p[p] - pos.p[p] ) > 1e-12 ) paramsMatch = false;
        }
        Check( paramsMatch, "negative and positive sides have identical fit parameters (mirrored)" );
    }

    CleanupFiles( trigPath, notrigPath, rootPath, prefix );
}

// [4] Eta ordering in the eta (independent, unmirrored) text file
void TestEtaFileOrdering(){
    std::cout << "\n[4] Eta ordering (full eta, independent fits)\n";

    const char* trigPath = "/tmp/tw_trig4.root";
    const char* notrigPath = "/tmp/tw_notrig4.root";
    const char* rootPath = "/tmp/tw_out4.root";
    const char* prefix = "test_tw4";
    CleanupFiles( trigPath, notrigPath, rootPath, prefix );

    MakeTrigNoTrig( trigPath, notrigPath, 1.0, 1.0 );
    runTextFile( trigPath, notrigPath, rootPath, prefix, "gauss", false );

    std::string header;
    auto lines = ReadJECFile( TxtPath( prefix, "eta", false ).Data(), header );
    if( ( int )lines.size() != 36 ){
        std::cout << "  SKIP  (line count wrong)\n";
    } else {
        Check( std::fabs( lines.front().etaLo - ( -5.191 )) < kEps, "first line etaLo = -5.191" );
        Check( std::fabs( lines.back().etaHi - 5.191 ) < kEps, "last line etaHi = 5.191" );

        bool ascending = true;
        for( int i = 0; i < 36; i++ )
            if( !( lines[i].etaLo < lines[i].etaHi )) ascending = false;
        for( int i = 1; i < 36; i++ )
            if( !( std::fabs( lines[i].etaLo - lines[i - 1].etaHi ) < kEps )) ascending = false;
        Check( ascending, "36 lines are contiguous and ascending in eta, no mirroring" );
    }

    CleanupFiles( trigPath, notrigPath, rootPath, prefix );
}

// [5] Unity fallback: fewer than kMinSlices (3) valid pT slices across the merge
void TestUnityFallback_TooFewSlices(){
    std::cout << "\n[5] Unity fallback — fewer than kMinSlices pT slices\n";

    const char* trigPath = "/tmp/tw_trig5.root";
    const char* notrigPath = "/tmp/tw_notrig5.root";
    const char* rootPath = "/tmp/tw_out5.root";
    const char* prefix = "test_tw5";
    CleanupFiles( trigPath, notrigPath, rootPath, prefix );

    // Only 2 valid slices total: Triggered supplies 100-175 and 175-250, NonTriggered supplies none.
    MakeResiduals( trigPath, 1.0, 0.001, {2, 3} );
    MakeResiduals( notrigPath, 1.0, 0.001, {} );
    runTextFile( trigPath, notrigPath, rootPath, prefix, "gauss", false );

    std::string header;
    auto lines = ReadJECFile( TxtPath( prefix, "abseta", false ).Data(), header );
    if( ( int )lines.size() != 36 ){
        std::cout << "  SKIP  (line count wrong)\n";
    } else {
        bool allUnity = true;
        for( const auto& l : lines )
            if( std::fabs( l.p[0] - 1.0 ) > kEps || std::fabs( l.p[1] ) > kEps || std::fabs( l.p[2] ) > kEps )
                allUnity = false;
        Check( allUnity, "all bins fall back to unity (p0=1, p1=0, p2=0)" );
    }

    CleanupFiles( trigPath, notrigPath, rootPath, prefix );
}

// [6] Unity fallback for empty histograms (all bins zero)
void TestUnityFallback_EmptyBins(){
    std::cout << "\n[6] Unity fallback — empty histograms\n";

    const char* trigPath = "/tmp/tw_trig6.root";
    const char* notrigPath = "/tmp/tw_notrig6.root";
    const char* rootPath = "/tmp/tw_out6.root";
    const char* prefix = "test_tw6";
    CleanupFiles( trigPath, notrigPath, rootPath, prefix );

    MakeTrigNoTrig( trigPath, notrigPath, 0.0, 0.0, 0.0 );
    runTextFile( trigPath, notrigPath, rootPath, prefix, "gauss", false );

    std::string header;
    auto lines = ReadJECFile( TxtPath( prefix, "abseta", false ).Data(), header );
    if( ( int )lines.size() != 36 ){
        std::cout << "  SKIP  (line count wrong)\n";
    } else {
        bool allUnity = true;
        for( const auto& l : lines )
            if( std::fabs( l.p[0] - 1.0 ) > kEps || std::fabs( l.p[1] ) > kEps || std::fabs( l.p[2] ) > kEps )
                allUnity = false;
        Check( allUnity, "empty bins produce unity correction (p0=1, p1=0, p2=0)" );
    }

    CleanupFiles( trigPath, notrigPath, rootPath, prefix );
}

// [7] Fit round-trip: unit input recovers f≈1.0 at all pT centers
void TestFitRoundTrip(){
    std::cout << "\n[7] Fit round-trip\n";

    const char* trigPath = "/tmp/tw_trig7.root";
    const char* notrigPath = "/tmp/tw_notrig7.root";
    const char* rootPath = "/tmp/tw_out7.root";
    const char* prefix = "test_tw7";
    CleanupFiles( trigPath, notrigPath, rootPath, prefix );

    MakeTrigNoTrig( trigPath, notrigPath, 1.0, 1.0, 1e-5 );
    runTextFile( trigPath, notrigPath, rootPath, prefix, "gauss", false );

    std::string header;
    auto lines = ReadJECFile( TxtPath( prefix, "abseta", false ).Data(), header );
    if( ( int )lines.size() != 36 ){
        std::cout << "  SKIP  (line count wrong)\n";
    } else {
        // slice midpoints: 30-70, 70-100, 100-175, 175-250, 250-500, 500-1000
        const double ptCenters[] = { 50.0, 85.0, 137.5, 212.5, 375.0, 750.0 };
        const JECLine& barrel = lines[18];  // innermost positive eta
        bool fitsClose = true;
        for( double pT : ptCenters ){
            double denom = barrel.p[0] + barrel.p[1] * std::log10( 0.01 * pT ) + barrel.p[2] / ( pT / 10.0 );
            double corr = ( denom != 0.0 ) ? 1.0 / denom : 0.0;
            if( std::fabs( corr - 1.0 ) > 0.01 ){ fitsClose = false; break; }
        }
        Check( fitsClose, "fit recovers f≈1.0 at all pT centers for unit input" );
    }

    CleanupFiles( trigPath, notrigPath, rootPath, prefix );
}

// [8] CMS JEC edge value: outermost eta bin ends at 5.191
void TestEtaExtent(){
    std::cout << "\n[8] CMS eta extent\n";

    const char* trigPath = "/tmp/tw_trig8.root";
    const char* notrigPath = "/tmp/tw_notrig8.root";
    const char* rootPath = "/tmp/tw_out8.root";
    const char* prefix = "test_tw8";
    CleanupFiles( trigPath, notrigPath, rootPath, prefix );

    MakeTrigNoTrig( trigPath, notrigPath, 1.0, 1.0 );
    runTextFile( trigPath, notrigPath, rootPath, prefix, "gauss", false );

    std::string header;
    auto lines = ReadJECFile( TxtPath( prefix, "abseta", false ).Data(), header );
    if( ( int )lines.size() != 36 ){
        std::cout << "  SKIP  (line count wrong)\n";
    } else {
        Check( std::fabs( lines[0].etaLo - ( -5.191 )) < kEps, "negative outer edge = -5.191" );
        Check( std::fabs( lines[35].etaHi - 5.191 ) < kEps, "positive outer edge =  5.191" );
    }

    CleanupFiles( trigPath, notrigPath, rootPath, prefix );
}

// ---- main ----

int main(){
    gROOT->SetBatch( true );

    std::cout << "=== TestTextFileWriter ===\n";

    TestFileStructure();
    TestMergeSourceSelection();
    TestNormSelection();
    TestSingleFileMode();
    TestEtaOrdering();
    TestEtaFileOrdering();
    TestUnityFallback_TooFewSlices();
    TestUnityFallback_EmptyBins();
    TestFitRoundTrip();
    TestEtaExtent();

    std::cout << "\n=== " << nPass << " passed, " << nFail << " failed ===\n";
    return nFail > 0 ? 1 : 0;
}
