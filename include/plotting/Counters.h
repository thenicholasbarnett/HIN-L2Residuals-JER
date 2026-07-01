#ifndef L2RESIDUALS_PLOTTING_COUNTERS_H
#define L2RESIDUALS_PLOTTING_COUNTERS_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "plotting/AsymmetryDistributions.h"   // kMinEntriesPlot
#include "plotting/Kinematics.h"               // kNKinematicsCollections, kKinematicsCollections, kNKinematicsPtMins
#include "Binning.h"
#include "Naming.h"

// ============================================================
// Progress accounting
// ============================================================

inline int CountAsymDistPlots( TFile* fIn, const TString& cone, const BinningConfig& bins ){
    TDirectory* dData = ( TDirectory* )fIn->Get( cone + "/QA_data" );
    TDirectory* dMC   = ( TDirectory* )fIn->Get( cone + "/QA_mc" );
    if( !dData || !dMC ) return 0;

    int n = 0;
    const int nEta = ( int )kAbsEtaEdges.size() - 1;
    for( const auto& ptSl : bins.ptavgSlices ){
        const TString ptKey = L2Name::PtKey( ptSl );
        for( const auto& aSl : bins.alphaSlices ){
            const TString alphaKey = L2Name::AlphaKey( aSl );
            for( int ie = 0; ie < nEta; ie++ ){
                const TString etaKey = L2Name::EtaKey( ie, false );
                const TString dname = L2Name::ObjectName( cone, "A_data",
                    {L2Name::EtaModeKey( false ), etaKey, ptKey, alphaKey} );
                const TString mname = L2Name::ObjectName( cone, "A_mc",
                    {L2Name::EtaModeKey( false ), etaKey, ptKey, alphaKey} );
                TH1D* hd = ( TH1D* )dData->Get( dname );
                TH1D* hm = ( TH1D* )dMC->Get( mname );
                if( hd && hm && hd->GetEntries() >= kMinEntriesPlot ) n += 3;
            }
        }
    }
    return n;
}

inline int CountROverlayPlots( TFile* fIn, const TString& cone, const BinningConfig& bins ){
    TDirectory* dRvals = ( TDirectory* )fIn->Get( cone + "/Rvals" );
    if( !dRvals ) return 0;

    int n = 0;
    for( int m = 0; m < kNMethods; m++ ){
        for( const auto& aSl : bins.alphaSlices ){
            const TString alphaKey = L2Name::AlphaKey( aSl );
            for( const auto& ptSl : bins.ptavgSlices ){
                const TString ptKey = L2Name::PtKey( ptSl );
                const TString rdName = L2Name::ObjectName( cone, "R_data",
                    {L2Name::EtaModeKey( false ), ptKey, alphaKey}, {kMethodKeys[m]} );
                const TString rmName = L2Name::ObjectName( cone, "R_mc",
                    {L2Name::EtaModeKey( false ), ptKey, alphaKey}, {kMethodKeys[m]} );
                if( HasH( dRvals, rdName ) && HasH( dRvals, rmName ) ) n++;
            }
        }
    }
    return n;
}

inline int CountAlphaFitPlots( TFile* fIn, const TString& cone, const BinningConfig& bins ){
    TDirectory* dGraphs = ( TDirectory* )fIn->Get( cone + "/graphs" );
    if( !dGraphs ) return 0;

    int n = 0;
    const int nEta = ( int )kAbsEtaEdges.size() - 1;
    for( int m = 0; m < kNMethods; m++ ){
        for( const auto& ptSl : bins.ptavgSlices ){
            const TString ptKey = L2Name::PtKey( ptSl );
            for( int ie = 0; ie < nEta; ie++ ){
                const TString etaKey = L2Name::EtaKey( ie, false );
                const TString gname = L2Name::ObjectName( cone, "R",
                    {L2Name::EtaModeKey( false ), etaKey, ptKey}, {kMethodKeys[m]} );
                if( HasGraphWithN( dGraphs, {gname + "_norm", gname}, 2 ) ) n++;
            }
        }
    }
    return n;
}

inline int CountPtFitPlots( TFile* fIn, const TString& cone ){
    TDirectory* dGraphs = ( TDirectory* )fIn->Get( cone + "/graphs" );
    if( !dGraphs ) return 0;

    int n = 0;
    for( int m = 0; m < kNMethods; m++ ){
        for( int mode = 0; mode < 2; mode++ ){
            const bool fullEta = ( mode == 1 );
            const int nEta = fullEta ? ( int )kEtaEdges.size() - 1 : ( int )kAbsEtaEdges.size() - 1;
            const TString etaMode = L2Name::EtaModeKey( fullEta );
            for( int ie = 0; ie < nEta; ie++ ){
                const TString etaKey = L2Name::EtaKey( ie, fullEta );
                const TString gname = L2Name::ObjectName( cone, "ptcorr",
                    {etaMode, etaKey}, {kMethodKeys[m]} );
                if( HasGraphWithN( dGraphs, {gname + "_norm", gname}, 2 ) ) n++;
            }
        }
    }
    return n;
}

inline int CountEtaSymPlots( TFile* fIn, const TString& cone, const BinningConfig& bins ){
    int n = 0;
    for( int m = 0; m < kNMethods; m++ ){
        for( const auto& sl : bins.ptavgSlices ){
            const TString ptKey = L2Name::PtKey( sl );
            const TString nameAbs = L2Name::ObjectName( cone, "intercept",
                {L2Name::EtaModeKey( false ), ptKey}, {kMethodKeys[m]} );
            const TString nameFull = L2Name::ObjectName( cone, "intercept",
                {L2Name::EtaModeKey( true ), ptKey}, {kMethodKeys[m]} );
            if( HasHAny( fIn, {cone + "/" + nameAbs} ) && HasHAny( fIn, {cone + "/" + nameFull} ) ) n++;
        }
    }
    return n;
}

inline int CountMethodCompPlots( TFile* fIn, const TString& cone, const BinningConfig& bins, bool fullEta ){
    int n = 0;
    const TString etaMode = L2Name::EtaModeKey( fullEta );
    for( const auto& sl : bins.ptavgSlices ){
        const TString name = L2Name::ObjectName( cone, "intercept",
            {etaMode, L2Name::PtKey( sl )}, {kMethodKeys[0]} );
        if( HasHAny( fIn, {cone + "/" + name} ) ) n++;
    }
    return n;
}

inline int CountFinalsPlots( TFile* fIn, const TString& cone, const BinningConfig& bins ){
    int n = 0;
    for( int m = 0; m < kNMethods; m++ ){
        for( int mode = 0; mode < 2; mode++ ){
            const bool fullEta = ( mode == 1 );
            const TString etaMode = L2Name::EtaModeKey( fullEta );

            bool anyValid = false;
            for( const auto& ptSl : bins.ptavgSlices ){
                const TString name = L2Name::ObjectName( cone, "intercept",
                    {etaMode, L2Name::PtKey( ptSl )}, {kMethodKeys[m]} );
                if( HasHAny( fIn, {cone + "/" + name} ) ){ anyValid = true; break; }
            }
            if( !anyValid ){
                const TString gridName = L2Name::ObjectName( cone, "corrfinal", {etaMode}, {kMethodKeys[m]} );
                anyValid = HasHAny( fIn, {cone + "/" + gridName + "_norm", gridName + "_norm",
                                          cone + "/" + gridName, gridName} );
            }
            if( anyValid ) n++;
        }
    }
    return n;
}

inline int CountNormCompPlots( TFile* fIn, const TString& cone, const BinningConfig& bins, bool fullEta ){
    int n = 0;
    const TString etaMode = L2Name::EtaModeKey( fullEta );
    for( int m = 0; m < kNMethods; m++ ){
        bool anyValid = false;
        for( const auto& sl : bins.ptavgSlices ){
            const TString ptKey = L2Name::PtKey( sl );
            const TString nameDirect = L2Name::ObjectName( cone, "intercept",
                {etaMode, ptKey}, {kMethodKeys[m]} );
            if( HasHAny( fIn, {cone + "/" + nameDirect} )
                && HasHAny( fIn, {cone + "/" + nameDirect + "_norm"} ) ){
                anyValid = true;
                break;
            }
        }
        if( anyValid ) n++;
    }
    return n;
}

inline int CountKinematicsPlots( TFile* fIn, const TString& cone ){
    if( !fIn->Get( cone + "/" + cone + "_incl" ) ) return 0;

    const int plotsPerCollection = 3 + kNKinematicsPtMins;
    int n = 0;
    bool anyCollection = false;
    for( int ic = 0; ic < kNKinematicsCollections; ic++ ){
        const TString coll = kKinematicsCollections[ic].inputKey;
        if( fIn->Get( cone + "/" + cone + "_" + coll ) ){
            n += plotsPerCollection;
            anyCollection = true;
        }
    }
    if( anyCollection ) n += 3;
    return n;
}

inline int CountEventPlots( TFile* fIn ){
    int n = 0;
    if( fIn->Get( "hvz_all" ) && fIn->Get( "hvz" ) ) n++;
    if( fIn->Get( "hfilt" ) ) n++;
    if( fIn->Get( "h_hlt_j80" ) ) n++;
    return n;
}

#endif
