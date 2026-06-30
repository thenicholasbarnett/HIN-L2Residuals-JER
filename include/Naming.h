#ifndef NAMING_H
#define NAMING_H

#include "TString.h"
#include "Binning.h"

#include <vector>

namespace L2Name {

inline TString CleanKey( const TString& s ){
    return s.BeginsWith( "_" ) ? s( 1, s.Length() - 1 ) : s;
}

inline TString EtaKey( int ieta ){
    return Form( "eta%02d", ieta );
}

inline TString PtKey( const RangeBin& ptSlice ){
    return CleanKey( ptSlice.shortName );
}

inline TString AlphaKey( const RangeBin& alphaSlice ){
    return CleanKey( alphaSlice.shortName );
}

inline TString EtaModeKey( bool fullEta ){
    return fullEta ? "fulleta" : "abseta";
}

inline TString Join( const std::vector<TString>& parts, const TString& sep = "_" ){
    TString out;
    for( const auto& part : parts ){
        if( part.IsNull() ) continue;
        if( !out.IsNull() ) out += sep;
        out += part;
    }
    return out;
}

inline TString ObjectName( const TString& cone,
                          const TString& objectKind,
                          const std::vector<TString>& orderedKeys = {},
                          const std::vector<TString>& detailKeys = {} ){
    std::vector<TString> parts = { cone, objectKind };
    parts.insert( parts.end(), orderedKeys.begin(), orderedKeys.end() );
    parts.insert( parts.end(), detailKeys.begin(), detailKeys.end() );
    return Join( parts );
}

}

#endif
