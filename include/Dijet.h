#ifndef DIJET_H
#define DIJET_H

#include "Rtypes.h"
#include "TMath.h"

#include <cmath>

static constexpr float kBarrelEtaCut = 1.3f;
static constexpr float kDphiCut = 2.7f;
static constexpr float kSubleadPtCut = 10.0f;

struct SortedJets {
    int lead = -1;
    int sublead = -1;
    int third = -1;
};

inline SortedJets FindLeadingJets( const float* pt, int n ){
    SortedJets s;
    for( int i = 0; i < n; i++ ){
        if( s.lead == -1 || pt[i] > pt[s.lead] ){
            s.third = s.sublead;
            s.sublead = s.lead;
            s.lead = i;
        }
        else if( s.sublead == -1 || pt[i] > pt[s.sublead] ){
            s.third = s.sublead;
            s.sublead = i;
        }
        else if( s.third == -1 || pt[i] > pt[s.third] ){
            s.third = i;
        }
    }
    return s;
}

struct DijetResult {
    bool valid = false;
    int tagIdx = -1;
    int probeIdx = -1;
    float ptavg = 0.0f;
    float A = 0.0f;
    float alpha = 0.0f;
    float dphi = 0.0f;
};

inline float DijetDPhi( float phi1, float phi2 ){
    float d = std::fabs( phi1 - phi2 );
    if( d > ( float )TMath::Pi() ){ d = 2.0f * ( float )TMath::Pi() - d; }
    return d;
}

inline DijetResult MakeDijet( SortedJets s, bool hasThird,
                              const float* pt, const float* eta, const float* phi,
                              ULong64_t eventNumber,
                              float subleadPtCut = kSubleadPtCut,
                              float dphiCut = kDphiCut ){
    DijetResult d;

    bool leadBarrel = std::fabs( eta[s.lead] ) < kBarrelEtaCut;
    bool subleadBarrel = std::fabs( eta[s.sublead] ) < kBarrelEtaCut;

    if( !leadBarrel && !subleadBarrel ){ return d; }
    if( pt[s.sublead] < subleadPtCut ){ return d; }

    float dphi = DijetDPhi( phi[s.lead], phi[s.sublead] );
    if( dphi < dphiCut ){ return d; }

    int tagIdx, probeIdx;
    if( leadBarrel && !subleadBarrel ){
        tagIdx = s.lead;
        probeIdx = s.sublead;
    }
    else if( !leadBarrel && subleadBarrel ){
        tagIdx = s.sublead;
        probeIdx = s.lead;
    }
    else {
        if( eventNumber % 2 == 0 ){
            tagIdx = s.lead;
            probeIdx = s.sublead;
        }
        else {
            tagIdx = s.sublead;
            probeIdx = s.lead;
        }
    }

    float ptsum = pt[tagIdx] + pt[probeIdx];

    d.valid = true;
    d.tagIdx = tagIdx;
    d.probeIdx = probeIdx;
    d.ptavg = 0.5f * ptsum;
    d.A = ( pt[probeIdx] - pt[tagIdx] ) / ptsum;
    d.alpha = hasThird ? pt[s.third] / d.ptavg : 0.0f;
    d.dphi = dphi;

    return d;
}

#endif
