#ifndef JETSMEARING_H
#define JETSMEARING_H

#include "jetmet_jer/JetResolution.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <string>

// CMS hybrid-method JER smearing (mirrors PhysicsTools/PatUtils's
// SmearedJetProducerT.h, stripped of the EDProducer/EventSetup scaffolding).
// Not yet wired into any pipeline step -- see JER smearing follow-up.
namespace JetSmearing {

struct Result {
  double smearFactor = 1.0;
  double resolution = 0.0; // sigma_JER actually used (getResolution() output)
  double scaleFactor = 1.0;
  bool matched = false;
};

// genPt < 0 means "no matched gen jet", the same sentinel convention already
// used for ref.pt in ConeHistograms' response matching -- this function does
// not re-derive a dR match itself, it trusts the caller's existing gen match.
inline Result ComputeSmearFactor(double recoPt, double eta, double rho,
                                 double genPt,
                                 const JME::JetResolution &resolution,
                                 const JME::JetResolutionScaleFactor &resolutionSF,
                                 std::mt19937 &rng,
                                 Variation variation = Variation::NOMINAL,
                                 const std::string &uncertaintySource = "",
                                 double dPtMaxFactor = 3.0) {
  Result r;
  r.resolution = resolution.getResolution(
      JME::JetParameters().setJetPt(recoPt).setJetEta(eta).setRho(rho));
  r.scaleFactor = resolutionSF.getScaleFactor(
      JME::JetParameters().setJetPt(recoPt).setJetEta(eta), variation,
      uncertaintySource);

  if (genPt >= 0 &&
      std::abs(recoPt - genPt) < dPtMaxFactor * r.resolution * recoPt) {
    // scaling method: gen-matched jets are corrected deterministically
    r.matched = true;
    r.smearFactor = 1.0 + (r.scaleFactor - 1.0) * (recoPt - genPt) / recoPt;
  } else if (r.scaleFactor > 1.0) {
    // stochastic method: no usable gen match, widen resolution by hand
    double sigma = r.resolution * std::sqrt(r.scaleFactor * r.scaleFactor - 1.0);
    std::normal_distribution<double> d(0.0, sigma);
    r.smearFactor = 1.0 + d(rng);
  }
  return r;
}

// pT floor mirrors SmearedJetProducerT's MIN_JET_ENERGY clamp (adapted to pT
// since JetStruct carries no jet energy) -- avoids a negative/flipped jet.
inline double SmearedPt(double recoPt, double smearFactor,
                       double minPt = 1e-2) {
  double smeared = recoPt * smearFactor;
  return (smeared < minPt) ? minPt : smeared;
}

} // namespace JetSmearing

#endif
