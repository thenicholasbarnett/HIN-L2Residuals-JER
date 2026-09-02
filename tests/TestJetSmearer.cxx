#include <cmath>
#include <iostream>

#include "TROOT.h"

#include "JetSmearer.h"

static int nPass = 0;
static int nFail = 0;
static constexpr double kEps = 1e-6;

void Check(bool cond, const char *msg) {
  if (cond) {
    std::cout << "  PASS  " << msg << std::endl;
    nPass++;
  } else {
    std::cout << "  FAIL  " << msg << std::endl;
    nFail++;
  }
}

int main() {
  gROOT->SetBatch(true);

  std::cout << "=== TestJetSmearer ===" << std::endl;

  const std::string kResFile =
      "data/jec/preliminary/2024ppref_V1_MC_PtResolution_ak4PFchs.txt";
  const std::string kSfFile =
      "data/jec/preliminary/2024ppref_V1_MC_SF_ak4PFchs.txt";

  std::cout << "\n[1] Smearer construction from real JER text files"
            << std::endl;
  {
    bool threw = false;
    try {
      JetSmearer smearer(kResFile, kSfFile);
      (void)smearer;
    } catch (...) {
      threw = true;
    }
    Check(!threw, "Smearer(resFile, sfFile) constructs without throwing");
  }

  const double recoPt = 100.0, eta = 0.3, rho = 1.5;

  std::cout << "\n[2] Gen-matched jet uses the deterministic scaling branch"
            << std::endl;
  {
    JetSmearer smearer(kResFile, kSfFile);
    const double genPt = 98.0; // close enough to recoPt to count as matched
    JetSmearing::Result r = smearer.Smear(recoPt, eta, rho, genPt);
    Check(r.matched, "matched == true for a close gen pT");
    Check(r.resolution > 0.0, "resolution > 0");

    // deterministic scaling formula from JetSmearer.h's JetSmearing::
    // namespace -- recompute
    // independently against the same resolution/SF objects and compare
    JetSmearerJME::JetResolution resolution(kResFile);
    JetSmearerJME::JetResolutionScaleFactor scaleFactor(kSfFile);
    double expectSf = scaleFactor.getScaleFactor(
        JetSmearerJME::JetParameters().setJetPt(recoPt).setJetEta(eta));
    double expectFactor =
        1.0 + (expectSf - 1.0) * (recoPt - genPt) / recoPt;
    Check(std::abs(r.smearFactor - expectFactor) < kEps,
          "smearFactor matches the hand-computed scaling formula");

    // deterministic branch draws no random number, so repeated calls (and
    // the SmearedPt() convenience wrapper) must reproduce it exactly
    JetSmearing::Result r2 = smearer.Smear(recoPt, eta, rho, genPt);
    Check(std::abs(r2.smearFactor - r.smearFactor) < kEps,
          "repeated Smear() calls agree on the deterministic branch");
    double smearedPt = smearer.SmearedPt(recoPt, eta, rho, genPt);
    Check(std::abs(smearedPt - recoPt * r.smearFactor) < kEps,
          "SmearedPt() == recoPt * smearFactor");
  }

  std::cout << "\n[3] Unmatched jet (genPt < 0) never uses the gen-matched "
               "branch"
            << std::endl;
  {
    JetSmearer smearer(kResFile, kSfFile);
    JetSmearing::Result r = smearer.Smear(recoPt, eta, rho, -1.0);
    Check(!r.matched, "matched == false when genPt < 0");
    Check(std::isfinite(r.smearFactor), "smearFactor is finite");
    Check(r.smearFactor > 0.0, "smearFactor > 0");
  }

  std::cout << "\n[4] SmearedPt respects the pT floor (JetSmearing::SmearedPt)"
            << std::endl;
  { Check(JetSmearing::SmearedPt(10.0, -5.0) > 0.0,
         "a wildly negative smear factor still returns a positive pT"); }

  std::cout << "\n=== " << nPass << " passed, " << nFail
            << " failed ===" << std::endl;
  return nFail > 0 ? 1 : 0;
}
