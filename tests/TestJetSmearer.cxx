#include <cmath>
#include <cstdint>
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
  const std::int64_t evt = 123456789;

  std::cout << "\n[2] Gen-matched jet uses the deterministic scaling branch"
            << std::endl;
  {
    JetSmearer smearer(kResFile, kSfFile);
    const double genPt = 98.0; // close enough to recoPt to count as matched
    JetSmearing::Result r = smearer.Smear(recoPt, eta, rho, genPt, evt);
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
    JetSmearing::Result r2 = smearer.Smear(recoPt, eta, rho, genPt, evt);
    Check(std::abs(r2.smearFactor - r.smearFactor) < kEps,
          "repeated Smear() calls agree on the deterministic branch");
    double smearedPt = smearer.SmearedPt(recoPt, eta, rho, genPt, evt);
    Check(std::abs(smearedPt - recoPt * r.smearFactor) < kEps,
          "SmearedPt() == recoPt * smearFactor");
  }

  std::cout << "\n[3] Unmatched jet (genPt < 0) never uses the gen-matched "
               "branch"
            << std::endl;
  {
    JetSmearer smearer(kResFile, kSfFile);
    JetSmearing::Result r = smearer.Smear(recoPt, eta, rho, -1.0, evt);
    Check(!r.matched, "matched == false when genPt < 0");
    Check(std::isfinite(r.smearFactor), "smearFactor is finite");
    Check(r.smearFactor > 0.0, "smearFactor > 0");
    // seeded from (pT, eta, rho, event), no hidden state
    Check(smearer.Smear(recoPt, eta, rho, -1.0, evt).smearFactor ==
              r.smearFactor,
          "same jet, same event: same stochastic smear");
    Check(smearer.Smear(recoPt, eta, rho, -1.0, evt + 1).smearFactor !=
              r.smearFactor,
          "next event: different stochastic smear");
  }

  std::cout << "\n[4] SmearedPt has no pT floor (as JERSmear)" << std::endl;
  {
    Check(JetSmearing::SmearedPt(10.0, -5.0) == -50.0,
          "SmearedPt(10, -5) == -50");
  }

  std::cout << "\n[5] JERSmear bit-identical to correctionlib 2.9.0 "
               "(jer_smear.json)"
            << std::endl;
  {
    // pt, eta, genPt, rho, event, jer, sf -> correctionlib output: 4
    // stochastic, 1 with SF < 1 (no smear), 3 scaling, genPt = 0
    struct Row {
      double pt, eta, gen, rho;
      long long eid;
      double jer, sf, out;
    };
    const Row rows[] = {
        {1901.175074170241, 1.2003216515760828, -1.7605236605147883,
         49.996041183624136, 498466514399LL, 0.27191330920905893,
         1.2697298083646846, 0.9158300351062498},
        {292.5984273756693, -2.9957116995516833, -1.3369702368047562,
         5.3796858018933875, 251184097034LL, 0.0813456794597646,
         1.5732070721663405, 0.9781077384817323},
        {1897.5556470388015, 2.2178430356851724, -0.8897957898451898,
         39.72878710804939, 73745862620LL, 0.31865445043870977,
         1.5394817815316748, 1.283631533141122},
        {627.1037467609185, -4.06746032951771, -1.968805762920349,
         24.383448444178285, 787931576372LL, 0.1514453582191491,
         1.1991746744491527, 0.9643934955832941},
        {1026.0841412770121, 0.5590372795786862, -2.4167598158391037,
         42.949664325946955, 67542619370LL, 0.1272446636891928,
         0.9017280067451869, 1.0},
        {219.7850830459586, -4.060994193998628, 150.43631497551985,
         57.052220298905816, 12914606801LL, 0.20807274131368877,
         1.4548616135213392, 1.1435224452136308},
        {329.18086494644473, 3.0548883742955906, 173.54084198906205,
         29.278166134811695, 839826462277LL, 0.17091251445441266,
         1.0084658294555655, 1.0040027292930038},
        {1393.7750427427804, 1.52252903549085, 1310.947289856692,
         29.3045044750378, 786162560040LL, 0.19510196312948058,
         1.1800754725466016, 1.0107013300450438},
        {1086.9413704608742, -1.9940269312924461, 0.0, 22.428722902239823,
         246936153312LL, 0.05889133843847026, 1.0070368556329086,
         1.0070368556329086},
    };
    int exact = 0;
    for (const Row &r : rows) {
      exact += JetSmearing::JERSmear(r.pt, r.eta, r.gen, r.rho, r.eid, r.jer,
                                     r.sf) == r.out;
    }
    Check(exact == (int)(sizeof(rows) / sizeof(rows[0])), "all 9 rows exact");
  }

  std::cout << "\n=== " << nPass << " passed, " << nFail
            << " failed ===" << std::endl;
  return nFail > 0 ? 1 : 0;
}
