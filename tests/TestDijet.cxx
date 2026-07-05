#include <iostream>
#include <cmath>

#include "TROOT.h"
#include "Dijet.h"

static int nPass = 0;
static int nFail = 0;
static constexpr float kEps = 1e-5f;

void Check(bool cond, const char *msg) {
  if (cond) {
    std::cout << "  PASS  " << msg << std::endl;
    nPass++;
  } else {
    std::cout << "  FAIL  " << msg << std::endl;
    nFail++;
  }
}

// ── FindLeadingJets ──────────────────────────────────────────────────────────

void TestFindLeadingJets() {
  std::cout << "\n[FindLeadingJets]\n";

  // 0 jets
  {
    SortedJets s = FindLeadingJets(nullptr, 0);
    Check(s.lead == -1 && s.sublead == -1 && s.third == -1,
          "0 jets: all indices -1");
  }

  // 1 jet
  {
    float pt[] = {50.0f};
    SortedJets s = FindLeadingJets(pt, 1);
    Check(s.lead == 0, "1 jet: lead = 0");
    Check(s.sublead == -1, "1 jet: sublead = -1");
    Check(s.third == -1, "1 jet: third = -1");
  }

  // 2 jets, already ordered
  {
    float pt[] = {100.0f, 80.0f};
    SortedJets s = FindLeadingJets(pt, 2);
    Check(s.lead == 0 && s.sublead == 1, "2 jets ordered: lead=0 sublead=1");
    Check(s.third == -1, "2 jets ordered: third=-1");
  }

  // 2 jets, reverse order
  {
    float pt[] = {60.0f, 90.0f};
    SortedJets s = FindLeadingJets(pt, 2);
    Check(s.lead == 1 && s.sublead == 0, "2 jets reversed: lead=1 sublead=0");
  }

  // 3 jets, mixed order: jet 2 is leading, jet 0 is sublead, jet 1 is third
  {
    float pt[] = {80.0f, 30.0f, 120.0f};
    SortedJets s = FindLeadingJets(pt, 3);
    Check(s.lead == 2, "3 jets mixed: lead=2 (pt=120)");
    Check(s.sublead == 0, "3 jets mixed: sublead=0 (pt=80)");
    Check(s.third == 1, "3 jets mixed: third=1 (pt=30)");
  }

  // 4 jets: only top 3 matter
  {
    float pt[] = {40.0f, 100.0f, 70.0f, 55.0f};
    SortedJets s = FindLeadingJets(pt, 4);
    Check(s.lead == 1, "4 jets: lead=1 (pt=100)");
    Check(s.sublead == 2, "4 jets: sublead=2 (pt=70)");
    Check(s.third == 3, "4 jets: third=3 (pt=55)");
  }
}

// ── MakeDijet ────────────────────────────────────────────────────────────────

// pi for a clean back-to-back dphi
static constexpr float kPi = 3.14159265f;

void TestMakeDijet() {
  std::cout << "\n[MakeDijet]\n";

  // base event: lead in barrel, sublead not
  // lead: pt=120, eta=0.3, phi=0
  // sublead: pt=100, eta=2.5, phi=pi (back-to-back)
  // third: pt=20, eta=1.0, phi=1.0
  {
    float pt[] = {120.0f, 100.0f, 20.0f};
    float eta[] = {0.3f, 2.5f, 1.0f};
    float phi[] = {0.0f, kPi, 1.0f};
    SortedJets s = {0, 1, 2};

    DijetResult d = MakeDijet(s, true, pt, eta, phi, 0ULL);
    Check(d.valid, "lead-only barrel: valid");
    Check(d.tagIdx == 0, "lead-only barrel: tagIdx = 0 (lead)");
    Check(d.probeIdx == 1, "lead-only barrel: probeIdx = 1 (sublead)");
  }

  // sublead in barrel, lead not
  {
    float pt[] = {120.0f, 100.0f};
    float eta[] = {2.5f, 0.3f};
    float phi[] = {0.0f, kPi};
    SortedJets s = {0, 1, -1};

    DijetResult d = MakeDijet(s, false, pt, eta, phi, 0ULL);
    Check(d.valid, "sublead-only barrel: valid");
    Check(d.tagIdx == 1, "sublead-only barrel: tagIdx = 1 (sublead)");
    Check(d.probeIdx == 0, "sublead-only barrel: probeIdx = 0 (lead)");
  }

  // neither in barrel: invalid
  {
    float pt[] = {120.0f, 100.0f};
    float eta[] = {2.0f, 2.5f};
    float phi[] = {0.0f, kPi};
    SortedJets s = {0, 1, -1};

    DijetResult d = MakeDijet(s, false, pt, eta, phi, 0ULL);
    Check(!d.valid, "neither in barrel: invalid");
  }

  // both in barrel, even event: lead is tag
  {
    float pt[] = {120.0f, 100.0f};
    float eta[] = {0.3f, -0.3f};
    float phi[] = {0.0f, kPi};
    SortedJets s = {0, 1, -1};

    DijetResult d = MakeDijet(s, false, pt, eta, phi, 4ULL);
    Check(d.valid, "both barrel even event: valid");
    Check(d.tagIdx == 0, "both barrel even event: tagIdx = 0 (lead)");
    Check(d.probeIdx == 1, "both barrel even event: probeIdx = 1 (sublead)");
  }

  // both in barrel, odd event: sublead is tag
  {
    float pt[] = {120.0f, 100.0f};
    float eta[] = {0.3f, -0.3f};
    float phi[] = {0.0f, kPi};
    SortedJets s = {0, 1, -1};

    DijetResult d = MakeDijet(s, false, pt, eta, phi, 7ULL);
    Check(d.valid, "both barrel odd event: valid");
    Check(d.tagIdx == 1, "both barrel odd event: tagIdx = 1 (sublead)");
    Check(d.probeIdx == 0, "both barrel odd event: probeIdx = 0 (lead)");
  }

  // dphi too small: invalid
  {
    float pt[] = {120.0f, 100.0f};
    float eta[] = {0.3f, 2.5f};
    float phi[] = {0.0f, 0.5f}; // dphi = 0.5, below 2.7
    SortedJets s = {0, 1, -1};

    DijetResult d = MakeDijet(s, false, pt, eta, phi, 0ULL);
    Check(!d.valid, "dphi < 2.7: invalid");
  }

  // sublead pT too low: invalid
  {
    float pt[] = {120.0f, 5.0f};
    float eta[] = {0.3f, 2.5f};
    float phi[] = {0.0f, kPi};
    SortedJets s = {0, 1, -1};

    DijetResult d = MakeDijet(s, false, pt, eta, phi, 0ULL);
    Check(!d.valid, "sublead pt < 10 GeV: invalid");
  }

  // ptavg and A values
  // lead (tag): pt=120, sublead (probe): pt=80
  // ptavg = 100, A = (80-120)/200 = -0.2
  {
    float pt[] = {120.0f, 80.0f};
    float eta[] = {0.3f, 2.5f};
    float phi[] = {0.0f, kPi};
    SortedJets s = {0, 1, -1};

    DijetResult d = MakeDijet(s, false, pt, eta, phi, 0ULL);
    Check(d.valid, "ptavg/A: valid");
    Check(std::fabs(d.ptavg - 100.0f) < kEps, "ptavg = 100");
    Check(std::fabs(d.A - (-0.2f)) < kEps, "A = -0.2 (tag is higher pT)");
  }

  // A is positive when probe is higher pT (sublead is tag via odd event)
  {
    float pt[] = {120.0f, 80.0f};
    float eta[] = {0.3f, -0.3f};
    float phi[] = {0.0f, kPi};
    SortedJets s = {0, 1, -1};

    // odd event → sublead (pt=80) is tag, lead (pt=120) is probe
    DijetResult d = MakeDijet(s, false, pt, eta, phi, 1ULL);
    Check(std::fabs(d.A - 0.2f) < kEps, "A = +0.2 (probe is higher pT)");
  }

  // alpha = 0 when no third jet
  {
    float pt[] = {120.0f, 80.0f};
    float eta[] = {0.3f, 2.5f};
    float phi[] = {0.0f, kPi};
    SortedJets s = {0, 1, -1};

    DijetResult d = MakeDijet(s, false, pt, eta, phi, 0ULL);
    Check(std::fabs(d.alpha) < kEps, "alpha = 0 with no third jet");
  }

  // alpha = third.pt / ptavg
  // lead pt=120, sublead pt=80 → ptavg=100, third pt=30 → alpha=0.3
  {
    float pt[] = {120.0f, 80.0f, 30.0f};
    float eta[] = {0.3f, 2.5f, 1.0f};
    float phi[] = {0.0f, kPi, 0.5f};
    SortedJets s = {0, 1, 2};

    DijetResult d = MakeDijet(s, true, pt, eta, phi, 0ULL);
    Check(std::fabs(d.alpha - 0.3f) < kEps,
          "alpha = 0.3 (third.pt=30 / ptavg=100)");
  }

  // hasThird = false suppresses alpha even if third index is set
  {
    float pt[] = {120.0f, 80.0f, 30.0f};
    float eta[] = {0.3f, 2.5f, 1.0f};
    float phi[] = {0.0f, kPi, 0.5f};
    SortedJets s = {0, 1, 2};

    DijetResult d = MakeDijet(s, false, pt, eta, phi, 0ULL);
    Check(std::fabs(d.alpha) < kEps,
          "hasThird=false: alpha = 0 regardless of third index");
  }

  // dphi stored in result
  {
    float pt[] = {120.0f, 80.0f};
    float eta[] = {0.3f, 2.5f};
    float phi[] = {0.0f, kPi};
    SortedJets s = {0, 1, -1};

    DijetResult d = MakeDijet(s, false, pt, eta, phi, 0ULL);
    Check(std::fabs(d.dphi - kPi) < kEps, "dphi stored in result");
  }

  // dphi wraps correctly: jets near ±pi straddle the boundary
  // phi1 = pi-0.1, phi2 = -(pi-0.1) → raw diff = 2pi-0.2 > pi → wraps to 0.2
  {
    float phi1 = kPi - 0.1f;
    float phi2 = -(kPi - 0.1f);
    float got = DijetDPhi(phi1, phi2);
    Check(std::fabs(got - 0.2f) < kEps,
          "DijetDPhi wraps: jets near ±pi → dphi = 0.2");
  }
}

int main() {
  gROOT->SetBatch(true);

  std::cout << "=== TestDijet ===" << std::endl;

  TestFindLeadingJets();
  TestMakeDijet();

  std::cout << "\n=== " << nPass << " passed, " << nFail
            << " failed ===" << std::endl;
  return nFail > 0 ? 1 : 0;
}
