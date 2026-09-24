#ifndef DIJETPROFILES_H
#define DIJETPROFILES_H

#include "TDirectory.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile2D.h"
#include "TString.h"

#include <cmath>
#include <vector>

// Input histograms for JME's dijet L2Res/JER framework (l2l3res-multijet,
// DijetHistosFill.C "Dijet2" folder), filled from HiForest jets instead of
// NanoAOD. Same object names, binning and HDM bisector definitions
// (JME-21-001), so histogram_scripts/DijetHistosJER.C and
// src/DijetHistosL2Res.C read the output directly (dir = "<cone>/Dijet2").
//
// Jets only: DB (p2m2*) and the jets-only FSR term (p2mn*) are filled. The
// MET-based MPF objects (p2m0*, p2mu*, p2mnu*, h3m0) are booked empty so
// the JME macros' asserts pass -- no MET in the forests yet.

namespace dijet2 {

// JME "isdijet2" (DESY) selection, DijetHistosFill.C
static constexpr double kTagEtaMax = 1.3;
static constexpr double kDphiMin = 2.7;
static constexpr double kAsymMax = 0.7; // |pT,tag - pT,probe| / sum
static constexpr double kJetPtMin = 15.0;
// HIN deviation: jets beyond this |eta| are exempt from allJetsGood and are
// never probes. HiForest's PF multiplicities exclude HF, so every |eta| > 3
// jet fails the ID -- JME's check would otherwise reject every event with a
// > 15 GeV HF jet. The L2Res/JER here target |eta| < 3.
static constexpr double kMaxIdEta = 3.0;

// L2Res |eta| binning (vxd) and pT binning (vptd), DijetHistosFill.C
inline const std::vector<double> &AbsEtaEdges() {
  static const std::vector<double> v = {
      0,     0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696,
      0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
      1.566, 1.653, 1.74,  1.83,  1.93,  2.043, 2.172, 2.322, 2.5,
      2.65,  2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
      4.191, 4.363, 4.538, 4.716, 4.889, 5.191};
  return v;
}

inline const std::vector<double> &PtEdges() {
  static const std::vector<double> v = {
      10,   15,   21,   28,   37,   49,   59,   86,   110,  132,  170,  204,
      236,  279,  302,  373,  460,  575,  638,  737,  846,  967,  1101, 1248,
      1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103, 3450, 4037, 5220};
  return v;
}

// per (tag, probe) ordering, transverse plane only
struct Balance {
  double absEtaProbe = 0;
  double ptavp = 0; // bisector-projected pT average
  double ptTag = 0;
  double ptProbe = 0;
  double m2b = 0, m2bx = 0, mnb = 0, mnbx = 0; // bisector axis (+ 90 deg)
  double m2c = 0, mnc = 0;                     // tag axis
  double m2f = 0, mnf = 0;                     // probe axis
};

// mnX/mnY: -(sum of every other jet above kJetPtMin)
inline Balance ComputeBalance(double ptT, double phiT, double ptP, double phiP,
                              double etaP, double mnX, double mnY) {
  Balance b;
  const double tx = std::cos(phiT), ty = std::sin(phiT);
  const double px = std::cos(phiP), py = std::sin(phiP);

  // bisector: unit(t - p), and rotated by +90 deg
  double bx = tx - px, by = ty - py;
  const double bn = std::hypot(bx, by);
  bx /= bn;
  by /= bn;
  const double bxx = -by, bxy = bx;

  // -(p_tag + p_probe)
  const double m2X = -(ptT * tx + ptP * px);
  const double m2Y = -(ptT * ty + ptP * py);

  b.absEtaProbe = std::fabs(etaP);
  b.ptTag = ptT;
  b.ptProbe = ptP;
  b.ptavp = 0.5 * (ptT * (tx * bx + ty * by) - ptP * (px * bx + py * by));

  b.m2b = 1 + (m2X * bx + m2Y * by) / b.ptavp;
  b.m2bx = 1 + (m2X * bxx + m2Y * bxy) / b.ptavp;
  b.mnb = (mnX * bx + mnY * by) / b.ptavp;
  b.mnbx = (mnX * bxx + mnY * bxy) / b.ptavp;

  b.m2c = 1 + (m2X * tx + m2Y * ty) / ptT;
  b.mnc = (mnX * tx + mnY * ty) / ptT;

  b.m2f = 1 + (m2X * -px + m2Y * -py) / ptP;
  b.mnf = (mnX * -px + mnY * -py) / ptP;
  return b;
}

} // namespace dijet2

struct DijetProfiles {

  // bisector (pT,avp) binning
  TH2D *h2pteta = nullptr;
  TProfile2D *p2m2 = nullptr, *p2m2x = nullptr;
  TProfile2D *p2mn = nullptr, *p2mnx = nullptr;
  TH3D *h3m2 = nullptr;
  // tag (pT,tag) and probe (pT,probe) binning
  TH2D *h2ptetatc = nullptr, *h2ptetapf = nullptr;
  TProfile2D *p2m2tc = nullptr, *p2mntc = nullptr;
  TProfile2D *p2m2pf = nullptr, *p2mnpf = nullptr;
  // MET-based, empty
  std::vector<TObject *> placeholders;

  void Init() {
    TDirectory::TContext detached(nullptr);
    const auto &ve = dijet2::AbsEtaEdges();
    const auto &vp = dijet2::PtEdges();
    const int ne = (int)ve.size() - 1, np = (int)vp.size() - 1;

    auto p2 = [&](const char *name, const char *axes, const char *z,
                  const char *opt = "") {
      return new TProfile2D(name, Form(";|#eta|;%s (GeV);%s", axes, z), ne,
                            ve.data(), np, vp.data(), opt);
    };
    auto h2 = [&](const char *name, const char *axes) {
      TH2D *h = new TH2D(name, Form(";|#eta|;%s (GeV);N_{events}", axes), ne,
                         ve.data(), np, vp.data());
      h->Sumw2();
      return h;
    };
    std::vector<double> vm(201);
    for (int i = 0; i <= 200; i++) {
      vm[i] = 0.01 * i;
    }
    auto h3 = [&](const char *name, const char *z) {
      TH3D *h = new TH3D(name, Form(";|#eta|;p_{T,avp} (GeV);%s", z), ne,
                         ve.data(), np, vp.data(), 200, vm.data());
      h->Sumw2();
      return h;
    };

    const char *avp = "p_{T,avp}";
    h2pteta = h2("h2pteta", avp);
    p2m2 = p2("p2m2", avp, "MPF2");
    p2m2x = p2("p2m2x", avp, "MPF2X (DBX)", "S");
    p2mn = p2("p2mn", avp, "MPFn");
    p2mnx = p2("p2mnx", avp, "MPFNX", "S");
    h3m2 = h3("h3m2", "MPF2");

    h2ptetatc = h2("h2ptetatc", "p_{T,tag}");
    p2m2tc = p2("p2m2tc", "p_{T,tag}", "MPF2");
    p2mntc = p2("p2mntc", "p_{T,tag}", "MPFn");
    h2ptetapf = h2("h2ptetapf", "p_{T,probe}");
    p2m2pf = p2("p2m2pf", "p_{T,probe}", "MPF2");
    p2mnpf = p2("p2mnpf", "p_{T,probe}", "MPFn");

    placeholders = {p2("p2m0", avp, "MPF0"),
                    p2("p2m0x", avp, "MPF0X (MPFX)", "S"),
                    p2("p2mu", avp, "MPFu"),
                    p2("p2mux", avp, "MPFUX", "S"),
                    p2("p2mnu", avp, "MPFnu"),
                    p2("p2mnux", avp, "MPFNUX", "S"),
                    h3("h3m0", "MPF0"),
                    p2("p2m0tc", "p_{T,tag}", "MPF0"),
                    p2("p2mutc", "p_{T,tag}", "MPFu"),
                    p2("p2m0pf", "p_{T,probe}", "MPF0"),
                    p2("p2mupf", "p_{T,probe}", "MPFu")};
  }

  // useX: whether this dataset owns the bin of that pT variable (trigger
  // split between the triggered and non-triggered samples)
  void Fill(const dijet2::Balance &b, double w, bool useAvp, bool useTag,
            bool useProbe) {
    const double eta = b.absEtaProbe;
    if (useAvp) {
      h2pteta->Fill(eta, b.ptavp, w);
      p2m2->Fill(eta, b.ptavp, b.m2b, w);
      p2m2x->Fill(eta, b.ptavp, b.m2bx, w);
      p2mn->Fill(eta, b.ptavp, b.mnb, w);
      p2mnx->Fill(eta, b.ptavp, b.mnbx, w);
      h3m2->Fill(eta, b.ptavp, b.m2b, w);
    }
    if (useTag) {
      h2ptetatc->Fill(eta, b.ptTag, w);
      p2m2tc->Fill(eta, b.ptTag, b.m2c, w);
      p2mntc->Fill(eta, b.ptTag, b.mnc, w);
    }
    if (useProbe) {
      h2ptetapf->Fill(eta, b.ptProbe, w);
      p2m2pf->Fill(eta, b.ptProbe, b.m2f, w);
      p2mnpf->Fill(eta, b.ptProbe, b.mnf, w);
    }
  }

  void Write(TDirectory *dir) const {
    TDirectory::TContext ctx(dir);
    for (TObject *o : std::vector<TObject *>{h2pteta, p2m2, p2m2x, p2mn, p2mnx,
                                             h3m2, h2ptetatc, p2m2tc, p2mntc,
                                             h2ptetapf, p2m2pf, p2mnpf}) {
      o->Write();
    }
    for (TObject *o : placeholders) {
      o->Write();
    }
  }
};

#endif
