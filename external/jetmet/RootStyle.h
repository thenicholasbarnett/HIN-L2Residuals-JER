// Vendored verbatim from JetMETAnalysis/JetUtilities/interface/RootStyle.h
// (https://github.com/lmartika/JetMETAnalysis) into L2Residuals-2024ppref on
// 2026-07-03. Dependency-free (just TROOT.h/TStyle.h). setTDRStyle() is the
// one used by external/jetmet/Style.h's tdrCanvas/tdrDiCanvas and by this
// repo's own pilot integration (include/plotting/EventQA.h,
// FinalCorrections.h) -- see Style.h for the CMS_lumi()/tdrDraw()/tdrLeg()
// helpers that go with it.

#ifndef JETUTILITIES_ROOTSTYLE_H
#define JETUTILITIES_ROOTSTYLE_H 1


void set_root_style();
void set_tdr_style();
void set_vasu_style();
void setTDRStyle();


#endif
