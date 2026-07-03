// Vendored verbatim from JetMETAnalysis/JetUtilities/interface/ProgressBar.hh
// (https://github.com/lmartika/JetMETAnalysis) into L2Residuals-2024ppref on
// 2026-07-03. Internal-only: used by exactly one call site inside
// HistogramUtilities.cc (a debug-gated loadbar2() in the median calculation).
// Not used anywhere else in this repo -- include/ProgressBar.h (richer: rate,
// ETA, color) remains the sole progress-bar system everywhere else, per the
// original JME compatibility survey's assessment that this file is a strict
// subset of what this repo already has.

#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>

using std::cout;
using std::endl;
using std::flush;
using std::setw;
using std::string;
using std::remainder;

const int r = 100;

inline void loadbar2(unsigned int x, unsigned int n, unsigned int w = 50, string prefix = "") {
   if ( (x != n) && ((n >= 100) && x % (n/100) != 0) ) return;

   float ratio  =  x/(float)n;
   int   c      =  ratio * w;

   cout << prefix << setw(3) << (int)(ratio*100) << "% [";
   for (int x=0; x<c; x++) cout << "=";
   for (unsigned int x=c; x<w; x++) cout << " ";
   cout << "] (" << x << "/" << n << ")\r" << flush;
}

inline void loadbar3(unsigned int x, unsigned int n, unsigned int w = 50, string prefix = "") {
   if ( (x != n) && ((n >= 10000) && x % (n/10000) != 0) ) return;
   //// round to nearest r
   //unsigned int ntmp = n/10000;
   //ntmp += r/2; // 1: Add r/2.
   //ntmp /= r;   // 2: Divide by r.
   //ntmp *= r;   // 3: Multiply by r
   //if ( (x != n) && (remainder(x,ntmp) != 0) ) return;

   float ratio  =  x/(float)n;
   int   c      =  ratio * w;

   cout << prefix << std::fixed << setw(8) << std::setprecision(2) << (ratio*100) << "% [";
   for (int x=0; x<c; x++) cout << "=";
   for (unsigned int x=c; x<w; x++) cout << " ";
   cout << "] (" << x << "/" << n << ")\r" << flush;
}
