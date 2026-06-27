#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"

#include "Binning.h"
#include "TextFileWriter.h"

static int nPass = 0;
static int nFail = 0;
static constexpr double kEps = 1e-9;

void Check(bool cond, const char* msg) {
    if (cond) { std::cout << "  PASS  " << msg << std::endl; nPass++; }
    else       { std::cout << "  FAIL  " << msg << std::endl; nFail++; }
}

// ---- helpers ----

// Build a synthetic residuals ROOT file with `nSlices` pT-slice intercept
// TH1Ds for ak4PF/gauss, each bin set to `corrValue` ± `corrErr`.
// Use corrValue=0, corrErr=0 to produce entirely empty histograms.

static TString MakeResiduals(const char* path, double corrValue, double corrErr,
                              int nSlices = 5) {
    TFile* f = new TFile(path, "recreate");
    BinningConfig bins;
    const int nEta = (int)kAbsEtaEdges.size() - 1;
    for (int ip = 0; ip < nSlices; ip++) {
        TString name = Form("ak4PF_intercept_gauss%s",
            bins.ptavgSlices[ip].shortName.Data());
        TH1D* h = new TH1D(name, "", nEta, kAbsEtaEdges.data());
        for (int ieta = 1; ieta <= nEta; ieta++) {
            h->SetBinContent(ieta, corrValue);
            h->SetBinError(ieta, corrErr);
        }
        h->Write();
        delete h;
    }
    f->Close();
    delete f;
    return TString(path);
}

// Parse one data line from the JEC text file.
// Returns false if parsing fails or field count is wrong.
struct JECLine {
    double etaLo, etaHi;
    int    npar;
    double ptLo, ptHi;
    double p[3];
};

static bool ParseJECLine(const std::string& line, JECLine& out) {
    std::istringstream ss(line);
    if (!(ss >> out.etaLo >> out.etaHi >> out.npar >> out.ptLo >> out.ptHi)) return false;
    for (int i = 0; i < 3; i++) if (!(ss >> out.p[i])) return false;
    return true;
}

// Read all JEC data lines (skip the header) from a text file.
static std::vector<JECLine> ReadJECFile(const char* path, std::string& header) {
    std::ifstream f(path);
    std::vector<JECLine> lines;
    bool firstLine = true;
    std::string line;
    while (std::getline(f, line)) {
        if (firstLine) { header = line; firstLine = false; continue; }
        if (line.empty()) continue;
        JECLine jl;
        if (ParseJECLine(line, jl)) lines.push_back(jl);
    }
    return lines;
}

// ---- test cases ----

// [1] Header line and overall structure
// Input: 5 pT slices with correction=1.0 (unity everywhere — fit should converge)
void TestFileStructure() {
    std::cout << "\n[1] File structure\n";

    const char* rootPath = "/tmp/tw_test1.root";
    const char* txtPath  = "/tmp/tw_test1.txt";
    std::remove(txtPath);

    MakeResiduals(rootPath, 1.0, 0.001, 5);
    runTextFile(rootPath, txtPath, "gauss", "ak4PF");

    std::string header;
    auto lines = ReadJECFile(txtPath, header);

    Check(header.rfind("{1 JetEta", 0) == 0,  "header starts with {1 JetEta");
    Check(header.find("L2Residual") != std::string::npos, "header contains L2Residual");
    Check(header.find("log10")      != std::string::npos, "header contains log10 formula");
    Check((int)lines.size() == 36,             "36 data lines total");

    if (!lines.empty()) {
        Check(lines.front().npar == 5,  "Npar = 5");
        Check(std::fabs(lines.front().ptLo - 40.0)   < kEps, "pT_lo = 40");
        Check(std::fabs(lines.front().ptHi - 1000.0) < kEps, "pT_hi = 1000");
    }

    std::remove(rootPath);
    std::remove(txtPath);
}

// [2] Eta ordering: negative outer→inner, then positive inner→outer
void TestEtaOrdering() {
    std::cout << "\n[2] Eta ordering\n";

    const char* rootPath = "/tmp/tw_test2.root";
    const char* txtPath  = "/tmp/tw_test2.txt";
    std::remove(txtPath);

    MakeResiduals(rootPath, 1.0, 0.001, 5);
    runTextFile(rootPath, txtPath, "gauss", "ak4PF");

    std::string header;
    auto lines = ReadJECFile(txtPath, header);
    if ((int)lines.size() != 36) {
        std::cout << "  SKIP  (line count wrong)\n";
        std::remove(rootPath); std::remove(txtPath);
        return;
    }

    // First line: most negative eta (outermost)
    Check(lines[0].etaLo < -5.0,   "first line etaLo < -5");
    Check(lines[0].etaHi < 0.0,    "first line etaHi < 0");

    // Line 18 (index 17): innermost negative eta
    Check(lines[17].etaLo < 0.0,   "line 18 etaLo < 0 (inner negative)");
    Check(std::fabs(lines[17].etaHi) < kEps, "line 18 etaHi = 0 (boundary)");

    // Line 19 (index 18): innermost positive eta
    Check(std::fabs(lines[18].etaLo) < kEps, "line 19 etaLo = 0 (boundary)");
    Check(lines[18].etaHi > 0.0,   "line 19 etaHi > 0 (inner positive)");

    // Last line: outermost positive eta
    Check(lines[35].etaLo > 0.0,   "last line etaLo > 0");
    Check(lines[35].etaHi > 5.0,   "last line etaHi > 5");

    // All negative-eta lines have etaLo < etaHi < 0 (or etaHi == 0 for innermost)
    bool negOK = true;
    for (int i = 0; i < 18; i++)
        if (!(lines[i].etaLo < lines[i].etaHi && lines[i].etaHi <= 0.0)) negOK = false;
    Check(negOK, "all negative-eta lines have etaLo < etaHi ≤ 0");

    // All positive-eta lines have 0 ≤ etaLo < etaHi
    bool posOK = true;
    for (int i = 18; i < 36; i++)
        if (!(lines[i].etaLo >= 0.0 && lines[i].etaLo < lines[i].etaHi)) posOK = false;
    Check(posOK, "all positive-eta lines have 0 ≤ etaLo < etaHi");

    // Negative and positive halves are mirror images (|eta| bin boundaries match)
    // line[17] = innermost negative, line[18] = innermost positive
    // line[0] = outermost negative, line[35] = outermost positive
    bool mirrorOK = true;
    for (int i = 0; i < 18; i++) {
        double absLoNeg = -lines[17 - i].etaHi;  // |eta_lo| of negative bin i from center
        double absHiNeg = -lines[17 - i].etaLo;
        double absloPOS = lines[18 + i].etaLo;
        double abshiPOS = lines[18 + i].etaHi;
        if (std::fabs(absLoNeg - absloPOS) > kEps || std::fabs(absHiNeg - abshiPOS) > kEps)
            mirrorOK = false;
    }
    Check(mirrorOK, "negative and positive eta bin edges are symmetric");

    std::remove(rootPath);
    std::remove(txtPath);
}

// [3] Symmetric corrections: negative and positive halves get the same fit params
void TestSymmetricCorrections() {
    std::cout << "\n[3] Symmetric corrections\n";

    const char* rootPath = "/tmp/tw_test3.root";
    const char* txtPath  = "/tmp/tw_test3.txt";
    std::remove(txtPath);

    MakeResiduals(rootPath, 1.0, 0.001, 5);
    runTextFile(rootPath, txtPath, "gauss", "ak4PF");

    std::string header;
    auto lines = ReadJECFile(txtPath, header);
    if ((int)lines.size() != 36) {
        std::cout << "  SKIP  (line count wrong)\n";
        std::remove(rootPath); std::remove(txtPath);
        return;
    }

    // Negative bin i (from center) = lines[17-i], positive bin i = lines[18+i]
    // They correspond to the same |eta| correction → parameters must be identical.
    bool paramsMatch = true;
    for (int i = 0; i < 18; i++) {
        const JECLine& neg = lines[17 - i];
        const JECLine& pos = lines[18 + i];
        for (int p = 0; p < 3; p++) {
            if (std::fabs(neg.p[p] - pos.p[p]) > 1e-12) { paramsMatch = false; break; }
        }
    }
    Check(paramsMatch, "negative and positive sides have identical fit parameters");

    std::remove(rootPath);
    std::remove(txtPath);
}

// [4] Unity fallback when fewer than kMinSlices (3) pT slices are available
void TestUnityFallback_TooFewSlices() {
    std::cout << "\n[4] Unity fallback — fewer than kMinSlices pT slices\n";

    const char* rootPath = "/tmp/tw_test4.root";
    const char* txtPath  = "/tmp/tw_test4.txt";
    std::remove(txtPath);

    // Only 2 pT slices — not enough to fit 3 parameters
    MakeResiduals(rootPath, 1.0, 0.001, 2);
    runTextFile(rootPath, txtPath, "gauss", "ak4PF");

    std::string header;
    auto lines = ReadJECFile(txtPath, header);
    if ((int)lines.size() != 36) {
        std::cout << "  SKIP  (line count wrong)\n";
        std::remove(rootPath); std::remove(txtPath);
        return;
    }

    bool allUnity = true;
    for (const auto& l : lines) {
        if (std::fabs(l.p[0] - 1.0) > kEps ||
            std::fabs(l.p[1] - 0.0) > kEps ||
            std::fabs(l.p[2] - 0.0) > kEps) {
            allUnity = false;
            break;
        }
    }
    Check(allUnity, "all bins fall back to unity (p0=1, p1=0, p2=0)");

    std::remove(rootPath);
    std::remove(txtPath);
}

// [5] Unity fallback for empty histograms (all bins zero)
void TestUnityFallback_EmptyBins() {
    std::cout << "\n[5] Unity fallback — empty histograms\n";

    const char* rootPath = "/tmp/tw_test5.root";
    const char* txtPath  = "/tmp/tw_test5.txt";
    std::remove(txtPath);

    // corrValue=0, corrErr=0 → all bins empty → no valid pT slices → unity
    MakeResiduals(rootPath, 0.0, 0.0, 5);
    runTextFile(rootPath, txtPath, "gauss", "ak4PF");

    std::string header;
    auto lines = ReadJECFile(txtPath, header);
    if ((int)lines.size() != 36) {
        std::cout << "  SKIP  (line count wrong)\n";
        std::remove(rootPath); std::remove(txtPath);
        return;
    }

    bool allUnity = true;
    for (const auto& l : lines) {
        if (std::fabs(l.p[0] - 1.0) > kEps ||
            std::fabs(l.p[1] - 0.0) > kEps ||
            std::fabs(l.p[2] - 0.0) > kEps) {
            allUnity = false;
            break;
        }
    }
    Check(allUnity, "empty bins produce unity correction (p0=1, p1=0, p2=0)");

    std::remove(rootPath);
    std::remove(txtPath);
}

// [6] Fit round-trip: set corrections from known params, check fit recovers them
// Uses params (p0=1.0, p1=0.0, p2=0.0) → constant correction=1.0 everywhere.
// The fit should converge to a solution very close to (1,0,0).
void TestFitRoundTrip() {
    std::cout << "\n[6] Fit round-trip\n";

    const char* rootPath = "/tmp/tw_test6.root";
    const char* txtPath  = "/tmp/tw_test6.txt";
    std::remove(txtPath);

    // Evaluate 1/(p0+p1*log10(0.01*pT)+p2*10/pT) at each pT center.
    // With p0=1, p1=0, p2=0 → correction=1.0 everywhere.
    // The fit should recover params that give f≈1.0 at all pT centers.
    MakeResiduals(rootPath, 1.0, 1e-5, 5);
    runTextFile(rootPath, txtPath, "gauss", "ak4PF");

    std::string header;
    auto lines = ReadJECFile(txtPath, header);
    if ((int)lines.size() != 36) {
        std::cout << "  SKIP  (line count wrong)\n";
        std::remove(rootPath); std::remove(txtPath);
        return;
    }

    // Evaluate the fit at each of the 5 pT centers for the innermost barrel bin.
    // The output correction should be close to 1.0 if the fit converged.
    const double ptCenters[] = {65.0, 105.0, 155.0, 225.0, 630.0};
    const JECLine& barrel = lines[18];  // innermost positive eta
    bool fitsClose = true;
    for (double pT : ptCenters) {
        double denom = barrel.p[0]
                     + barrel.p[1] * std::log10(0.01 * pT)
                     + barrel.p[2] / (pT / 10.0);
        double corr = (denom != 0.0) ? 1.0 / denom : 0.0;
        if (std::fabs(corr - 1.0) > 0.01) { fitsClose = false; break; }
    }
    Check(fitsClose, "fit recovers f≈1.0 at all pT centers for unit input");

    std::remove(rootPath);
    std::remove(txtPath);
}

// [7] CMS JEC edge value: outermost eta bin ends at 5.191
void TestEtaExtent() {
    std::cout << "\n[7] CMS eta extent\n";

    const char* rootPath = "/tmp/tw_test7.root";
    const char* txtPath  = "/tmp/tw_test7.txt";
    std::remove(txtPath);

    MakeResiduals(rootPath, 1.0, 0.001, 5);
    runTextFile(rootPath, txtPath, "gauss", "ak4PF");

    std::string header;
    auto lines = ReadJECFile(txtPath, header);
    if ((int)lines.size() != 36) {
        std::cout << "  SKIP  (line count wrong)\n";
        std::remove(rootPath); std::remove(txtPath);
        return;
    }

    Check(std::fabs(lines[0].etaLo  - (-5.191)) < kEps, "negative outer edge = -5.191");
    Check(std::fabs(lines[35].etaHi -   5.191)  < kEps, "positive outer edge =  5.191");

    std::remove(rootPath);
    std::remove(txtPath);
}

// ---- main ----

int main() {
    gROOT->SetBatch(true);

    std::cout << "=== TestTextFileWriter ===\n";

    TestFileStructure();
    TestEtaOrdering();
    TestSymmetricCorrections();
    TestUnityFallback_TooFewSlices();
    TestUnityFallback_EmptyBins();
    TestFitRoundTrip();
    TestEtaExtent();

    std::cout << "\n=== " << nPass << " passed, " << nFail << " failed ===\n";
    return nFail > 0 ? 1 : 0;
}
