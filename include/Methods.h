#ifndef METHODS_H
#define METHODS_H

// Shared method identifiers for every fitting/extraction routine offering
// Gauss/doubleGauss/trunc90/trunc95/crystalball as parallel alternatives
// (CalibrationExtractor.cxx's A-distribution fitting, ResponseExtractor.cxx's
// JES/JER extraction). Kept separate from plotting/Style.h's
// presentation-only arrays (labels/colors/marker styles) so non-plotting
// extraction code doesn't need to pull in plotting headers just for these
// two identifiers -- previously duplicated independently in both places.

static const char *const kMethodKeys[] = {"gauss", "doubleGauss", "trunc90",
                                          "trunc95", "crystalball"};
static constexpr int kNMethods = 5;

#endif
