#ifndef COLORS_H
#define COLORS_H

#include "TColor.h"

// "Sailing Boats Returning to Yabase", Lake Biwa, 1835, Utagawa Hiroshige https://github.com/BlakeRMills/MetBrewer?tab=readme-ov-file#hiroshige

inline Color_t HiroshigeLightRed() {
  static const Color_t c = TColor::GetColor("#e76254");
  return c;
}
inline Color_t HiroshigeOrange() {
  static const Color_t c = TColor::GetColor("#ef8a47");
  return c;
}
inline Color_t HiroshigeLightOrange() {
  static const Color_t c = TColor::GetColor("#f7aa58");
  return c;
}
inline Color_t HiroshigeYellow() {
  static const Color_t c = TColor::GetColor("#ffd06f");
  return c;
}
inline Color_t HiroshigeIceBlue() {
  static const Color_t c = TColor::GetColor("#aadce0");
  return c;
}
inline Color_t HiroshigeLightBlue() {
  static const Color_t c = TColor::GetColor("#72bcd5");
  return c;
}
inline Color_t HiroshigeBlue() {
  static const Color_t c = TColor::GetColor("#528fad");
  return c;
}
inline Color_t HiroshigeGrayBlue() {
  static const Color_t c = TColor::GetColor("#376795");
  return c;
}
inline Color_t HiroshigeNightBlue() {
  static const Color_t c = TColor::GetColor("#1e466e");
  return c;
}

// "Mäda Primavesi", 1912, Gustav Klimt https://github.com/BlakeRMills/MetBrewer?tab=readme-ov-file#klimt

inline Color_t KlimtPink() {
  static const Color_t c = TColor::GetColor("#df9ed4");
  return c;
}
inline Color_t KlimtRed() {
  static const Color_t c = TColor::GetColor("#c93f55");
  return c;
}
inline Color_t KlimtYellow() {
  static const Color_t c = TColor::GetColor("#eacc62");
  return c;
}
inline Color_t KlimtGreen() {
  static const Color_t c = TColor::GetColor("#469d76");
  return c;
}
inline Color_t KlimtBlue() {
  static const Color_t c = TColor::GetColor("#3c4b99");
  return c;
}
inline Color_t KlimtPurple() {
  static const Color_t c = TColor::GetColor("#924099");
  return c;
}

#endif
