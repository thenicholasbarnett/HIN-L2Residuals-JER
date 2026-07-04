#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <sys/ioctl.h>
#include <unistd.h>

// ProgressBar pb("Label:", total, ProgressBar::kRandom);
// for(...){ doWork(); pb.Update(); }
// pb.Finish();

// Real terminal width via ioctl, falling back when stdout isn't a tty (e.g.
// Condor job output redirected to a log file).
inline int TerminalWidth(int fallback = 80) {
  struct winsize w;
  if (isatty(STDOUT_FILENO) && ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0 &&
      w.ws_col > 0) {
    return w.ws_col;
  }
  return fallback;
}

// Colors: kBlack, kGreen, kBlue, kRed, kPink, kPurple, kOrange, kYellow, kCyan, kWhite, kRandom

struct ProgressBar {

  enum Color {
    kBlack = 0,
    kGreen,
    kBlue,
    kRed,
    kPink,
    kPurple,
    kOrange,
    kYellow,
    kCyan,
    kWhite,
    kRandom
  };

  std::string label;
  int total;
  int current = 0;
  Color color;
  time_t startTime;

  ProgressBar(const std::string &label, int total, Color color = kRandom)
      : label(label), total(total), color(color), startTime(time(nullptr)) {
    if (this->color == kRandom) {
      static const Color all[] = {kGreen,  kBlue,   kRed,  kPink, kPurple,
                                  kOrange, kYellow, kCyan, kWhite};
      static const int N = (int)(sizeof(all) / sizeof(all[0]));
      static bool seeded = false;
      static Color lastColor = kRandom;
      static Color deck[N];
      static int remaining = 0;
      if (!seeded) {
        srand((unsigned int)(time(nullptr) ^ (unsigned int)getpid()));
        seeded = true;
      }
      if (remaining == 0) {
        for (int i = 0; i < N; i++) {
          deck[i] = all[i];
        }
        for (int i = N - 1; i > 0; i--) {
          int j = rand() % (i + 1);
          Color t = deck[i];
          deck[i] = deck[j];
          deck[j] = t;
        }
        if (deck[0] == lastColor && N > 1) {
          Color t = deck[0];
          deck[0] = deck[1];
          deck[1] = t;
        }
        remaining = N;
      }
      this->color = deck[N - remaining--];
      lastColor = this->color;
    }
  }

  void Update() {
    current++;
    Draw();
  }

  void Finish() {
    current = total;
    Draw();
    printf("\n\n");
    fflush(stdout);
  }

private:
  const char *AnsiColor() const {
    switch (color) {
    case kGreen:
      return "\033[32m";
    case kBlue:
      return "\033[34m";
    case kRed:
      return "\033[31m";
    case kPink:
      return "\033[35m";
    case kPurple:
      return "\033[95m";
    case kOrange:
      return "\033[33m"; // standard terminals render dim yellow as orange
    case kYellow:
      return "\033[93m"; // bright yellow
    case kCyan:
      return "\033[96m"; // bright cyan
    case kWhite:
      return "\033[97m";
    case kBlack:
    default:
      return "\033[30m";
    }
  }

  static constexpr int kMinBarWidth = 10;

  void Draw() const {
    if (total <= 0) {
      printf("\r  %s  [%s\033[0m]  0/0  100%%  00:00\033[K", label.c_str(),
             AnsiColor());
      fflush(stdout);
      return;
    }
    int pct = (current >= total) ? 100 : (current * 100 / total);

    char countBuf[32];
    snprintf(countBuf, sizeof(countBuf), "%d/%d", current, total);

    long elapsed = (long)(time(nullptr) - startTime);

    char rateBuf[16] = "";
    if (current > 0 && elapsed >= 1)
      snprintf(rateBuf, sizeof(rateBuf), "  %4.1f/s",
               (double)current / elapsed);

    char timeBuf[24] = "  --:-- ETA";
    if (current >= total && elapsed > 0)
      snprintf(timeBuf, sizeof(timeBuf), "  %02ld:%02ld", elapsed / 60,
               elapsed % 60);
    else if (current > 0 && pct >= 5 && elapsed >= 1) {
      long eta = elapsed * (long)(total - current) / (long)current;
      snprintf(timeBuf, sizeof(timeBuf), "  %02ld:%02ld ETA", eta / 60,
               eta % 60);
    }

    // everything on the line except the bar itself: "  label  [" + "]" +
    // "  " + counts + "  " + "100%" + rate + time, plus a 1-col safety margin
    int overhead = 2 + (int)label.size() + 2 + 2 + 1 + 2 +
                   (int)strlen(countBuf) + 2 + 4 + (int)strlen(rateBuf) +
                   (int)strlen(timeBuf) + 1;
    int barWidth = TerminalWidth() - overhead;
    if (barWidth < kMinBarWidth)
      barWidth = kMinBarWidth;

    int filled = (current >= total) ? barWidth : (current * barWidth / total);
    int empty = barWidth - filled;

    std::string bar, gap;
    for (int i = 0; i < filled; i++)
      bar += "\xe2\x96\x88"; // █
    for (int i = 0; i < empty; i++)
      gap += "\xe2\x96\x91"; // ░

    printf("\r  %s  [%s%s\033[90m%s\033[0m]  %s  %3d%%%s%s\033[K",
           label.c_str(), AnsiColor(), bar.c_str(), gap.c_str(), countBuf, pct,
           rateBuf, timeBuf);
    fflush(stdout);
  }
};

#endif
