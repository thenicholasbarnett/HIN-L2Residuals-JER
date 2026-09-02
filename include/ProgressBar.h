#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>
#include <sys/ioctl.h>
#include <unistd.h>

// ProgressBar pb("Label:", total, ProgressBar::kRandom);
// for(...){ doWork(); pb.Update(); }
// pb.Finish();

// Real terminal width via ioctl, falling back when stdout isn't a tty (e.g.
// Condor job output redirected to a log file).
inline int TerminalWidth(int fallback = 80) {
  struct winsize w;
  int cols = fallback;
  if (isatty(STDOUT_FILENO) && ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0 &&
      w.ws_col > 0) {
    cols = w.ws_col;
  }
  // reserve one column -- writing all the way to the last column leaves the
  // terminal in a pending-wrap state that \r doesn't reliably clear, which
  // walks the bar down the screen over many redraws
  return cols - 1;
}

// Thrown by ProgressBar::Update() once a sample cap set via EnableSample()
// is hit. Unwinds straight out of a plotting loop regardless of how many
// nested loops or separate Plot* calls sit in between -- used by
// runPlotting's -sample flag to preview cosmetics without a full run.
struct SampleLimitReached {};

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
  int sampleLimit = -1;      // -1 == no cap; set via EnableSample()
  time_t sampleDeadline = 0; // absolute time(nullptr) cutoff; 0 == no cap

  // trailing-window rate/ETA state -- see Draw()'s comment
  std::vector<int> histCur;
  std::vector<time_t> histTime;
  std::string lastRateStr;
  std::string lastEtaStr;
  time_t lastEtaTime = 0;

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
    if (sampleLimit >= 0 &&
        (current >= sampleLimit || time(nullptr) >= sampleDeadline))
      throw SampleLimitReached{};
  }

  // Caps this bar to whichever comes first: maxCount total Update() calls,
  // or maxSeconds of wall-clock time since construction. Update() throws
  // SampleLimitReached the instant either limit is hit -- the caller must
  // wrap its plotting loop in a try/catch to stop early.
  void EnableSample(int maxCount, int maxSeconds) {
    sampleLimit = maxCount;
    sampleDeadline = startTime + maxSeconds;
  }

  // Call once per candidate item, before doing its expensive work, and skip
  // (continue) the iteration if this returns false. Unlike the maxCount cap
  // in Update() -- which only stops the run once that many items have
  // already been fully built, biasing a sample toward whatever comes first
  // in the plotting loops' fixed nesting order (one cone, one family, low
  // bin indices) -- this is an independent Bernoulli trial per candidate
  // with keep-probability sampleLimit/total, so a sample touches every
  // cone/eta range/pT slice instead of just a prefix of them.
  //
  // Deliberately not an exact evenly-spaced/stride selection: that needs
  // `seen` to land exactly on `total`, which breaks the moment one candidate
  // yields more than one Update() (e.g. adist's 3 sub-plots per eta bin) --
  // `seen` would only ever reach total/3, never crossing the stride's later
  // thresholds. A per-candidate independent probability doesn't care how
  // many Update() calls one kept candidate produces: since `total` already
  // has that multiplicity baked in (it's summed straight from the same
  // Count*Plots functions used to size the bar), the expected number of
  // kept plots works out to sampleLimit regardless. Also enforces the
  // wall-clock deadline here, since a long run of skipped (cheap)
  // candidates never reaches Update() to check it there.
  bool ShouldKeep() {
    if (sampleLimit < 0)
      return true;
    if (time(nullptr) >= sampleDeadline)
      throw SampleLimitReached{};
    if (total <= 0)
      return true;
    return (rand() % total) < sampleLimit;
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

  void Draw() {
    if (total <= 0) {
      printf("\r  %s  [%s\033[0m]  0/0  100%%  00:00\033[K", label.c_str(),
             AnsiColor());
      fflush(stdout);
      return;
    }
    int pct = (current >= total) ? 100 : (current * 100 / total);

    char countBuf[32];
    snprintf(countBuf, sizeof(countBuf), "%d/%d", current, total);

    time_t now = time(nullptr);
    long elapsed = (long)(now - startTime);

    char timeBuf[24] = "  --:-- ETA";
    if (current >= total)
      timeBuf[0] = '\0';

    if (current > 0) {
      // rate is a trailing window over the last ~5% of total, not a
      // since-start running average -- a real slowdown would get smoothed
      // away by a full-run average instead of showing up in the number
      histCur.push_back(current);
      histTime.push_back(now);
      int window = total * 5 / 100;
      if (window < 1)
        window = 1;
      int threshold = current - window;
      while (histCur.size() > 1 && histCur.front() < threshold) {
        histCur.erase(histCur.begin());
        histTime.erase(histTime.begin());
      }
      long winElapsed = (long)(now - histTime.front());
      int winDelta = current - histCur.front();
      if (winElapsed >= 1 && winDelta > 0) {
        char rateBuf[16];
        snprintf(rateBuf, sizeof(rateBuf), "  %4.1f/s",
                 (double)winDelta / (double)winElapsed);
        lastRateStr = rateBuf;
      }

      if (elapsed >= 1) {
        if (current >= total) {
          snprintf(timeBuf, sizeof(timeBuf), "  %02ld:%02ld", elapsed / 60,
                   elapsed % 60);
        } else if (pct >= 5) {
          // ETA text only refreshes once a second -- the windowed rate
          // above still updates every call, this just stops the countdown
          // from jittering with it
          if (lastEtaTime == 0 || now - lastEtaTime >= 1) {
            long eta = (winElapsed >= 1 && winDelta > 0)
                          ? winElapsed * (long)(total - current) / winDelta
                          : elapsed * (long)(total - current) / (long)current;
            char etaBuf[24];
            snprintf(etaBuf, sizeof(etaBuf), "  %02ld:%02ld ETA", eta / 60,
                     eta % 60);
            lastEtaStr = etaBuf;
            lastEtaTime = now;
          }
          snprintf(timeBuf, sizeof(timeBuf), "%s", lastEtaStr.c_str());
        }
      }
    }

    // drop rate, then time, then label as the terminal narrows, so the bar
    // is last to go; below kMinBarWidth, drop the bar too
    const char *useLabel = label.c_str();
    const char *useRate = lastRateStr.c_str();
    const char *useTime = timeBuf;
    int termWidth = TerminalWidth();
    int barWidth = 0;
    for (int stage = 0; stage < 4; stage++) {
      if (stage == 1)
        useRate = "";
      if (stage == 2)
        useTime = "";
      if (stage == 3)
        useLabel = "";
      int overhead = 2 + (int)strlen(useLabel) + 2 + 2 + 1 + 2 +
                     (int)strlen(countBuf) + 2 + 4 + (int)strlen(useRate) +
                     (int)strlen(useTime) + 1;
      barWidth = termWidth - overhead;
      if (barWidth >= kMinBarWidth)
        break;
    }

    if (barWidth < kMinBarWidth) {
      printf("\r  %s  %3d%%%s%s\033[K", countBuf, pct, useRate, useTime);
      fflush(stdout);
      return;
    }

    int filled = (current >= total) ? barWidth : (current * barWidth / total);
    int empty = barWidth - filled;

    std::string bar, gap;
    for (int i = 0; i < filled; i++)
      bar += "\xe2\x96\x88"; // █
    for (int i = 0; i < empty; i++)
      gap += "\xe2\x96\x91"; // ░

    printf("\r  %s  [%s%s\033[90m%s\033[0m]  %s  %3d%%%s%s\033[K",
           useLabel, AnsiColor(), bar.c_str(), gap.c_str(), countBuf, pct,
           useRate, useTime);
    fflush(stdout);
  }
};

#endif
