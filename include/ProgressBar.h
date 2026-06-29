#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <unistd.h>

// ProgressBar pb("Label:", total, ProgressBar::kRandom);  // color optional, default kBlack
// for(...){ doWork(); pb.Update(); }
// pb.Finish();

// Colors: kBlack, kGreen, kBlue, kRed, kPink, kPurple, kOrange, kYellow, kCyan, kWhite, kRandom

struct ProgressBar {

    enum Color {
        kBlack  = 0,
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
    int width = 40;
    Color color;

    ProgressBar(const std::string& label, int total, Color color = kRandom)
        : label(label), total(total), color(color)
    {
        if(this->color == kRandom){
            static const Color choices[] = {kGreen, kBlue, kRed, kPink, kPurple, kOrange, kYellow, kCyan, kWhite};
            static constexpr int nChoices = (int)(sizeof(choices) / sizeof(choices[0]));
            static bool seeded = false;
            static Color lastColor = kRandom;  // kRandom not in choices, so first pick is unconstrained
            if(!seeded){ srand((unsigned int)(time(nullptr) ^ (unsigned int)getpid())); seeded = true; }
            Color c; int tries = 0;
            do { c = choices[rand() % nChoices]; } while(c == lastColor && ++tries < 20);
            lastColor = c;
            this->color = c;
        }
    }

    void Update(){current++; Draw();}
    void Finish(){current = total; Draw(); printf("\n"); fflush(stdout);}

private:
    const char* AnsiColor() const {
        switch(color){
            case kGreen: return "\033[32m";
            case kBlue: return "\033[34m";
            case kRed: return "\033[31m";
            case kPink: return "\033[35m";
            case kPurple: return "\033[95m";
            case kOrange: return "\033[33m";  // standard terminals render dim yellow as orange
            case kYellow: return "\033[93m";  // bright yellow
            case kCyan: return "\033[96m";  // bright cyan
            case kWhite: return "\033[97m";
            case kBlack:
            default: return "\033[30m";
        }
    }

    void Draw() const {
        int filled = (total > 0 && current >= total) ? width : (current * width / total);
        int empty = width - filled;
        int pct = (total > 0 && current >= total) ? 100 : (current * 100 / total);

        std::string bar, gap;
        for(int i = 0; i < filled; i++) bar += "\xe2\x96\x88";  // █
        for(int i = 0; i < empty;  i++) gap += "\xe2\x96\x91";  // ░

        printf("\r  %-24s [%s%s\033[90m%s\033[0m] %d/%d (%d%%)",
               label.c_str(),
               AnsiColor(), bar.c_str(),
               gap.c_str(),
               current, total, pct);
        fflush(stdout);
    }
};

#endif
