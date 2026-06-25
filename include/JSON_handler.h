#ifndef JSON_HANDLER_H
#define JSON_HANDLER_H

#include <string>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;


// Most recent version of the JSON file is in the public EOS HI area:
// /eos/cms/store/group/phys_heavyions/2026_DCS_JSON_Preliminary/DCS2026_latest.json

class JSON_handler {
private:
    json dcsJson;
    void loadJSON(const std::string &jsonFilePath) {
        std::ifstream dcsFile(jsonFilePath);
        if (!dcsFile.is_open()) {
            throw std::runtime_error("Error: file JSON impossible to be opened: " + jsonFilePath);
        }
        dcsFile >> dcsJson;
        dcsFile.close();

        // colore verde ANSI
        std::cout << "\033[1;32mLoaded DCS JSON from: " << jsonFilePath << "\033[0m" << std::endl;
    }

public:

    /// Default → file standard on EOS
    JSON_handler() {
        std::string defaultPath = "/eos/cms/store/group/phys_heavyions/2026_DCS_JSON_Preliminary/DCS2026_latest.json";
        loadJSON(defaultPath);
    }

    // Custom constructor
    JSON_handler(const std::string &jsonFilePath){ loadJSON(jsonFilePath);}

    //For testing if run is included before attempting to apply json
    //Useful for provisional jsons being updated intermittently
    bool isGoodRun(unsigned int run) const {
        std::string runStr = std::to_string(run);
		if (!dcsJson.contains(runStr)){
		  std::cout << __PRETTY_FUNCTION__ << " WARNING: JSON does not contain run '" << run << "'. return false" << std::endl;
		  return false;
		}
    	return false;
    }

    bool isGood(unsigned int run, unsigned int lumi) const {
        std::string runStr = std::to_string(run);

        if (!dcsJson.contains(runStr)) return false;
        for (const auto &interval : dcsJson[runStr]) if (lumi >= interval[0] && lumi <= interval[1]) return true;

        return false;
    }
};

#endif


// How to use it:

// #include "JSON_handler.h"
// JSON_handler dcs("/path/al/file/dcs.json");
// OR
// JSON_handler dcs;

// if (dcs.isGood(run, lumi)) {
//     // Do your stuff
// }



