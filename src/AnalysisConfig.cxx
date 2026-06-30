#include "AnalysisConfig.h"

#include "TSystem.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

std::string Trim( const std::string& s ){
    size_t first = 0;
    while( first < s.size() && std::isspace( ( unsigned char )s[first] ) ) first++;
    size_t last = s.size();
    while( last > first && std::isspace( ( unsigned char )s[last - 1] ) ) last--;
    return s.substr( first, last - first );
}

std::string StripComment( const std::string& line ){
    bool inString = false;
    bool escape = false;
    for( size_t i = 0; i < line.size(); i++ ){
        char c = line[i];
        if( inString ){
            if( escape ) escape = false;
            else if( c == '\\' ) escape = true;
            else if( c == '"' ) inString = false;
        } else {
            if( c == '"' ) inString = true;
            else if( c == '#' ) return line.substr( 0, i );
        }
    }
    return line;
}

std::string ParseString( const std::string& value ){
    std::string v = Trim( value );
    if( v.size() < 2 || v.front() != '"' || v.back() != '"' )
        throw std::runtime_error( "Expected quoted string, got: " + value );
    return v.substr( 1, v.size() - 2 );
}

std::vector<std::string> ParseStringArray( const std::string& value ){
    std::vector<std::string> out;
    std::string v = Trim( value );
    if( v.size() < 2 || v.front() != '[' || v.back() != ']' )
        throw std::runtime_error( "Expected string array, got: " + value );

    bool inString = false;
    bool escape = false;
    std::string token;
    for( size_t i = 1; i + 1 < v.size(); i++ ){
        char c = v[i];
        if( inString ){
            token += c;
            if( escape ) escape = false;
            else if( c == '\\' ) escape = true;
            else if( c == '"' ){
                inString = false;
                out.push_back( token.substr( 1, token.size() - 2 ) );
                token.clear();
            }
        } else if( c == '"' ){
            inString = true;
            token = "\"";
        }
    }
    return out;
}

std::vector<std::vector<std::string>> ParseNestedStringArray( const std::string& value ){
    std::vector<std::vector<std::string>> out;
    std::string v = Trim( value );
    if( v.size() < 2 || v.front() != '[' || v.back() != ']' )
        throw std::runtime_error( "Expected nested string array, got: " + value );

    int depth = 0;
    bool inString = false;
    bool escape = false;
    size_t start = std::string::npos;
    for( size_t i = 1; i + 1 < v.size(); i++ ){
        char c = v[i];
        if( inString ){
            if( escape ) escape = false;
            else if( c == '\\' ) escape = true;
            else if( c == '"' ) inString = false;
            continue;
        }
        if( c == '"' ){
            inString = true;
        } else if( c == '[' ){
            if( depth == 0 ) start = i;
            depth++;
        } else if( c == ']' ){
            depth--;
            if( depth == 0 && start != std::string::npos ){
                out.push_back( ParseStringArray( v.substr( start, i - start + 1 ) ) );
                start = std::string::npos;
            }
        }
    }
    return out;
}

TString ToTString( const std::string& s ){
    return TString( s.c_str() );
}

std::vector<TString> ToTStrings( const std::vector<std::string>& values ){
    std::vector<TString> out;
    out.reserve( values.size() );
    for( const auto& value : values ) out.push_back( ToTString( value ) );
    return out;
}

void Require( bool condition, const std::string& message ){
    if( !condition ) throw std::runtime_error( "Invalid analysis config: " + message );
}

}  // namespace

std::string DefaultConfigPath(){
    const char* envConfig = std::getenv( "L2RESIDUALS_CONFIG" );
    if( envConfig && envConfig[0] ) return envConfig;

    const char* envHome = std::getenv( "L2RESIDUALS_HOME" );
    if( envHome && envHome[0] ) return std::string( envHome ) + "/cfg/2024ppRef.toml";

    return "cfg/2024ppRef.toml";
}

AnalysisConfig LoadAnalysisConfig( const std::string& path ){
    const std::string configPath = path.empty() ? DefaultConfigPath() : path;
    std::ifstream in( configPath );
    if( !in.is_open() )
        throw std::runtime_error( "Cannot open analysis config: " + configPath );

    AnalysisConfig cfg;
    std::string section;
    std::string pendingKey;
    std::string pendingValue;
    int bracketDepth = 0;
    std::string rawLine;

    auto assign = [&]( const std::string& fullKey, const std::string& rawValue ){
        std::string value = Trim( rawValue );
        if( fullKey == "paths.veto_map" ) cfg.vetoMapPath = ParseString( value );
        else if( fullKey == "paths.golden_json" ) cfg.jsonPath = ToTString( ParseString( value ) );
        else if( fullKey == "trees.hi" ) cfg.hiTreePath = ToTString( ParseString( value ) );
        else if( fullKey == "trees.skim" ) cfg.skimTreePath = ToTString( ParseString( value ) );
        else if( fullKey == "trees.trigger" ) cfg.trigTreePath = ToTString( ParseString( value ) );
        else if( fullKey == "trees.jets" ) cfg.jetTreePaths = ToTStrings( ParseStringArray( value ) );
        else if( fullKey == "trigger.branch" ) cfg.hltJ80Branch = ToTString( ParseString( value ) );
        else if( fullKey == "trigger.threshold" ) cfg.hltJ80Thresh = std::stof( value );
        else if( fullKey == "trigger.cone" ) cfg.trigCone = ToTString( ParseString( value ) );
        else if( fullKey == "cones.labels" ) cfg.coneLabels = ToTStrings( ParseStringArray( value ) );
        else if( fullKey == "jec.files" ) cfg.jecFilesPerCone = ParseNestedStringArray( value );
        else throw std::runtime_error( "Unknown analysis config key: " + fullKey );
    };

    while( std::getline( in, rawLine ) ){
        std::string line = Trim( StripComment( rawLine ) );
        if( line.empty() ) continue;

        if( !pendingKey.empty() ){
            pendingValue += " " + line;
            bracketDepth += std::count( line.begin(), line.end(), '[' );
            bracketDepth -= std::count( line.begin(), line.end(), ']' );
            if( bracketDepth <= 0 ){
                assign( pendingKey, pendingValue );
                pendingKey.clear();
                pendingValue.clear();
            }
            continue;
        }

        if( line.front() == '[' && line.back() == ']' ){
            section = Trim( line.substr( 1, line.size() - 2 ) );
            continue;
        }

        size_t eq = line.find( '=' );
        if( eq == std::string::npos )
            throw std::runtime_error( "Invalid analysis config line: " + rawLine );

        std::string key = Trim( line.substr( 0, eq ) );
        std::string fullKey = section.empty() ? key : section + "." + key;
        std::string value = Trim( line.substr( eq + 1 ) );
        bracketDepth = std::count( value.begin(), value.end(), '[' ) -
                       std::count( value.begin(), value.end(), ']' );
        if( bracketDepth > 0 ){
            pendingKey = fullKey;
            pendingValue = value;
        } else {
            assign( fullKey, value );
        }
    }

    Require( !cfg.jecFilesPerCone.empty(), "jec.files is empty" );
    Require( !cfg.coneLabels.empty(), "cones.labels is empty" );
    Require( cfg.jecFilesPerCone.size() == cfg.coneLabels.size(), "jec.files and cones.labels have different lengths" );
    Require( cfg.jetTreePaths.size() == cfg.coneLabels.size(), "trees.jets and cones.labels have different lengths" );
    Require( !cfg.trigCone.IsNull(), "trigger.cone is missing" );
    return cfg;
}

const AnalysisConfig& Config(){
    static const AnalysisConfig cfg = LoadAnalysisConfig();
    return cfg;
}
