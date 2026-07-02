#ifndef CLITOKENS_H
#define CLITOKENS_H

#include <cctype>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>

// Generic KEY=value CLI token parser shared by every compiled entry point
// (macros/*.C). There are no positional arguments anywhere in this repo's
// CLI: every value is a self-identifying KEY=value token, so nothing
// depends on argument order. This removes an entire class of silent
// argument-order bugs (e.g. swapping data.root/mc.root, or triggered/
// non-triggered residual files) — a malformed or misspelled token is a
// hard, immediate CLI error, never a silently-wrong run.
namespace L2Cli {

struct Tokens {
    std::map<std::string, std::string> kv;

    bool Has( const std::string& key ) const {
        return kv.find( key ) != kv.end();
    }

    std::string Get( const std::string& key, const std::string& def = "" ) const {
        auto it = kv.find( key );
        return ( it != kv.end() ) ? it->second : def;
    }

    // For required tokens (CONFIG, INPUT, OUTPUT, ...): errors and exits
    // immediately if the token wasn't given, printing usage.
    std::string Require( const std::string& key, const char* usage ) const {
        auto it = kv.find( key );
        if( it == kv.end() ){
            std::cerr << "ERROR: missing required " << key << "=... argument\n" << usage;
            std::exit( 1 );
        }
        return it->second;
    }
};

// Splits "KEY=value" on the first '='. The value may itself contain '='
// (e.g. a path never would, but this keeps the split unambiguous either
// way). Returns false if there's no '=', or the key part isn't a valid
// identifier ([A-Za-z_][A-Za-z0-9_]*) — that's what flags an argument as
// not being a token at all, rather than silently misparsing it.
inline bool SplitToken( const std::string& arg, std::string& key, std::string& value ){
    const size_t eq = arg.find( '=' );
    if( eq == std::string::npos || eq == 0 ) return false;
    key = arg.substr( 0, eq );
    if( std::isdigit( ( unsigned char )key[0] ) ) return false;
    for( char c : key ){
        if( !( std::isalnum( ( unsigned char )c ) || c == '_' ) ) return false;
    }
    value = arg.substr( eq + 1 );
    return true;
}

// Parses argv[1..] into KEY=value tokens. knownKeys is the exact set of
// keys this binary understands — any token whose key isn't in knownKeys,
// or any argument that isn't a valid KEY=value token at all, is a usage
// error. Nothing is ever silently ignored, so a typo (e.g. OUTPT=) fails
// loudly instead of quietly running with a missing/default value.
inline Tokens ParseTokens( int argc, char* argv[], const std::set<std::string>& knownKeys, const char* usage ){
    Tokens t;
    for( int i = 1; i < argc; i++ ){
        std::string key, value;
        if( !SplitToken( argv[i], key, value ) ){
            std::cerr << "ERROR: argument \"" << argv[i] << "\" is not a KEY=value token\n" << usage;
            std::exit( 1 );
        }
        if( knownKeys.find( key ) == knownKeys.end() ){
            std::cerr << "ERROR: unrecognized token \"" << key << "=\"\n" << usage;
            std::exit( 1 );
        }
        t.kv[key] = value;
    }
    return t;
}

}

#endif
