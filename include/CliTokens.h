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

inline std::string Trim( const std::string& s ){
    const size_t b = s.find_first_not_of( " \t" );
    if( b == std::string::npos ) return "";
    const size_t e = s.find_last_not_of( " \t" );
    return s.substr( b, e - b + 1 );
}

inline std::string ToUpper( const std::string& s ){
    std::string r = s;
    for( char& c : r ) c = ( char )std::toupper( ( unsigned char )c );
    return r;
}

// Splits "KEY=value" on the first '='. The value may itself contain '='
// (e.g. a path never would, but this keeps the split unambiguous either
// way). Returns false if there's no '=', or the key part (after trimming
// whitespace) isn't a valid identifier ([A-Za-z_][A-Za-z0-9_]*) — that's
// what flags an argument as not being a token at all, rather than silently
// misparsing it. The key is case-normalized to uppercase and both key and
// value have surrounding whitespace trimmed, so "config = cfg/x.toml" and
// "CONFIG=cfg/x.toml" parse identically -- this is about tolerating
// formatting, not about accepting a misspelled key: "COFNIG=" still fails,
// exactly as it should. Note a genuinely *unquoted* "KEY = value" typed at
// a shell is three separate argv words ("KEY", "=", "value"), which this
// can't reassemble -- only a single already-quoted argument benefits from
// the whitespace tolerance.
inline bool SplitToken( const std::string& arg, std::string& key, std::string& value ){
    const size_t eq = arg.find( '=' );
    if( eq == std::string::npos || eq == 0 ) return false;
    key = Trim( arg.substr( 0, eq ) );
    if( key.empty() ) return false;
    if( std::isdigit( ( unsigned char )key[0] ) ) return false;
    for( char c : key ){
        if( !( std::isalnum( ( unsigned char )c ) || c == '_' ) ) return false;
    }
    key = ToUpper( key );
    value = Trim( arg.substr( eq + 1 ) );
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
