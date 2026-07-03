#ifndef CLITOKENS_H
#define CLITOKENS_H

#include <cctype>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

// Generic command parser shared by every compiled entry point (macros/*.C).
// It supports this repo's original KEY=value shell tokens and the JetMET-style
// interface:
//   runX params.config -input file.root -output out.root
// where params.config contains lines like:
//   input = file.root
//   flags = finals etasym
// Command-line values override earlier config-file values.
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

inline bool IsIdentifier( const std::string& key ){
    if( key.empty() ) return false;
    if( std::isdigit( ( unsigned char )key[0] ) ) return false;
    for( char c : key ){
        if( !( std::isalnum( ( unsigned char )c ) || c == '_' ) ) return false;
    }
    return true;
}

inline void StoreToken( Tokens& t, std::string key, const std::string& value,
                       const std::set<std::string>& knownKeys, const char* usage,
                       const std::string& source ){
    key = ToUpper( Trim( key ) );
    if( !IsIdentifier( key ) ){
        std::cerr << "ERROR: invalid token key \"" << key << "\" in " << source << "\n" << usage;
        std::exit( 1 );
    }
    if( knownKeys.find( key ) == knownKeys.end() ){
        std::cerr << "ERROR: unrecognized token \"" << key << "=\" in " << source << "\n" << usage;
        std::exit( 1 );
    }
    t.kv[key] = Trim( value );
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
    if( !IsIdentifier( key ) ) return false;
    key = ToUpper( key );
    value = Trim( arg.substr( eq + 1 ) );
    return true;
}

inline std::vector<std::string> TokenizeConfigFile( const std::string& path ){
    std::ifstream in( path );
    if( !in.is_open() ){
        std::cerr << "ERROR: cannot open configuration argument file: " << path << "\n";
        std::exit( 1 );
    }

    std::vector<std::string> tokens;
    std::string cur;
    bool inString = false;
    bool inComment = false;
    char c = '\0';
    while( in.get( c ) ){
        if( !inString && ( c == '$' || c == '#' ) ){
            inComment = true;
            continue;
        }
        if( inComment ){
            if( c == '\n' ) inComment = false;
            continue;
        }
        if( c == '"' ){
            inString = !inString;
            continue;
        }
        if( !inString && c == '=' ){
            if( !cur.empty() ){ tokens.push_back( cur ); cur.clear(); }
            tokens.push_back( "=" );
            continue;
        }
        if( !inString && std::isspace( ( unsigned char )c ) ){
            if( !cur.empty() ){ tokens.push_back( cur ); cur.clear(); }
            continue;
        }
        cur.push_back( c );
    }
    if( !cur.empty() ) tokens.push_back( cur );
    return tokens;
}

inline void ParseConfigFile( Tokens& t, const std::string& path,
                            const std::set<std::string>& knownKeys, const char* usage ){
    std::deque<std::string> tokens;
    for( const auto& tok : TokenizeConfigFile( path ) ) tokens.push_back( tok );

    while( !tokens.empty() ){
        std::string key = tokens.front();
        tokens.pop_front();
        if( key == "include" || key == "INCLUDE" ){
            if( tokens.empty() ){
                std::cerr << "ERROR: include without file in " << path << "\n" << usage;
                std::exit( 1 );
            }
            std::string included = tokens.front();
            tokens.pop_front();
            ParseConfigFile( t, included, knownKeys, usage );
            continue;
        }
        if( tokens.empty() || tokens.front() != "=" ){
            std::cerr << "ERROR: expected \"" << key << " = value\" in " << path << "\n" << usage;
            std::exit( 1 );
        }
        tokens.pop_front(); // =
        if( tokens.empty() ){
            std::cerr << "ERROR: missing value for \"" << key << "\" in " << path << "\n" << usage;
            std::exit( 1 );
        }

        std::string value = tokens.front();
        tokens.pop_front();
        while( !tokens.empty() ){
            if( tokens.size() >= 2 && tokens[1] == "=" ) break;
            value += " " + tokens.front();
            tokens.pop_front();
        }
        StoreToken( t, key, value, knownKeys, usage, path );
    }
}

inline std::string JoinValues( const std::vector<std::string>& values ){
    std::string out;
    for( const auto& v : values ){
        if( !out.empty() ) out += " ";
        out += v;
    }
    return out;
}

// Parses argv[1..]. knownKeys is the exact set of keys this binary understands.
// Accepted forms:
//   CONFIG=cfg/2024ppRef.toml
//   -config cfg/2024ppRef.toml
//   params.config        (before dashed options; contents use key = value)
inline Tokens ParseTokens( int argc, char* argv[], const std::set<std::string>& knownKeys, const char* usage ){
    Tokens t;
    bool allowConfigFiles = true;
    for( int i = 1; i < argc; i++ ){
        const std::string arg = argv[i];
        std::string key, value;

        if( SplitToken( arg, key, value ) ){
            StoreToken( t, key, value, knownKeys, usage, "command line" );
            continue;
        }

        if( !arg.empty() && arg[0] == '-' ){
            allowConfigFiles = false;
            key = arg.substr( 1 );
            if( key.empty() || !IsIdentifier( key ) ){
                std::cerr << "ERROR: invalid option \"" << arg << "\"\n" << usage;
                std::exit( 1 );
            }

            std::vector<std::string> values;
            while( i + 1 < argc ){
                const std::string next = argv[i + 1];
                std::string dummyKey, dummyValue;
                if( !next.empty() && next[0] == '-' ) break;
                if( SplitToken( next, dummyKey, dummyValue ) ) break;
                values.push_back( next );
                i++;
            }
            if( values.empty() ){
                std::cerr << "ERROR: option \"" << arg << "\" requires a value\n" << usage;
                std::exit( 1 );
            }
            StoreToken( t, key, JoinValues( values ), knownKeys, usage, "command line" );
            continue;
        }

        if( allowConfigFiles ){
            ParseConfigFile( t, arg, knownKeys, usage );
            continue;
        }

        std::cerr << "ERROR: argument \"" << arg
                  << "\" is not KEY=value, -key value, or an initial config file\n"
                  << usage;
        std::exit( 1 );
    }
    return t;
}

inline void PrintUsageForms( std::ostream& os ){
    os << "Argument forms:\n"
       << "  KEY=value                      original shell-token form\n"
       << "  -key value                     JetMET-style command-line form\n"
       << "  params.config -key override    config file plus command-line overrides\n"
       << "Config file syntax:\n"
       << "  key = value\n"
       << "  flags = finals etasym methods\n";
}

inline void RequireNoUnknownPositionalExample(){
    // Link-time no-op anchor for old cling sessions that may cache this header.
}

}

#endif
