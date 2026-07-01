#ifndef CONFIGCLI_H
#define CONFIGCLI_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace L2ConfigCli {

inline bool IsConfigArg( const std::string& arg ){
    return arg.rfind( "CONFIG=", 0 ) == 0;
}

inline void ApplyConfigArgument( int argc, char* argv[] ){
    for( int i = 1; i < argc; i++ ){
        const std::string arg = argv[i];
        if( !IsConfigArg( arg ) ) continue;
        const std::string path = arg.substr( 7 );
        if( path.empty() ){
            std::cerr << "ERROR: CONFIG= requires a TOML path\n";
            std::exit( 1 );
        }
        setenv( "L2RESIDUALS_CONFIG", path.c_str(), 1 );
    }
}

inline std::vector<std::string> PositionalArgs( int argc, char* argv[], int first = 1 ){
    std::vector<std::string> args;
    for( int i = first; i < argc; i++ ){
        const std::string arg = argv[i];
        if( IsConfigArg( arg ) ) continue;
        args.push_back( arg );
    }
    return args;
}

inline const char* ConfigUsage(){
    return " [CONFIG=path]";
}

}

#endif
