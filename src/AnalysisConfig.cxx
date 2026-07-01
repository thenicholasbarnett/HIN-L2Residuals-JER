#include "AnalysisConfig.h"

#include "toml.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

namespace fs = std::filesystem;

bool IsAbsolutePath( const std::string& path ){
    return !path.empty() && fs::path( path ).is_absolute();
}

std::string RepoRoot(){
    const char* envHome = std::getenv( "L2RESIDUALS_HOME" );
    if( envHome && envHome[0] ) return envHome;

#ifdef L2RESIDUALS_SOURCE_DIR
    return L2RESIDUALS_SOURCE_DIR;
#else
    return ".";
#endif
}

std::string ResolvePath( const std::string& path, const std::string& root ){
    if( path.empty() || IsAbsolutePath( path ) ) return path;
    return ( fs::path( root ) / path ).lexically_normal().string();
}

std::string ResolveConfigPath( const std::string& path ){
    if( path.empty() ) return path;
    return ResolvePath( path, RepoRoot() );
}

std::string ConfigRoot( const std::string& configPath ){
    const char* envHome = std::getenv( "L2RESIDUALS_HOME" );
    if( envHome && envHome[0] ) return envHome;

    if( IsAbsolutePath( configPath ) ){
        fs::path configDir = fs::path( configPath ).parent_path();
        if( configDir.filename() == "cfg" ) return configDir.parent_path().string();
        return configDir.string();
    }

    return RepoRoot();
}

TString ResolvePath( const TString& path, const std::string& root ){
    return TString( ResolvePath( std::string( path.Data() ), root ).c_str() );
}

}

std::string DefaultConfigPath(){
    const char* envConfig = std::getenv( "L2RESIDUALS_CONFIG" );
    if( envConfig && envConfig[0] ) return ResolveConfigPath( envConfig );

    const char* envHome = std::getenv( "L2RESIDUALS_HOME" );
    if( envHome && envHome[0] ) return std::string( envHome ) + "/cfg/2024ppRef.toml";

    return ResolvePath( std::string( "cfg/2024ppRef.toml" ), RepoRoot() );
}

AnalysisConfig LoadAnalysisConfig( const std::string& path ){
    const bool explicitConfig = !path.empty() || ( std::getenv( "L2RESIDUALS_CONFIG" ) && std::getenv( "L2RESIDUALS_CONFIG" )[0] );
    const std::string configPath = path.empty() ? DefaultConfigPath() : ResolveConfigPath( path );
    const auto doc = toml::parse_file( configPath );
    const std::string repoRoot = ConfigRoot( configPath );

    AnalysisConfig cfg;
    cfg.configPath = configPath;
    cfg.repoRoot = repoRoot;
    cfg.usedDefaultConfig = !explicitConfig;

    cfg.vetoMapPath  = ResolvePath( doc["paths"]["veto_map"].value_or( std::string{} ), repoRoot );
    cfg.jsonPath     = ResolvePath( TString( doc["paths"]["golden_json"].value_or( std::string{} ).c_str() ), repoRoot );

    cfg.vetoMapHist  = doc["jet_id"]["veto_map_histogram"].value_or( std::string{} );

    cfg.hiTreePath   = TString( doc["trees"]["hi"].value_or( std::string{} ).c_str() );
    cfg.skimTreePath = TString( doc["trees"]["skim"].value_or( std::string{} ).c_str() );
    cfg.trigTreePath = TString( doc["trees"]["trigger"].value_or( std::string{} ).c_str() );
    cfg.filterBranch = TString( doc["trees"]["filter"].value_or( std::string{} ).c_str() );

    for( const auto& v : *doc["trees"]["jets"].as_array() )
        cfg.jetTreePaths.push_back( TString( v.value_or( std::string{} ).c_str() ) );

    cfg.hltJ80Branch = TString( doc["trigger"]["branch"].value_or( std::string{} ).c_str() );
    cfg.hltJ80Thresh = doc["trigger"]["threshold"].value_or( 0.0f );
    cfg.trigCone     = TString( doc["trigger"]["cone"].value_or( std::string{} ).c_str() );

    for( const auto& v : *doc["cones"]["labels"].as_array() )
        cfg.coneLabels.push_back( TString( v.value_or( std::string{} ).c_str() ) );

    for( const auto& row : *doc["jec"]["files"].as_array() ){
        std::vector<std::string> files;
        for( const auto& v : *row.as_array() )
            files.push_back( ResolvePath( v.value_or( std::string{} ), repoRoot ) );
        cfg.jecFilesPerCone.push_back( std::move( files ) );
    }

    // Optional — absent entirely, or empty per-cone, both mean "no residual
    // correction chained in for this cone yet."
    if( auto* residualArr = doc["jec"]["residual_files"].as_array() ){
        for( const auto& row : *residualArr ){
            std::vector<std::string> files;
            for( const auto& v : *row.as_array() )
                files.push_back( ResolvePath( v.value_or( std::string{} ), repoRoot ) );
            cfg.residualFilesPerCone.push_back( std::move( files ) );
        }
    }

    for( const auto& v : *doc["binning"]["ptavg_edges"].as_array() )
        cfg.ptavgEdges.push_back( v.value_or( 0.0f ) );

    // min_jet_pt is a single global floor used both for the Step 1 inclusive-
    // jet histogram and as the dijet subleading-jet validity cut — the third
    // jet's pT (used for alpha) is deliberately NOT subject to this.
    cfg.minJetPt        = doc["cuts"]["min_jet_pt"].value_or( 0.0f );
    cfg.dphiCut         = doc["cuts"]["dphi"].value_or( 0.0f );
    cfg.maxAbsA         = doc["cuts"]["max_abs_a"].value_or( 0.0f );
    cfg.minEntriesPerBin = doc["cuts"]["min_entries_per_bin"].value_or( 0 );
    cfg.residualGausFitHalfWidth = doc["cuts"]["residual_gaus_fit_half_width"].value_or( 0.5 );
    cfg.residualAlphaFitHi = doc["cuts"]["residual_alpha_fit_hi"].value_or( 0.31 );

    if( cfg.jecFilesPerCone.empty() ) throw std::runtime_error( "jec.files is empty" );
    if( cfg.coneLabels.empty() )      throw std::runtime_error( "cones.labels is empty" );
    if( cfg.jecFilesPerCone.size() != cfg.coneLabels.size() )
        throw std::runtime_error( "jec.files and cones.labels have different lengths" );
    if( !cfg.residualFilesPerCone.empty() && cfg.residualFilesPerCone.size() != cfg.coneLabels.size() )
        throw std::runtime_error( "jec.residual_files and cones.labels have different lengths" );
    if( cfg.jetTreePaths.size() != cfg.coneLabels.size() )
        throw std::runtime_error( "trees.jets and cones.labels have different lengths" );
    if( cfg.trigCone.IsNull() ) throw std::runtime_error( "trigger.cone is missing" );
    if( cfg.ptavgEdges.size() < 2 ) throw std::runtime_error( "binning.ptavg_edges needs at least 2 edges" );
    for( size_t i = 1; i < cfg.ptavgEdges.size(); i++ )
        if( cfg.ptavgEdges[i] <= cfg.ptavgEdges[i - 1] )
            throw std::runtime_error( "binning.ptavg_edges must be strictly ascending" );

    return cfg;
}

const AnalysisConfig& Config(){
    static const AnalysisConfig cfg = LoadAnalysisConfig();
    return cfg;
}

void PrintConfigSummary( const AnalysisConfig& cfg ){
    static constexpr const char* kBoldYellow = "\033[1;33m";
    static constexpr const char* kReset = "\033[0m";

    std::cout << "Config: " << cfg.configPath << "\n";
    if( cfg.usedDefaultConfig ){
        std::cout << kBoldYellow
                  << "WARNING: no CONFIG=path or L2RESIDUALS_CONFIG was provided; "
                  << "using the default cfg/2024ppRef.toml. Pass CONFIG="
                  << cfg.configPath
                  << " to make this run explicit."
                  << kReset << "\n";
    }
}
