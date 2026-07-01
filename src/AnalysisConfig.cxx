#include "AnalysisConfig.h"

#include "toml.hpp"

#include <cstdlib>
#include <stdexcept>

std::string DefaultConfigPath(){
    const char* envConfig = std::getenv( "L2RESIDUALS_CONFIG" );
    if( envConfig && envConfig[0] ) return envConfig;

    const char* envHome = std::getenv( "L2RESIDUALS_HOME" );
    if( envHome && envHome[0] ) return std::string( envHome ) + "/cfg/2024ppRef.toml";

    return "cfg/2024ppRef.toml";
}

AnalysisConfig LoadAnalysisConfig( const std::string& path ){
    const std::string configPath = path.empty() ? DefaultConfigPath() : path;
    const auto doc = toml::parse_file( configPath );

    AnalysisConfig cfg;

    cfg.vetoMapPath  = doc["paths"]["veto_map"].value_or( std::string{} );
    cfg.jsonPath     = TString( doc["paths"]["golden_json"].value_or( std::string{} ).c_str() );

    cfg.hiTreePath   = TString( doc["trees"]["hi"].value_or( std::string{} ).c_str() );
    cfg.skimTreePath = TString( doc["trees"]["skim"].value_or( std::string{} ).c_str() );
    cfg.trigTreePath = TString( doc["trees"]["trigger"].value_or( std::string{} ).c_str() );

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
            files.push_back( v.value_or( std::string{} ) );
        cfg.jecFilesPerCone.push_back( std::move( files ) );
    }

    for( const auto& v : *doc["binning"]["ptavg_edges"].as_array() )
        cfg.ptavgEdges.push_back( v.value_or( 0.0f ) );

    cfg.minPt = doc["cuts"]["min_pt"].value_or( 0.0f );

    if( cfg.jecFilesPerCone.empty() ) throw std::runtime_error( "jec.files is empty" );
    if( cfg.coneLabels.empty() )      throw std::runtime_error( "cones.labels is empty" );
    if( cfg.jecFilesPerCone.size() != cfg.coneLabels.size() )
        throw std::runtime_error( "jec.files and cones.labels have different lengths" );
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
