//
// Created by konstantin on 1/16/21.
//

#ifndef MRTPHASESIMULATION_SIMCONFIG_H
#define MRTPHASESIMULATION_SIMCONFIG_H

#include <string>
#include <list>
#include <memory>
#include <filesystem>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "IsoSurface.h"

namespace fs = std::filesystem;
namespace pt = boost::property_tree;
using namespace std;

struct Config {
    struct Paths {
        /*
         * store as absolute system paths
         */
        fs::path input = fs::path();
        fs::path output = fs::path();
    } paths;
    struct Simulation {
        double dt = 0.0;
        double t0 = 0.0;
        double T = 0.0;
        size_t EnsembleSize = 0;
        size_t SampleSize = 0;
        bool operator == (const Simulation& other) const;
    } simulation;
    string modelName = string();
    list<ParameterSet> pSetList;
    string domainName = string();
    list<ParameterSet> domainDimList;
    //Paths paths;
    //Simulation simulation;
    Config() = default;
    bool operator == (const Config& other);
};

class SimConfigFile {
private:
    Config config = Config();
    bool paths_are_not_absolute = true;
public:
    SimConfigFile() = default;
    SimConfigFile(Config& config): config(config) {};
    /*
     * absolute and relative paths are accepted
     * however, the Config structure is supposed to hold
     * absolute system paths -> a conversion is implemented
     */
    void read(const fs::path& filepath);
    void write(const fs::path& filepath);
    void write(const fs::path& filepath, bool write_relative_paths);
    Config::Simulation get_simConfig() const;
    fs::path get_inPath() const;
    fs::path get_outPath() const;
    string get_modelName() const;
    bool operator == (const SimConfigFile& other);
};


#endif //MRTPHASESIMULATION_SIMCONFIG_H
