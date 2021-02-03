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
        fs::path input = fs::path();
        fs::path output = fs::path();
    } paths;
    struct Simulation {
        double dt = 0.0;
        double t0 = 0.0;
        double T = 0.0;
        size_t EnsembleSize = 0;
        size_t SampleSize = 0;
        bool operator == (const Simulation& other);
    } simulation;
    string modelName = string();
    list<ParameterSet> pSetList;
    //Paths paths;
    //Simulation simulation;
    Config() = default;
    bool operator == (const Config& other);
};

class SimConfigFile {
public:
    class paramSetIterator {
        using pSetIt = list<ParameterSet>::iterator;
    private:
        pSetIt it;
    public:
        paramSetIterator(pSetIt listIt): it(listIt) {};
        bool operator != (const paramSetIterator& other);
        ParameterSet& operator ++ ();
        ParameterSet& operator * ();
    };
private:
    Config config = Config();
public:
    SimConfigFile() = default;
    SimConfigFile(Config& config): config(config) {};
    void read(const fs::path& filepath);
    void write(const fs::path& filepath);
    paramSetIterator pSet_begin();
    paramSetIterator pSet_end();
    Config::Simulation get_simConfig() const;
    fs::path get_inPath() const;
    fs::path get_outPath() const;
    string get_modelName() const;
    bool operator == (const SimConfigFile& other);
};


#endif //MRTPHASESIMULATION_SIMCONFIG_H
