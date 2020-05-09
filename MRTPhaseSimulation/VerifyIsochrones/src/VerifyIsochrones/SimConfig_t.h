//
// Created by konstantin on 5/2/20.
//

#ifndef NEWBYSCHWEMMER_SIMCONFIG_T_H
#define NEWBYSCHWEMMER_SIMCONFIG_T_H

#include <string>
#include <filesystem>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace fs = std::filesystem;
namespace pt = boost::property_tree;

struct SimConfig_t {
    std::string modelname;
    struct Paths {
        fs::path input;
        fs::path output;
    } Paths;
    struct Simulation {
        double dt;
        double t0;
        double T;
        unsigned int EnsembleSize;
        unsigned int SampleSize;
    } Simulation;
    void load_from_file(fs::path filepath);
};


#endif //NEWBYSCHWEMMER_SIMCONFIG_T_H
