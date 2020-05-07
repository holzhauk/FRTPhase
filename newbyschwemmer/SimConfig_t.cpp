//
// Created by konstantin on 5/2/20.
//

#include "SimConfig_t.h"

void SimConfig_t::load_from_file(fs::path filepath) {
    pt::ptree root;
    pt::read_json(filepath, root);
    modelname = root.get<std::string>("Model Name");

    pt::ptree inPaths = root.get_child("Paths.In");
    Paths.input = filepath.remove_filename() /
            fs::path(inPaths.get<std::string>("filepath")) /
            fs::path(inPaths.get<std::string>("filename"));

    pt::ptree outPaths = root.get_child("Paths.Out");
    Paths.output = filepath.remove_filename() /
            fs::path(outPaths.get<std::string>("filepath")) /
            fs::path(outPaths.get<std::string>("filename"));

    pt::ptree simulation_configuration = root.get_child("Simulation");
    Simulation.dt = simulation_configuration.get<double>("dt");
    Simulation.t0 = simulation_configuration.get<double>("t0");
    Simulation.T = simulation_configuration.get<double>("T");
    Simulation.EnsembleSize = simulation_configuration.get<double>("Ensemble Size");
    Simulation.SampleSize = simulation_configuration.get<double>("Sample Size");
}
