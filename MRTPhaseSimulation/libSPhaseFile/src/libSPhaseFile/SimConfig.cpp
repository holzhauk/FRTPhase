//
// Created by konstantin on 1/16/21.
//

#include "SimConfig.h"

bool Config::operator==(const Config &other) {
    bool is_equal = true;
    is_equal = is_equal && (this->modelName == other.modelName);
    is_equal = is_equal && (this->paths.input == other.paths.input);
    is_equal = is_equal && (this->paths.output == other.paths.output);
    is_equal == is_equal && (this->simulation == other.simulation);
    return is_equal;
}

bool Config::Simulation::operator== (const Config::Simulation& other) {
    bool is_equal = true;
    is_equal = is_equal && (this->dt == other.dt);
    is_equal = is_equal && (this->t0 == other.t0);
    is_equal = is_equal && (this->T == other.T);
    is_equal = is_equal && (this->EnsembleSize == other.EnsembleSize);
    is_equal = is_equal && (this->SampleSize == other.SampleSize);
    return is_equal;
}

void SimConfigFile::read(const fs::path& filepath) {
    pt::ptree root;
    pt::read_json(filepath, root);
    config.modelName = root.get<std::string>("Model Name");

    pt::ptree inPaths = root.get_child("Paths.In");
    config.paths.input = fs::path(inPaths.get<std::string>("filepath")) /
                  fs::path(inPaths.get<std::string>("filename"));

    pt::ptree outPaths = root.get_child("Paths.Out");
    config.paths.output = fs::path(outPaths.get<std::string>("filepath")) /
                   fs::path(outPaths.get<std::string>("filename"));

    pt::ptree simulation_configuration = root.get_child("Simulation");
    config.simulation.dt = simulation_configuration.get<double>("dt");
    config.simulation.t0 = simulation_configuration.get<double>("t0");
    config.simulation.T = simulation_configuration.get<double>("T");
    config.simulation.EnsembleSize = simulation_configuration.get<double>("Ensemble Size");
    config.simulation.SampleSize = simulation_configuration.get<double>("Sample Size");
}

void SimConfigFile::write(const fs::path& filepath) {
    pt::ptree root;
    root.put("Model Name", config.modelName);

    pt::ptree inPaths;
    inPaths.put("filepath", string(config.paths.input.parent_path()));
    inPaths.put("filename", string(config.paths.input.filename()));
    root.put_child("Paths.In", inPaths);

    pt::ptree outPaths;
    outPaths.put("filepath", string(config.paths.output.parent_path()));
    outPaths.put("filename", string(config.paths.output.filename()));
    root.put_child("Paths.Out", outPaths);

    pt::ptree simulation;
    simulation.put("dt", config.simulation.dt);
    simulation.put("t0", config.simulation.t0);
    simulation.put("T", config.simulation.T);
    simulation.put("Ensemble Size", config.simulation.EnsembleSize);
    simulation.put("Sample Size", config.simulation.SampleSize);
    root.put_child("Simulation", simulation);

    pt::write_json(filepath, root);
}

Config::Simulation SimConfigFile::get_simConfig() const{
    return Config::Simulation(config.simulation);
}

fs::path SimConfigFile::get_inPath() const {
    return fs::path(config.paths.input);
}

fs::path SimConfigFile::get_outPath() const {
    return fs::path(config.paths.output);
}

string SimConfigFile::get_modelName() const {
    return string(config.modelName);
}

bool SimConfigFile::operator == (const SimConfigFile& other) {
    return this->config == other.config;
}