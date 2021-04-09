//
// Created by konstantin on 1/16/21.
//

#include <algorithm>
#include "SimConfig.h"

bool Config::operator==(const Config &other) {
    bool is_equal = (this->modelName == other.modelName);

    auto lhsPSetIt = this->pSetList.begin();
    auto rhsPSetIt = other.pSetList.begin();
    while (lhsPSetIt != this->pSetList.end() && rhsPSetIt != other.pSetList.end()){
        is_equal = is_equal && ((lhsPSetIt->size() == rhsPSetIt->size()) &&
                (std::equal(lhsPSetIt->begin(), lhsPSetIt->end(), rhsPSetIt->begin())));
        ++lhsPSetIt;
        ++rhsPSetIt;
    }

    is_equal = is_equal && (this->domainName == other.domainName);
    auto lhsDomDimIt = this->domainDimList.begin();
    auto rhsDomDimIt = other.domainDimList.begin();
    while (lhsDomDimIt != this->domainDimList.end() && rhsDomDimIt != other.domainDimList.end()){
        is_equal = is_equal && ((lhsDomDimIt->size() == rhsDomDimIt->size()) &&
                (std::equal(lhsDomDimIt->begin(), lhsDomDimIt->end(), rhsDomDimIt->begin())));
        ++lhsDomDimIt;
        ++rhsDomDimIt;
    }

    is_equal = is_equal && (this->paths.input == other.paths.input);
    is_equal = is_equal && (this->paths.output == other.paths.output);
    is_equal = is_equal && (this->simulation == other.simulation);
    return is_equal;
}

bool Config::Simulation::operator== (const Config::Simulation& other) const {
    bool is_equal = (this->dt == other.dt);
    is_equal = is_equal && (this->t0 == other.t0);
    is_equal = is_equal && (this->T == other.T);
    is_equal = is_equal && (this->EnsembleSize == other.EnsembleSize);
    is_equal = is_equal && (this->SampleSize == other.SampleSize);
    return is_equal;
}

void SimConfigFile::read(const fs::path& filepath) {
    pt::ptree root;
    pt::read_json(filepath, root);
    pt::ptree modelG = root.get_child("Model");
    config.modelName = modelG.get<std::string>("Name");

    for (auto& ptPSet: modelG.get_child("Parameterizations")){
        ParameterSet pSet;
        for (auto& parameter: ptPSet.second){
            pSet[parameter.first.data()] = stod(parameter.second.data());
        }
        config.pSetList.push_back(pSet);
    }

    pt::ptree domainG = root.get_child("Domain");
    config.domainName = domainG.get<std::string>("Name");

    for (auto& ptPSet: domainG.get_child("Dimensions")){
        ParameterSet dimSet;
        for (auto& dim: ptPSet.second){
            dimSet[dim.first.data()] = stod(dim.second.data());
        }
        config.domainDimList.push_back(dimSet);
    }

    /*
     * implement conversions such that eventually
     * the paths stored in the config structure are absolute
     */
    pt::ptree inPaths = root.get_child("Paths.In");
    fs::path filepath_input = fs::path(inPaths.get<std::string>("filepath"));
    if (filepath_input == fs::path(""))
        filepath_input = fs::path(".");
    config.paths.input = filepath_input /
                  fs::path(inPaths.get<std::string>("filename"));

    pt::ptree outPaths = root.get_child("Paths.Out");
    fs::path filepath_output = fs::path(outPaths.get<std::string>("filepath"));
    if (filepath_output == fs::path(""))
        filepath_output = fs::path(".");
    config.paths.output = filepath_output /
                          fs::path(outPaths.get<std::string>("filename"));

    paths_are_not_absolute =
            !(config.paths.input.is_absolute() && config.paths.output.is_absolute());

    if (config.paths.input.is_relative()) {
        config.paths.input = fs::canonical(filepath).parent_path() / config.paths.input;
        config.paths.input = fs::weakly_canonical(config.paths.input);
    }

    if (config.paths.output.is_relative()) {
        config.paths.output = fs::canonical(filepath).parent_path() / config.paths.output;
        config.paths.output = fs::weakly_canonical(config.paths.output);
    }

    pt::ptree simulation_configuration = root.get_child("Simulation");
    config.simulation.dt = simulation_configuration.get<double>("dt");
    config.simulation.t0 = simulation_configuration.get<double>("t0");
    config.simulation.T = simulation_configuration.get<double>("T");
    config.simulation.EnsembleSize = simulation_configuration.get<double>("Ensemble Size");
    config.simulation.SampleSize = simulation_configuration.get<double>("Sample Size");
}

void SimConfigFile::write(const fs::path& filepath) {
    pt::ptree root;

    pt::ptree modelG;
    modelG.put("Name", config.modelName);

    pt::ptree paramG;
    for (auto& pSet: config.pSetList){
        pt::ptree ptPSet;
        for (auto& param: pSet){
            ptPSet.put(param.first.data(), param.second);
        }
        paramG.push_back(make_pair("", ptPSet));
    }
    modelG.put_child("Parameterizations", paramG);
    root.put_child("Model", modelG);

    pt::ptree domainG;
    domainG.put("Name", config.domainName);

    pt::ptree dimensionsG;
    for (auto& dimSet: config.domainDimList){
        pt::ptree ptDSet;
        for (auto& dim: dimSet){
            ptDSet.put(dim.first.data(), dim.second);
        }
        dimensionsG.push_back(make_pair("", ptDSet));
    }
    domainG.put_child("Dimensions", dimensionsG);
    root.put_child("Domain", domainG);

    pt::ptree inPaths;
    fs::path json_input_path, json_output_path;

    if (paths_are_not_absolute) {
        fs::path base = fs::weakly_canonical(filepath).parent_path();
        json_input_path = fs::relative(config.paths.input, base);
        json_output_path = fs::relative(config.paths.output, base);
    } else {
        json_input_path = config.paths.input;
        json_output_path = config.paths.output;
    }
    inPaths.put("filepath", string(json_input_path.parent_path()));
    inPaths.put("filename", string(json_input_path.filename()));
    root.put_child("Paths.In", inPaths);

    pt::ptree outPaths;
    outPaths.put("filepath", string(json_output_path.parent_path()));
    outPaths.put("filename", string(json_output_path.filename()));
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

void SimConfigFile::write(const fs::path &filepath, bool write_relative_paths) {
    this->paths_are_not_absolute = write_relative_paths;
    this->write(filepath);
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
    bool are_equal = (this->config == other.config);
    are_equal = are_equal && (this->paths_are_not_absolute && other.paths_are_not_absolute);
    return are_equal;
}