//
// Created by konstantin on 4/2/21.
//

#include <iostream>
#include <filesystem>
#include <memory>
#include <cmath>

#include <libSDEToolbox/libSDEToolbox.h>
#include <libSPhaseFile/libSPhaseFile.h>

using namespace std;
namespace fs = std::filesystem;

int main(int argc, char* argv[]){

    if (argc != 2) {
        std::cerr << "Argument missing! Usage: "
                  << argv[0] << " <ConfigurationFilePath>.json" << std::endl;
        return EXIT_FAILURE;
    }
    fs::path config_file_path = fs::path(argv[1]);

    SimConfigFile config;
    config.read(config_file_path);
    const string modelName = config.get_modelName();
    Config::Simulation simConfig;
    simConfig = config.get_simConfig();
    IsoSurfaceFile isoSurfaceFile(modelName);
    fs::path isoSurfaceFilePath = config.get_inPath();
    isoSurfaceFile.read(isoSurfaceFilePath);

    TbarDataFile tbarDataFile(modelName,
                              config_file_path,
                              isoSurfaceFilePath);
    ModelFactory theModelFactory;

    for (auto isoSurface: isoSurfaceFile){

        TbarData& tbarData = tbarDataFile.createDataSet(isoSurface.get_name());

        auto [rho_min, rho_max] = isoSurface.get_extensions();
        ReflectiveAnnulus domain(rho_min, rho_max);

        unique_ptr<IsotropicPlanarSSDE> modelPtr =
                theModelFactory.createModel(config.get_modelName(),
                                       isoSurface.get_parameterSet());

        ItoEulerIntegrator integrator(&domain, modelPtr.get());
        SDEIntegrator::config_t integratorConfig = SimConfig2IntegratorConfig(simConfig);
        integrator.configure(integratorConfig);

        Pos_t x = isoSurface.get_random_point();
        double t = 0.0;
        integrator.reset(x, t);
        tie(x, t) = integrator.integrate();
        tbarData.Tbar = integratorConfig.T * 2 * M_PI / x[1];

    }

    tbarDataFile.write(config.get_outPath());

    return EXIT_SUCCESS;
}