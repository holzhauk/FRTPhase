//
// Created by konstantin on 1/18/21.
//
#include <iostream>
#include <memory>
#include <mpi.h>
#include <filesystem>
#include <cmath>
#include <libSDEToolbox/libSDEToolbox.h>
#include <libSPhaseFile/libSPhaseFile.h>
#include <libMPIFunctions/libMPIFunctions.h>

using namespace std;
namespace fs = std::filesystem;

int main(int argc, char* argv[]){

    if (argc != 2) {
        std::cerr << "Argument missing! Usage: "
                  << argv[0] << " <ConfigurationFilePath>.json" << std::endl;
        return EXIT_FAILURE;
    }
    fs::path config_file_path = fs::path(argv[1]);
    /*
     * initialize MPI
     */
    int world_rank, world_size;
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    bool is_master = (world_rank == 0);
    string modelName;

    SimConfigFile config;
    Config::Simulation simConfig;
    if (is_master){
        config.read(config_file_path);
        modelName = config.get_modelName();
        simConfig = config.get_simConfig();
    }
    MPI_Share(world_rank, modelName);
    MPI_Share(world_rank, simConfig);

    IsoSurfaceFile isoSurfaceFile(modelName);
    if (is_master)
        isoSurfaceFile.read(config.get_inPath());
    MPI_Share(world_rank, isoSurfaceFile);

    FRTDataFile frtDataFile(modelName);
    if (is_master){
        fs::path in_path(config.get_inPath());
        frtDataFile = FRTDataFile(modelName, config_file_path, in_path);
    }

    EquidistantSampler sampler;
    ModelFactory theModelFactory;

    for (auto isoSurface: isoSurfaceFile) {
        auto [rho_min, rho_max] = isoSurface.get_extensions();
        if (is_master)
            cout << "model: " << modelName << endl;
            cout << "rho_min: " << rho_min << "; rho_max: " << rho_max << endl;
        ReflectiveAnnulus domain(rho_min, rho_max);
        if (is_master) {
            ParameterSet pSet = isoSurface.get_parameterSet();
            cout << "Parameters: " << endl;
            cout << "D: " << pSet["D"] << endl;
            cout << "omega: " << pSet["omega"] << endl;
            cout << "gamma: " << pSet["gamma"] << endl;
            cout << "c: " << pSet["c"] << endl;
        }
        unique_ptr<IsotropicPlanarSSDE> theModelPtr =
                theModelFactory.createModel(modelName, isoSurface.get_parameterSet());
        ItoEulerIntegrator integrator(&domain, theModelPtr.get());
        MPI::FRTDetector frtDetector(world_rank, world_size);
        FRTData& frtData = frtDataFile.createDataSet(isoSurface.get_name());
        frtData = frtDetector.run(simConfig, isoSurface, &sampler, &integrator);
    }
    if (is_master)
        frtDataFile.write(config.get_outPath());
    MPI_Finalize();
    return 0;
}

