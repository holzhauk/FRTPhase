//
// Created by konstantin on 1/18/21.
//
#include <iostream>
#include <mpi.h>
#include <filesystem>
#include <cmath>
#include <memory>
#include <libSDEToolbox/libSDEToolbox.h>
#include <libSPhaseFile/libSPhaseFile.h>
#include <libMPIFunctions/libMPIFunctions.h>

const string modelName = "NewbySchwemmer";

class NewbySchwemmer: public IsotropicPlanarSSDEAdditiveNoise{
public:
    NewbySchwemmer(ParameterSet pSet):
        IsotropicPlanarSSDEAdditiveNoise(pSet["D"], pSet["D"]) {
        this->pSet = pSet;
    };

    double g(double& rho) override{
        return -pSet["gamma"]*rho*(pow(rho, 2.0) - 1.0) + pSet["D"] / rho;
    };

    double f(double& rho) override{
        return -pSet["omega"]*(1 + pSet["gamma"]*pSet["c"]*pow((rho - 1.0), 2.0));
    };
};

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
    IsoSurfaceFile isoSurfaceFile(modelName);
    SimConfigFile config;
    Config::Simulation simConfig;
    if (is_master){
        config.read(config_file_path);
        simConfig = config.get_simConfig();
        isoSurfaceFile.read(config.get_inPath());
    }
    MPI_Share(world_rank, simConfig);
    MPI_Share(world_rank, isoSurfaceFile);

    FRTDataFile frtDataFile(modelName);
    EquidistantSampler sampler;
    for (auto isoSurface: isoSurfaceFile) {
        auto [rho_min, rho_max] = isoSurface.get_extensions();
        unique_ptr<Domain> domain_ptr(new ReflectiveAnnulus(rho_min, rho_max));
        unique_ptr<IsotropicPlanarSSDE> theModel_ptr(new NewbySchwemmer(isoSurface.get_parameterSet()));
        ItoEulerIntegrator integrator(domain_ptr, theModel_ptr);
        MPI::FRTDetector frtDetector(world_rank, world_size);
        FRTData& frtData = frtDataFile.createDataSet(isoSurface.get_name());
        frtData = frtDetector.run(simConfig, isoSurface, &sampler, &integrator);
    }
    if (is_master)
        frtDataFile.write(config.get_outPath());
    MPI_Finalize();
    return 0;
}

