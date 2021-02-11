//
// Created by konstantin on 1/12/21.
//

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE libMPIFunctions test_module

#include <boost/test/unit_test.hpp>
#include <mpi.h>
#include <iostream>
#include <random>
#include <libSPhaseFile/libSPhaseFile.h>
#include <libMPIFunctions/libMPIFunctions.h>

int world_rank, world_size;

struct MPIConfig {
    MPIConfig() {
        MPI_Init(nullptr, nullptr);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    };
    ~MPIConfig() {
        MPI_Finalize();
    };
};

BOOST_GLOBAL_FIXTURE(MPIConfig);

BOOST_AUTO_TEST_CASE(MPI_share_IsoSurfaceFile_test){

    string modelName = "default";
    IsoSurfaceFile isoSurfaceFile(modelName);
    if (world_rank == 0){
        // setting up test object to be transmitted
        // only the master process knows the composition rules
        string another_name = "another";
        isoSurfaceFile = IsoSurfaceFile(another_name);

        InterpolatedCurve& c1 = isoSurfaceFile.createInterpolatedCurve("Surface0");
        c1.add_parameter("D", 0.1);
        c1.add_parameter("omega", 1.0);
        c1.add_parameter("gamma", 15.0);
        c1.add_parameter("c", -15.0);
        c1.set_omegaBar(-5.934);
        vector<double> rhos = {0.1, 0.2, 0.5, 100.0};
        vector<double> phis = {1.5, 0.1, 0.0004, -15.0};
        c1.set_nodes(rhos, phis);

        InterpolatedCurve& c2 = isoSurfaceFile.createInterpolatedCurve("Surface1");
        c2.add_parameter("sigma", 0.15);
        c2.add_parameter("delta", 2.71);
        c2.add_parameter("gamma", 3.141592);
        c2.add_parameter("p0", 1.0);
        c2.set_omegaBar(0.564);
        rhos = {0.3, 0.8, 0.9623, 100.1};
        phis = {1.8, 3.0, 0.0004, -100.0};
        c2.set_nodes(rhos, phis);
    }
    MPI_Share(world_rank, isoSurfaceFile);

    BOOST_REQUIRE(isoSurfaceFile.get_NoSurfaces() == 2);

    // setting up reference object
    // both processes know its composition rule
    if (world_rank != 0) {
        string testString = "another";
        IsoSurfaceFile testFile(testString);

        InterpolatedCurve &c1 = testFile.createInterpolatedCurve("Surface0");
        c1.add_parameter("D", 0.1);
        c1.add_parameter("omega", 1.0);
        c1.add_parameter("gamma", 15.0);
        c1.add_parameter("c", -15.0);
        c1.set_omegaBar(-5.934);
        vector<double> rhos = {0.1, 0.2, 0.5, 100.0};
        vector<double> phis = {1.5, 0.1, 0.0004, -15.0};
        c1.set_nodes(rhos, phis);

        InterpolatedCurve &c2 = testFile.createInterpolatedCurve("Surface1");
        c2.add_parameter("sigma", 0.15);
        c2.add_parameter("delta", 2.71);
        c2.add_parameter("gamma", 3.141592);
        c2.add_parameter("p0", 1.0);
        c2.set_omegaBar(0.564);
        rhos = {0.3, 0.8, 0.9623, 100.1};
        phis = {1.8, 3.0, 0.0004, -100.0};
        c2.set_nodes(rhos, phis);

        // compare transmitted object to the reference one
        BOOST_REQUIRE(isoSurfaceFile == testFile);
    }
};

BOOST_AUTO_TEST_CASE(MPI_Share_ConfigSimulation_test) {
    bool is_master = (world_rank == 0);

    Config::Simulation simConfig;
    simConfig.t0 = 0.0;
    simConfig.dt = 0.025;
    simConfig.T = 1000.0;
    simConfig.EnsembleSize = 13;
    simConfig.SampleSize = 20;

    MPI_Share(world_rank, simConfig);

    Config::Simulation simConfig_reference;
    simConfig_reference.t0 = 0.0;
    simConfig_reference.dt = 0.025;
    simConfig_reference.T = 1000.0;
    simConfig_reference.EnsembleSize = 13;
    simConfig_reference.SampleSize = 20;
    BOOST_REQUIRE(simConfig == simConfig_reference);
    unsigned int e = SubEnsembleSize(world_rank, world_size, simConfig.EnsembleSize);
    if (is_master){
        BOOST_CHECK(e == 6);
    } else {
        BOOST_CHECK(e == 7);
    }

    unsigned int E;
    MPI_Allreduce(&e, &E, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    BOOST_REQUIRE(E == 13);
}

BOOST_AUTO_TEST_CASE(MPI_Dist_StatsFromData_test, * boost::unit_test::tolerance(0.01)){
    //random number generator
    unsigned seed = 976142536; // define seed for reproducibility
    std::mt19937 rn_engine(seed);
    std::normal_distribution<double> n_dist(150.0, 0.4);
    unsigned int eSize = 50000;

    Config::Simulation simConfig;
    simConfig.T = 100.0;

    // StationaryStatsTest
    Data<StationaryStats> sData(simConfig);
    double dummy;
    for (unsigned int e = 0; e < eSize; e++){
        dummy = n_dist(rn_engine);
        sData.add_data_point(dummy);
    }

    StationaryStats sStats = MPI_Dist_StatsFromData(sData, world_size*eSize);
    BOOST_TEST(sStats.mPhiT == 150.0);
    BOOST_REQUIRE(sStats.varPhiT >= 0.0);
    BOOST_TEST(sStats.varPhiT == pow(0.4, 2));
    BOOST_TEST(sStats.mT >= 1 / 150.0);
    BOOST_REQUIRE(sStats.varT >= 0.0);
    BOOST_TEST(sStats.mT == 4.18);
    BOOST_TEST(sStats.varT == 0.0001258);

    // FRTStatsTest
    Data<FRTStats> frtData(simConfig);
    for (unsigned int e = 0; e < eSize; e++){
        dummy = n_dist(rn_engine);
        frtData.add_data_point(dummy);
    }

    FRTStats frtStats = MPI_Dist_StatsFromData(frtData, world_size*eSize);
    BOOST_TEST(frtStats.mFRT == 150.0);
    BOOST_REQUIRE(frtStats.varFRT >= 0.0);
    BOOST_TEST(frtStats.varFRT == pow(0.4, 2));
}