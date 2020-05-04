//MPI
#include <mpi.h>

// CXX standard libraries
#include <iostream>
#include <filesystem>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <memory>
#include <array>
#include <vector>

// Tiz Master Library
#include "Isochrone.h"
#include "IsochroneSet.h"
#include "IsoPlanarSOsc.h"
#include "Domain.h"
#include "MFPTs.h"
#include "MFPTSet.h"
#include "NormMoments.h"
#include "ModelFactory.h"
#include "SimConfig_t.h"
#include "Adapters.h"

namespace fs = std::filesystem;

struct bcast_dimensions_t {
    unsigned int noIsochrones;
    unsigned int sizeOfModelname;
};

struct bcast_config_t {
    unsigned int noCoordinates = 0;
    unsigned int noSamples = 0;
    unsigned int ensembleSize = 0;
    bcast_config_t() = default;
    bcast_config_t(SimConfig_t& config, Isochrone& isochrone){
        noCoordinates = isochrone.get_Rho().size();
        noSamples = config.Simulation.SampleSize;
        ensembleSize = config.Simulation.EnsembleSize;
    }
    MPI_Datatype mpiType(){
        MPI_Datatype mpi_bcast_config_t;
        int blengths[3] = {1, 1, 1};
        MPI_Datatype types[3] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED};
        MPI_Aint offsets[3];
        offsets[0] = offsetof(bcast_config_t, noCoordinates);
        offsets[1] = offsetof(bcast_config_t, noSamples);
        offsets[2] = offsetof(bcast_config_t, ensembleSize);
        MPI_Type_create_struct(3, blengths, offsets, types, &mpi_bcast_config_t);
        MPI_Type_commit(&mpi_bcast_config_t);
        return mpi_bcast_config_t;
    }
};

void gen_GathervCounts(int data_size, int world_size,
        std::vector<int>& SendCounts, std::vector<int>& Displs);

std::unique_ptr<std::array<double, 4>> run_AsymptoticAngle(IsoPlanarSOsc& Model,
        IsoPlanarSOsc::config_t& sim_config, IsoPlanarSOsc::pSet_t& pSet,
        Domain& Disk, int EnsembleSize);

std::unique_ptr<std::array<double, 2>> run_FP(IsoPlanarSOsc& Model,
        IsoPlanarSOsc::config_t& config, IsoPlanarSOsc::pSet_t& pSet,
        Domain& Disk, int EnsembleSize, double MeanPhiT,
        std::vector<double>& IsoRho, std::vector<double>& IsoPhi);

int main(int argc, char* argv[]) {

    if (argc != 2){
        std::cerr << "Argument missing! Usage: "
            << argv[0] << " <ConfigurationFilePath>.json" << std::endl;
        return EXIT_FAILURE;
    }

    /*
     * initialize MPI
     */
    MPI_Init(nullptr, nullptr);
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Datatype mpi_bcast_dimensions_t;
    MPI_Type_contiguous(sizeof(bcast_dimensions_t), MPI_BYTE, &mpi_bcast_dimensions_t);
    MPI_Type_commit(&mpi_bcast_dimensions_t);

    ModelFactory TheFactory;

    if (world_rank == 0){

        /*
         * master process block
         */
        fs::path configPath = fs::path(argv[1]);
        SimConfig_t sim_config;
        sim_config.load_from_file(configPath);

        IsochroneSet IsoSet;
        IsoSet.load(sim_config.Paths.input);

        /*
         * Broadcast dimensions of communication types
         */
        bcast_dimensions_t bcast_dim;
        bcast_dim.noIsochrones = std::distance(IsoSet.begin(), IsoSet.end());
        bcast_dim.sizeOfModelname = sim_config.modelname.length() + 1;
        MPI_Bcast(&bcast_dim, 1, mpi_bcast_dimensions_t, 0, MPI_COMM_WORLD);

        /*
         * Broadcast model name
         */
        MPI_Bcast(sim_config.modelname.data(), (int) bcast_dim.sizeOfModelname,
                MPI_CHAR, 0, MPI_COMM_WORLD);

        auto Model_ptr = TheFactory.create_Model(sim_config.modelname);
        auto pSet_ptr = TheFactory.create_pSet(sim_config.modelname);
        IsoPlanarSOsc::config_t config;

        MFPTSet mfptSet(sim_config.modelname, configPath, sim_config.Paths.input);

        for (auto Iso_ptr : IsoSet){

            pSet_ptr->load(Iso_ptr->get_parameterMap());
            config = Adapters::Config2SimConfig(sim_config);

            auto mfpts_ptr = mfptSet.create(Iso_ptr->get_name());

            std::cout << "Isochrone with D: " << Iso_ptr->get_parameterMap()["D"] << std::endl;

            std::vector<double> IsoRho = Iso_ptr->get_Rho();
            std::vector<double> IsoPhi = Iso_ptr->get_Phi();

            int n = (int) IsoRho.size() / sim_config.Simulation.SampleSize; // 10 sample points per trajectory
            std::vector<double> IsoSample_rho = std::vector<double>(sim_config.Simulation.SampleSize);
            std::vector<double> IsoSample_phi = std::vector<double>(sim_config.Simulation.SampleSize);

            int j = 0;
            for (int i = 0; i < IsoRho.size(); i += n){
                IsoSample_rho[j] = IsoRho[i];
                IsoSample_phi[j] = IsoPhi[i];
                j++;
            }

            boundary_vals_t boundaryVals {  *std::min_element(IsoRho.begin(),IsoRho.end()),
                                            *std::max_element(IsoRho.begin(), IsoRho.end())
                                            };
            Domain Disk(boundaryVals);

            MPI_Bcast(&config, 1, config.mpiType(), 0, MPI_COMM_WORLD);
            MPI_Bcast(pSet_ptr.get(), 1, pSet_ptr->mpiType(), 0, MPI_COMM_WORLD);

            bcast_config_t bcast_config(sim_config, *Iso_ptr);
            MPI_Bcast(&bcast_config, 1, bcast_config.mpiType(), 0, MPI_COMM_WORLD);

            std::vector<int> RcvCounts, Displs;
            int SubEnsembleSize;
            gen_GathervCounts(sim_config.Simulation.EnsembleSize, world_size, RcvCounts, Displs);
            SubEnsembleSize = RcvCounts[world_rank];

            auto* SendBuffer = new double [2*(bcast_config.noCoordinates + bcast_config.noSamples)];
            for (int i = 0; i < bcast_config.noCoordinates; i++){
                SendBuffer[i] = IsoRho[i];
                SendBuffer[i + bcast_config.noCoordinates] = IsoPhi[i];
            }
            for (int i = 0; i < bcast_config.noSamples; i++ ){
                SendBuffer[i + 2*bcast_config.noCoordinates] = IsoSample_rho[i];
                SendBuffer[i + 2*bcast_config.noCoordinates + bcast_config.noSamples] = IsoSample_phi[i];
            }
            MPI_Bcast(SendBuffer, 2*(bcast_config.noCoordinates + bcast_config.noSamples),
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
            delete[] SendBuffer;

            std::vector<double> MFPT = std::vector<double>(bcast_config.noSamples, 0.0);
            std::vector<double> VarFPT = std::vector<double>(bcast_config.noSamples, 0.0);
            std::vector<double> Tbar = std::vector<double>(bcast_config.noSamples, 0.0);

            for (int s = 0; s < bcast_config.noSamples; s++){

                config.x0[0] = IsoSample_rho[s];
                config.x0[1] = 0.0;

                // AsympStatistics = [MPhiT, VarPhiT, Tbar, VarT]
                std::unique_ptr<std::array<double, 4>> AsympStatistics =
                        run_AsymptoticAngle(*Model_ptr, config, *pSet_ptr,
                                Disk, SubEnsembleSize);
                auto* TotAsympStatistics = new double [4*world_size];
                MPI_Gather(AsympStatistics->data(), 4, MPI_DOUBLE,
                        TotAsympStatistics, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                double MPhiT = 0.0;
                for (int k = 0; k < 4*world_size; k += 4){
                    MPhiT += TotAsympStatistics[k];
                    Tbar[s] += TotAsympStatistics[k + 2];
                }
                delete[] TotAsympStatistics;

                MPhiT = MPhiT / world_size;
                Tbar[s] = Tbar[s] / world_size;
                MPI_Bcast(&MPhiT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                config.x0[1] = IsoSample_phi[s];
                std::unique_ptr<std::array<double, 2>> FPTStatistics =
                        run_FP(*Model_ptr, config, *pSet_ptr, Disk,
                                SubEnsembleSize, MPhiT, IsoRho, IsoPhi);
                auto* TotFPTStatistics = new double [2*world_size];
                MPI_Gather(FPTStatistics->data(), 2, MPI_DOUBLE,
                        TotFPTStatistics, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MFPT[s] = 0.0;
                VarFPT[s] = 0.0;
                for (int k = 0; k < 2*world_size; k += 2){
                    MFPT[s] += TotFPTStatistics[k];
                    VarFPT[s] += pow(TotFPTStatistics[k + 1], 2.0);
                }
                delete[] TotFPTStatistics;

                MFPT[s] = MFPT[s] / world_size;
                VarFPT[s] = sqrt(VarFPT[s] / world_size);
            }
            mfpts_ptr->add(IsoSample_rho, IsoSample_phi, Tbar, MFPT, VarFPT);
        }
        mfptSet.write_to_file(sim_config.Paths.output);

    } else {

        /*
         * slave process block
         */
        bcast_dimensions_t bcast_dim;
        MPI_Bcast(&bcast_dim, 1, mpi_bcast_dimensions_t, 0, MPI_COMM_WORLD);

        auto* dummy_bf = new char[bcast_dim.sizeOfModelname];
        MPI_Bcast(dummy_bf, bcast_dim.sizeOfModelname, MPI_CHAR, 0, MPI_COMM_WORLD);
        std::string modelName = std::string(dummy_bf);

        auto Model_ptr = TheFactory.create_Model(modelName);
        auto pSet_ptr = TheFactory.create_pSet(modelName);
        IsoPlanarSOsc::config_t config;

        for (int i = 0; i < bcast_dim.noIsochrones; i++){

            MPI_Bcast(&config, 1, config.mpiType(), 0, MPI_COMM_WORLD);
            MPI_Bcast(pSet_ptr.get(), 1, pSet_ptr->mpiType(), 0, MPI_COMM_WORLD);

            bcast_config_t bcast_config;
            MPI_Bcast(&bcast_config, 1, bcast_config.mpiType(), 0, MPI_COMM_WORLD);

            std::vector<int> RcvCounts;
            std::vector<int> Displs;
            int SubEnsembleSize;
            gen_GathervCounts(bcast_config.ensembleSize, world_size, RcvCounts, Displs);
            SubEnsembleSize = RcvCounts[world_rank];

            std::vector<double> IsoRho, IsoPhi, IsoSample_rho, IsoSample_phi;
            IsoRho = std::vector<double>(bcast_config.noCoordinates, 0.0);
            IsoPhi = std::vector<double>(bcast_config.noCoordinates, 0.0);
            IsoSample_rho = std::vector<double>(bcast_config.noSamples, 0.0);
            IsoSample_phi = std::vector<double>(bcast_config.noSamples, 0.0);
            auto* rcBuffer = new double [2*(bcast_config.noCoordinates + bcast_config.noSamples)];
            MPI_Bcast(rcBuffer, 2*(bcast_config.noCoordinates + bcast_config.noSamples),
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
            for (int k = 0; k < bcast_config.noCoordinates; k++){
                IsoRho[k] = rcBuffer[k];
                IsoPhi[k] = rcBuffer[bcast_config.noCoordinates + k];
            }
            for (int k = 0; k < bcast_config.noSamples; k++){
                IsoSample_rho[k] = rcBuffer[2*bcast_config.noCoordinates + k];
                IsoSample_phi[k] = rcBuffer[2*bcast_config.noCoordinates + bcast_config.noSamples + k];
            }
            delete[] rcBuffer;

            boundary_vals_t boundaryVals {  *std::min_element(IsoRho.begin(),IsoRho.end()),
                                            *std::max_element(IsoRho.begin(), IsoRho.end())
            };
            auto Disk = Domain(boundaryVals);

            for (int s = 0; s < bcast_config.noSamples; ++s) {

                config.x0[0] = IsoSample_rho[s];
                config.x0[1] = 0.0;
                std::unique_ptr<std::array<double, 4>> AsympStatistics =
                        run_AsymptoticAngle(*Model_ptr, config,
                                *pSet_ptr, Disk, SubEnsembleSize);
                auto* dummy_buffer = new double [4*world_size];
                MPI_Gather(AsympStatistics->data(), 4, MPI_DOUBLE,
                        dummy_buffer, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                delete[] dummy_buffer;

                double MPhiT = 0.0;
                MPI_Bcast(&MPhiT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                config.x0[1] = IsoSample_phi[s];
                std::unique_ptr<std::array<double, 2>> FPTStatistics =
                        run_FP(*Model_ptr, config, *pSet_ptr,
                                Disk, SubEnsembleSize, MPhiT, IsoRho, IsoPhi);

                dummy_buffer = new double [2*world_size];
                MPI_Gather(FPTStatistics->data(), 2, MPI_DOUBLE,
                        dummy_buffer, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                delete[] dummy_buffer;

            }
        }
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}

void gen_GathervCounts(int data_size, int world_size,
                        std::vector<int>& SendCounts, std::vector<int>& Displs){

    std::div_t int_div_res = div(data_size, world_size);
    SendCounts = std::vector<int>(world_size, int_div_res.quot);
    for (int i = 0; i < int_div_res.rem; i++){
        SendCounts[world_size-i] += 1;
    }

    Displs = std::vector<int>(world_size);
    int d = 0;
    for (auto& D : Displs){
        D = d;
        d += D;
    }
}

std::unique_ptr<std::array<double, 4>> run_AsymptoticAngle(IsoPlanarSOsc& Model,
        IsoPlanarSOsc::config_t& sim_config, IsoPlanarSOsc::pSet_t& pSet,
        Domain& Disk, int EnsembleSize){
    std::vector<double> PhiT = std::vector<double>(EnsembleSize);
    std::vector<double> Tbar = std::vector<double>(EnsembleSize);
    for (int e = 0; e < EnsembleSize; e++){
        Model.configure(Disk, sim_config, pSet);
        while (Model.in_time()) {
            Model.evolve();
        }
        auto [Rho, Phi, T] = Model.get_state();
        PhiT[e] = Phi;
        Tbar[e] = (2*M_PI*T) / Phi;
    }
    auto PhiTMoments = NormMoments(PhiT);
    auto TbarMoments = NormMoments(Tbar);
    std::unique_ptr<std::array<double, 4>> Results_ptr =
            std::make_unique<std::array<double, 4>>(std::array<double, 4>{ PhiTMoments.get_mean(),
                                                                             PhiTMoments.get_variance(),
                                                                             TbarMoments.get_mean(),
                                                                             TbarMoments.get_variance() });
    return std::move(Results_ptr);
}

std::unique_ptr<std::array<double, 2>> run_FP(IsoPlanarSOsc& Model,
        IsoPlanarSOsc::config_t& sim_config, IsoPlanarSOsc::pSet_t& pSet,
        Domain& Disk, int EnsembleSize, double MeanPhiT,
        std::vector<double>& IsoRho, std::vector<double>& IsoPhi){

    std::vector<double> FPT = std::vector<double>(EnsembleSize);

    double Phi_FP;
    double rho_min = Disk.get_rho_min();
    double rho_max = Disk.get_rho_max();
    double Phi0 = sim_config.x0[1];

    for (int e = 0; e < EnsembleSize; e++) {

        double t = 0.0;
        Model.configure(Disk, sim_config, pSet);

        bool not_passed = true;
        if (MeanPhiT > Phi0 + 2 * M_PI) {
            while (not_passed && (t <= sim_config.T)) {
                auto[rho_e, phi_e, t_e] = Model.evolve();
                // rho_i equally spaced on rho->axis
                // -> we get the index as follows
                // this could also be solved more elegantly by using the STL
                int k = floor(((double)IsoRho.size() - 1) * (rho_e - rho_min) / (rho_max - rho_min));
                Phi_FP = IsoPhi[k] + 2 * M_PI;
                not_passed = (phi_e < Phi_FP);
                t = t_e;
            }
        } else {
            if (MeanPhiT < Phi0 - 2 * M_PI) {
                while (not_passed && (t <= sim_config.T)) {
                    auto[rho_e, phi_e, t_e] = Model.evolve();
                    // rho_i equally spaced on rho->axis
                    // -> we get the index as follows
                    // this could also be solved more elegantly by using the STL
                    int k = floor(((double)IsoRho.size() - 1) * (rho_e - rho_min) /
                                  (rho_max - rho_min));
                    Phi_FP = IsoPhi[k] - 2 * M_PI;
                    not_passed = (phi_e > Phi_FP);
                    t = t_e;
                }
            } else {
                t = sim_config.T;
            }
        }
        FPT[e] = t;
    }

    std::unique_ptr<std::array<double, 2>> Results =
            std::make_unique<std::array<double, 2>>();
    NormMoments FPTMoments(FPT);
    (*Results)[0] = FPTMoments.get_mean();
    (*Results)[1] = FPTMoments.get_variance();
    return std::move(Results);
}

