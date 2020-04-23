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
#include "NewSOsc.h"
#include "Domain.h"
#include "MFPTs.h"
#include "MFPTSet.h"
#include "Trajectory.h"
#include "TrajectorySet.h"
#include "NormMoments.h"

namespace fs = std::filesystem;
#define ISOPATH "../../../Isochrones"
#define MFPT_PATH "../../../SimData/MFPTs"
#define TRAJECTORIES_PATH "../../../SimData/Trajectories"
#define MODELNAME "NewbySchwemmer"
#define FILE_NAME "AntirotatingNewbySchwemmerIsochroneSet.h5"

#define ENSEMBLE_SIZE 10
#define SAMPLE_SIZE 14
#define TIME_STEP 0.01
#define T_INF 100.0

void gen_ScattervCounts(int& NoSamples, int& world_size,
        std::vector<int>& SendCounts, std::vector<int>& Displs);

std::unique_ptr<std::array<std::vector<double>, 3>> run(NewSOsc& Model, sim_config_C_t& config, pset_NewS_t& pSet,
        Domain& Disk, int EnsembleSize, std::vector<double>& Rho0, std::vector<double>& Phi0);

int main() {
    fs::path iso_path = ISOPATH;
    fs::path iso_model = MODELNAME;
    std::string model_name_tag = MODELNAME;
    fs::path iso_file = FILE_NAME;
    fs::path mfpt_path = MFPT_PATH;
    fs::path mfpt_file = "NewbySchwemmerAntirotating.h5";
    fs::path t_path = TRAJECTORIES_PATH;
    fs::path t_filename = "NewbySchwemmerAntirotating.h5";

    /*
     * initialize MPI
     */
    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    /*
     * replicate data structures for MPI inter process communication
     */
    // struct pSet_NewS_t = { double D; double omega; double gamma; double c };
    MPI_Datatype mpi_pSet_NewS_t;
    int pSet_NewS_blengths[4] = {1, 1, 1, 1};
    MPI_Datatype pSet_NewS_types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint pSet_News_offsets[4];
    pSet_News_offsets[0] = offsetof(pset_NewS_t, D);
    pSet_News_offsets[1] = offsetof(pset_NewS_t, omega);
    pSet_News_offsets[2] = offsetof(pset_NewS_t, gamma);
    pSet_News_offsets[3] = offsetof(pset_NewS_t, c);
    MPI_Type_create_struct(4, pSet_NewS_blengths,
            pSet_News_offsets, pSet_NewS_types, &mpi_pSet_NewS_t);
    MPI_Type_commit(&mpi_pSet_NewS_t);

    // struct sim_config_C_t = { double dt; double T; double x0[3]; };
    MPI_Datatype mpi_sim_config_C_t;
    int sim_config_blengths[3] = {1, 1, 3};
    MPI_Datatype sim_config_types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint sim_config_offsets[3];
    sim_config_offsets[0] = offsetof(sim_config_C_t, dt);
    sim_config_offsets[1] = offsetof(sim_config_C_t, T);
    sim_config_offsets[2] = offsetof(sim_config_C_t, x0);
    MPI_Type_create_struct(3, sim_config_blengths,
            sim_config_offsets, sim_config_types, &mpi_sim_config_C_t);
    MPI_Type_commit(&mpi_sim_config_C_t);

    if (world_rank == 0){

        /*
         * master process block
         */
        int NoIsochrones = 0;
        pset_NewS_t pSet_NewS;
        sim_config_C_t sim_config;
        double rho_boundary[2]; // [rho_min, rho_max]

        NewSOsc Model;

        IsochroneSet IsoSet;
        fs::path IsoFilePath = iso_path / iso_model / iso_file;
        IsoSet.load(IsoFilePath);
        NoIsochrones = IsoSet.get_NoIsochrones();
        MPI_Bcast(&NoIsochrones, 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::string modelname = MODELNAME;
        MFPTSet mfptSet(modelname, IsoFilePath);
        MFPTs* mfpts_ptr;

        for (Isochrone* Iso_ptr : IsoSet){
            mfpts_ptr = new MFPTs(Iso_ptr->get_name(), ENSEMBLE_SIZE, T_INF);

            std::cout << "Isochrone with D: " << Iso_ptr->get_Parameter("D") << std::endl;

            std::vector<double>* IsoRho = new std::vector<double>;
            std::vector<double>* IsoPhi = new std::vector<double>;
            *IsoRho = Iso_ptr->get_Rho();
            *IsoPhi = Iso_ptr->get_Phi();
            rho_boundary[0] = *std::min_element(IsoRho->begin(), IsoRho->end());
            rho_boundary[1] = *std::max_element(IsoRho->begin(), IsoRho->end());
            int n = IsoRho->size() / SAMPLE_SIZE; // 10 sample points per trajectory
            std::vector<double> IsoSample_rho = std::vector<double>(n);
            std::vector<double> IsoSample_phi = std::vector<double>(n);
            int j = 0;
            for (int i = 0; i < IsoRho->size(); i += n){
                IsoSample_rho[j] = (*IsoRho)[i];
                IsoSample_phi[j] = (*IsoPhi)[i];
                j++;
            }
            delete IsoRho;
            delete IsoPhi;

            MPI_Bcast(&rho_boundary, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            Domain Disk(rho_boundary[0], rho_boundary[1]);

            pSet_NewS = pset_NewS_t {Iso_ptr->get_Parameter("D"),
                                     Iso_ptr->get_Parameter("omega"),
                                     Iso_ptr->get_Parameter("gamma"),
                                     Iso_ptr->get_Parameter("c") };
            MPI_Bcast(&pSet_NewS, 1, mpi_pSet_NewS_t, 0, MPI_COMM_WORLD);

            sim_config = sim_config_C_t { TIME_STEP,
                                        T_INF,
                                          {1.0, 0.0, 0.0} };
            MPI_Bcast(&sim_config, 1, mpi_sim_config_C_t, 0, MPI_COMM_WORLD);

            int sample_size = IsoSample_rho.size();
            MPI_Bcast(&sample_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

            std::vector<int> SendCounts;
            std::vector<int> Displs;
            gen_ScattervCounts(sample_size, world_size, SendCounts, Displs);

            std::vector<double> IsoSample_rho_rcbf = std::vector<double>(SendCounts[world_rank]);
            MPI_Scatterv(IsoSample_rho.data(), SendCounts.data(), Displs.data(), MPI_DOUBLE,
                    IsoSample_rho_rcbf.data(), SendCounts[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
            std::vector<double> IsoSample_phi_rcbf = std::vector<double>(SendCounts[world_rank]);
            MPI_Scatterv(IsoSample_phi.data(), SendCounts.data(), Displs.data(), MPI_DOUBLE,
                    IsoSample_phi_rcbf.data(), SendCounts[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

            std::unique_ptr<std::array<std::vector<double>, 3>> Results =
                    run(Model, sim_config, pSet_NewS, Disk, ENSEMBLE_SIZE,
                            IsoSample_rho_rcbf, IsoSample_phi_rcbf);

            std::vector<double> MFPT_TotalSet = std::vector<double>(sample_size);
            MPI_Gatherv((*Results)[0].data(), SendCounts[world_rank], MPI_DOUBLE,
                    MFPT_TotalSet.data(), SendCounts.data(), Displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            std::vector<double> VarMFPT_TotalSet = std::vector<double>(sample_size);
            MPI_Gatherv((*Results)[1].data(), SendCounts[world_rank], MPI_DOUBLE,
                        VarMFPT_TotalSet.data(), SendCounts.data(), Displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            std::vector<double> Tbar_TotalSet = std::vector<double>(sample_size);
            MPI_Gatherv((*Results)[2].data(), SendCounts[world_rank], MPI_DOUBLE,
                        Tbar_TotalSet.data(), SendCounts.data(), Displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

            mfpts_ptr = new MFPTs(Iso_ptr->get_name(), ENSEMBLE_SIZE, T_INF);
            mfpts_ptr->add(IsoSample_rho,
                    IsoSample_phi,
                    Tbar_TotalSet,
                    MFPT_TotalSet,
                    VarMFPT_TotalSet);
            mfptSet.push_back(mfpts_ptr);
        }

        mfptSet.write_to_file(mfpt_path / mfpt_file);

    } else {

        /*
         * slave process block
         */
        int NoIsochrones = 0;
        pset_NewS_t pSet_NewS;
        sim_config_C_t sim_config;
        double rho_boundary[2]; // [rho_min, rho_max]

        NewSOsc Model;

        MPI_Bcast(&NoIsochrones, 1, MPI_INT, 0, MPI_COMM_WORLD);

        for (int i = 0; i < NoIsochrones; i++){

            MPI_Bcast(&rho_boundary, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            Domain Disk(rho_boundary[0], rho_boundary[1]);

            MPI_Bcast(&pSet_NewS, 1, mpi_pSet_NewS_t, 0, MPI_COMM_WORLD);
            MPI_Bcast(&sim_config, 1, mpi_sim_config_C_t, 0, MPI_COMM_WORLD);

            int sample_size;
            MPI_Bcast(&sample_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

            std::vector<int> SendCounts;
            std::vector<int> Displs;
            gen_ScattervCounts(sample_size, world_size, SendCounts, Displs);

            std::vector<double> DummySendBf = std::vector<double>(sample_size);
            std::vector<double> IsoSample_rho_rcbf = std::vector<double>(SendCounts[world_rank]);
            MPI_Scatterv(DummySendBf.data(), SendCounts.data(), Displs.data(), MPI_DOUBLE,
                    IsoSample_rho_rcbf.data(), SendCounts[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
            std::vector<double> IsoSample_phi_rcbf = std::vector<double>(SendCounts[world_rank]);
            MPI_Scatterv(DummySendBf.data(), SendCounts.data(), Displs.data(), MPI_DOUBLE,
                    IsoSample_phi_rcbf.data(), SendCounts[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

            std::unique_ptr<std::array<std::vector<double>, 3>> Results =
                    run(Model, sim_config, pSet_NewS, Disk, ENSEMBLE_SIZE,
                        IsoSample_rho_rcbf, IsoSample_phi_rcbf);

            std::vector<double> DummyRcBf = std::vector<double>(sample_size);
            MPI_Gatherv((*Results)[0].data(), SendCounts[world_rank], MPI_DOUBLE,
                        DummyRcBf.data(), SendCounts.data(), Displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gatherv((*Results)[1].data(), SendCounts[world_rank], MPI_DOUBLE,
                        DummyRcBf.data(), SendCounts.data(), Displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gatherv((*Results)[2].data(), SendCounts[world_rank], MPI_DOUBLE,
                        DummyRcBf.data(), SendCounts.data(), Displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

    }


    MPI_Finalize();




    /*
    std::map<std::string, double> parameters;
    parameters["D"] = 0.198;
    parameters["omega"] = 1.0;
    parameters["gamma"] = 15.0;
    parameters["c"] = -15.0;
    pset_NewS_t pset_NewbySchwemmer_antirotating = {
            parameters["D"],
            parameters["omega"],
            parameters["gamma"],
            parameters["c"]
    };
    Domain disk(0.3, 1.6);
    std::array<double, 3> initial_state = {1.5, M_PI / 2, 0.0};
    sim_config_t config = {
            0.01,
            1.0,
            initial_state
    };
    TrajectorySet TSet(parameters, model_name_tag);
    Trajectory* T_ptr = new Trajectory;
    T_ptr->push_back(initial_state);
    NewSOsc NewbySchwemmerAntirotatingModel(disk, pset_NewbySchwemmer_antirotating, config);
    */
    /*
    while(NewbySchwemmerAntirotatingModel.in_time()){
        auto [rho, phi, t] = NewbySchwemmerAntirotatingModel.evolve();
        std::cout << "[Rho, Phi](" << t << ") = [ " << rho << ", " << phi << " ]" << std::endl;
        T_ptr->push_back(rho, phi, t);
    }
    */
    /*
    for (auto& Model_ref : NewbySchwemmerAntirotatingModel){
        auto [rho, phi, t] = Model_ref.get_state();
        std::cout << "[Rho, Phi](" << t << ") = [ " << rho << ", " << phi << " ]" << std::endl;
        T_ptr->push_back(rho, phi, t);
    }

    TSet.push_back(T_ptr);
    TSet.write_to_file(t_path / t_filename);
    */

    /*
    IsochroneSet IsoSet;
    NewSOsc NewS_AntiR;
    std::string modelname = MODELNAME;
    fs::path IsoFilePath = iso_path / iso_model / iso_file;
    MFPTSet mfptSet(modelname, IsoFilePath);

    MPI_Init(NULL, NULL);
    int world_rank, world_size;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // # processes
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // process id
    //std::cout << "World_rank: " << world_rank << std::endl;

    if (world_rank == 0) {
        //std::cout << "IsochroneDataFile: " << IsoFilePath << std::endl;
        IsoSet.load(IsoFilePath);
    }

    Isochrone* Iso_ptr;
    //std::list<Isochrone*>::iterator it;
    //if (world_rank == 0){
    //    it = IsoSet.begin();
    //}

    for(int k = 0; k < 3; k++) {
        std::cout << world_rank;
        MFPTs *mfpts_ptr;
        pset_NewS_t pset_NewS_antirotating;
        std::vector<int> NoS;
        std::vector<int> Dspl;
        std::vector<double> IsoSample_rho;
        std::vector<double> IsoSample_phi;
        std::vector<double> IsoRhoPBf;

        std::div_t int_div_res = div(IsoSample_rho.size(), world_size);
        NoS = std::vector<int>(world_size, int_div_res.quot);
        for (int i = 0; i < int_div_res.rem; i++){
            NoS[world_size-i] += 1;
        }
        Dspl = std::vector<int>(world_size);
        int d = 0;
        for (int i = 0; i < Dspl.size(); i++){
            Dspl[i] = d;
            d += NoS[i];
        }

        if(world_rank==0) {
            Iso_ptr = *IsoSet.begin();
            mfpts_ptr = new MFPTs(Iso_ptr->get_name(), ENSEMBLE_SIZE, T_INF);

            pset_NewS_antirotating = pset_NewS_t {Iso_ptr->get_Parameter("D"),
                                                  Iso_ptr->get_Parameter("omega"),
                                                  Iso_ptr->get_Parameter("gamma"),
                                                  Iso_ptr->get_Parameter("c") };

            //std::cout << "Isochrone with D: " << Iso_ptr->get_Parameter("D") << std::endl;

            std::vector<double> &IsoRho = Iso_ptr->get_Rho();
            std::vector<double> &IsoPhi = Iso_ptr->get_Phi();

            int n = IsoRho.size() / SAMPLE_SIZE; // 10 sample points per trajectory
            IsoSample_rho = std::vector<double>(n);
            IsoSample_phi = std::vector<double>(n);
            int j = 0;
            for (int i = 0; i < IsoRho.size(); i += n){
                IsoSample_rho[j] = IsoRho[i];
                IsoSample_phi[j] = IsoRho[j];
                j++;
            }
        }
        //std::cout << "OpenMPI World of size: " << world_size << ", rank: " << world_rank << std::endl;

        std::cout << " " << NoS[world_rank] << " " << Dspl[world_rank] << std::endl;
        //IsoRhoPBf = std::vector<double>(NoS[world_rank]);
        //int NoS_rcCount;
        //MPI_Scatterv(IsoSample_rho.data(), NoS.data(), Dspl.data(), MPI_DOUBLE,
        //        IsoRhoPBf.data(), NoS[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //std::cout << "Rank: " << world_rank << " NoS: " << NoS_rcCount << std::endl;

    }

    MPI_Finalize();
*/
        /*
        for (int i = 0; i < IsoRho.size(); i += n) {

            double Rho0 = IsoRho[i];
            double Phi0 = IsoPhi[i];

            std::array<double, 3> x0 = { Rho0,
                                        Phi0,
                                        0.0 };
            sim_config_t sim_config = { TIME_STEP,
                                       T_INF,
                                       x0 };
            double rho_min = *std::min_element(IsoRho.begin(), IsoRho.end());
            double rho_max = *std::max_element(IsoRho.begin(), IsoRho.end());
            Domain disk(rho_min, rho_max);

            std::vector<double> FPTs;
            std::vector<double> Tbars;

            for (int e = 0; e < ENSEMBLE_SIZE; e++) {

                /*
                 * mean periods
                 */ /*
                NewS_AntiR.configure(disk, pset_NewS_antirotating, sim_config);
                while (NewS_AntiR.in_time()) {
                    NewS_AntiR.evolve();
                }
                auto[Rho_T, Phi_T, T] = NewS_AntiR.get_state();
                Tbars.push_back(T / (Phi_T / 2 * M_PI));

                double Phi_FP;
                double t = 0.0;

                NewS_AntiR.configure(disk, pset_NewS_antirotating, sim_config);

                bool not_passed = true;
                if (Phi_T > Phi0 + 2 * M_PI) {
                    while (not_passed && (t <= T_INF)) {
                        auto[rho_e, phi_e, t_e] = NewS_AntiR.evolve();
                        // rho_i equally spaced on rho->axis
                        // -> we get the index as follows
                        // this could also be solved more elegantly by using the STL
                        int k = floor((IsoRho.size() - 1) * (rho_e - rho_min) / (rho_max - rho_min));
                        Phi_FP = IsoPhi[k] + 2 * M_PI;
                        not_passed = (phi_e < Phi_FP);
                        t = t_e;
                    }
                } else {
                    if (Phi_T < Phi0 - 2 * M_PI) {
                        while (not_passed && (t <= T_INF)) {
                            auto[rho_e, phi_e, t_e] = NewS_AntiR.evolve();
                            // rho_i equally spaced on rho->axis
                            // -> we get the index as follows
                            // this could also be solved more elegantly by using the STL
                            int k = floor((IsoRho.size() - 1) * (rho_e - rho_min) / (rho_max - rho_min));
                            Phi_FP = IsoPhi[k] - 2 * M_PI;
                            not_passed = (phi_e > Phi_FP);
                            t = t_e;
                        }
                    } else {
                        t = T_INF;
                    }
                }
                //std::cout << "\t Rho0: " << Rho0 << " | FPT: " << t << std::endl;
                FPTs.push_back(t);
            }
            NormMoments FPTs_moments(FPTs);
            NormMoments Tbar_moments(Tbars);
            mfpts_ptr->add(Rho0, Phi0, Tbar_moments.get_mean(),
                           FPTs_moments.get_mean(), FPTs_moments.get_variance());
            FPTs.clear();
            Tbars.clear();
        }
        mfptSet.push_back(mfpts_ptr);

        MPI_Finalize();
    }
    //std::cout << "everything fine up to here" << std::endl;
    mfptSet.write_to_file(mfpt_path / mfpt_file);
    //std::cout << "still fine" << std::endl;

    */

    return EXIT_SUCCESS;
}

void gen_ScattervCounts(int& data_size, int& world_size,
                        std::vector<int>& SendCounts, std::vector<int>& Displs){

    std::div_t int_div_res = div(data_size, world_size);
    SendCounts = std::vector<int>(world_size, int_div_res.quot);
    for (int i = 0; i < int_div_res.rem; i++){
        SendCounts[world_size-i] += 1;
    }

    Displs = std::vector<int>(world_size);
    int d = 0;
    for (int i = 0; i < Displs.size(); i++){
        Displs[i] = d;
        d += Displs[i];
    }
}

std::unique_ptr<std::array<std::vector<double>, 3>> run(NewSOsc& Model, sim_config_C_t& sim_config, pset_NewS_t& pSet,
        Domain& Disk, int EnsembleSize, std::vector<double>& Rho0, std::vector<double>& Phi0){

    std::unique_ptr<std::array<std::vector<double>, 3>> Results =
            std::unique_ptr<std::array<std::vector<double>, 3>>(new std::array<std::vector<double>, 3>);
    for (std::vector<double>& V : *Results){
        V = std::vector<double>(Rho0.size());
    }

    for (int i = 0; i < Rho0.size(); i++){

        sim_config.x0[0] = Rho0[i];
        sim_config.x0[1] = Phi0[i];
        sim_config.x0[2] = 0.0; // t0

        std::vector<double> Tbar = std::vector<double>(EnsembleSize);
        std::vector<double> Phi_T = std::vector<double>(EnsembleSize);
        for (int e = 0; e < EnsembleSize; e++){
            Model.configure(Disk, pSet, sim_config);
            while (Model.in_time()) {
                Model.evolve();
            }
            auto [Rho, Phi, T] = Model.get_state();
            Phi_T[e] = Phi;
            Tbar[e] = (T / (Phi / 2 * M_PI));
        }
        NormMoments TbarMoments(Tbar);
        NormMoments PhiTMoments(Phi_T);
        (*Results)[2][i] = TbarMoments.get_mean();

        std::vector<double> FPT = std::vector<double>(EnsembleSize);
        std::vector<double> VarFPT = std::vector<double>(EnsembleSize);

        for (int e = 0; e < EnsembleSize; e++) {

            double Phi_FP;
            double t = 0.0;
            double rho_min = Disk.get_rho_min();
            double rho_max = Disk.get_rho_max();
            Model.configure(Disk, pSet, sim_config);

            bool not_passed = true;
            if (PhiTMoments.get_mean() > Phi0[i] + 2 * M_PI) {
                while (not_passed && (t <= T_INF)) {
                    auto[rho_e, phi_e, t_e] = Model.evolve();
                    // rho_i equally spaced on rho->axis
                    // -> we get the index as follows
                    // this could also be solved more elegantly by using the STL
                    int k = floor((Rho0.size() - 1) * (rho_e - rho_min) / (rho_max - rho_min));
                    Phi_FP = Phi0[k] + 2 * M_PI;
                    not_passed = (phi_e < Phi_FP);
                    t = t_e;
                }
            } else {
                if (PhiTMoments.get_mean() < Phi0[i] - 2 * M_PI) {
                    while (not_passed && (t <= T_INF)) {
                        auto[rho_e, phi_e, t_e] = Model.evolve();
                        // rho_i equally spaced on rho->axis
                        // -> we get the index as follows
                        // this could also be solved more elegantly by using the STL
                        int k = floor((Rho0.size() - 1) * (rho_e - rho_min) /
                                      (rho_max - rho_min));
                        Phi_FP = Phi0[k] - 2 * M_PI;
                        not_passed = (phi_e > Phi_FP);
                        t = t_e;
                    }
                } else {
                    t = T_INF;
                }
            }
            FPT[e] = t;
        }
        NormMoments MFPT(FPT);
        (*Results)[0][i] = MFPT.get_mean();
        (*Results)[1][i] = MFPT.get_variance();
    }

    return std::move(Results);
}

