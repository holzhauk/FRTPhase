// CXX standard libraries
#include <iostream>
#include <filesystem>
#include <cmath>

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
#define SAMPLE_SIZE 10
#define TIME_STEP 0.01
#define T_INF 1000.0

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

    IsochroneSet IsoSet;
    std::string modelname = MODELNAME;
    fs::path IsoFilePath = iso_path / iso_model / iso_file;
    std::cout << "IsochroneDataFile: " << IsoFilePath << std::endl;
    IsoSet.load(IsoFilePath);

    MFPTSet mfptSet(modelname, IsoFilePath);
    NewSOsc NewS_AntiR;

    for (auto Iso_ptr: IsoSet) {
        MFPTs *mfpts_ptr = new MFPTs(Iso_ptr->get_name(), ENSEMBLE_SIZE, T_INF);

        pset_NewS_t pset_NewS_antirotating = {Iso_ptr->get_Parameter("D"),
                                              Iso_ptr->get_Parameter("omega"),
                                              Iso_ptr->get_Parameter("gamma"),
                                              Iso_ptr->get_Parameter("c")
        };

        std::cout << "Isochrone with D: " << Iso_ptr->get_Parameter("D") << std::endl;

        std::vector<double> &IsoRho = Iso_ptr->get_Rho();
        std::vector<double> &IsoPhi = Iso_ptr->get_Phi();

        // normalize the Phi-coordinates for positive values
        /*
        double min = *std::min_element(IsoPhi.begin(), IsoPhi.end());
        int m = abs(floor(min / (2 * M_PI)));
        for (auto &x: IsoPhi) {
            x += m * 2 * M_PI;
        }
         */

        int n = IsoRho.size() / SAMPLE_SIZE; // 10 sample points per trajectory
        for (int i = 0; i < IsoRho.size(); i += n) {

            double Rho0 = IsoRho[i];
            double Phi0 = IsoPhi[i];

            std::array<double, 3> x0 = {Rho0,
                                        Phi0,
                                        0.0
            };
            sim_config_t sim_config = {TIME_STEP,
                                       T_INF,
                                       x0
            };
            double rho_min = *std::min_element(IsoRho.begin(), IsoRho.end());
            double rho_max = *std::max_element(IsoRho.begin(), IsoRho.end());
            Domain disk(rho_min, rho_max);

            std::vector<double> FPTs;
            std::vector<double> Tbars;

            for (int e = 0; e < ENSEMBLE_SIZE; e++) {

                /*
                 * mean periods
                 */
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
                if (Phi_T > Phi0 + 2*M_PI){
                    while(not_passed && (t <= T_INF)){
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
    }
    //std::cout << "everything fine up to here" << std::endl;
    mfptSet.write_to_file(mfpt_path / mfpt_file);
    //std::cout << "still fine" << std::endl;

    return EXIT_SUCCESS;
}
