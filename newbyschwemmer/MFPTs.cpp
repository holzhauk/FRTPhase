//
// Created by konstantin on 4/8/20.
//

#include "MFPTs.h"

MFPTs::MFPTs(std::string isochrone_group_name) {
    isochrone_name = isochrone_group_name;
}

std::string MFPTs::get_isochrone_name() {
    return isochrone_name;
}

int MFPTs::get_sample_size() {
    return MFPT.size();
}

const double* MFPTs::get_initial_Rhos_buf_ptr() {
    return initial_Rho.data();
}

const double* MFPTs::get_initial_Phis_buf_ptr() {
    return initial_Phi.data();
}

const double* MFPTs::get_MFPTs_buf_ptr() {
    return MFPT.data();
}

const double* MFPTs::get_VarFPTs_buf_ptr() {
    return VarFPT.data();
}

const double* MFPTs::get_Tbars_buf_ptr() {
    return Tbar.data();
}

void MFPTs::add(double rho_init, double phi_init, double mean_period, double mean_fpt, double var_fpt) {
    initial_Rho.push_back(rho_init);
    initial_Phi.push_back(phi_init);
    MFPT.push_back(mean_fpt);
    VarFPT.push_back(var_fpt);
    Tbar.push_back(mean_period);
}

void MFPTs::add(std::vector<double>& rho_init,
        std::vector<double>& phi_init,
        std::vector<double>& mean_period,
        std::vector<double>& mean_fpt,
        std::vector<double>& var_fpt){
    initial_Rho = rho_init;
    initial_Phi = phi_init;
    MFPT = mean_fpt;
    VarFPT = var_fpt;
    Tbar = mean_period;
}