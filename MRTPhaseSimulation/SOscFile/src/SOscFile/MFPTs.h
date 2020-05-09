//
// Created by konstantin on 4/8/20.
//

#ifndef NEWBYSCHWEMMER_MFPTS_H
#define NEWBYSCHWEMMER_MFPTS_H

#include <string>
#include <vector>

class MFPTs {
private:
    std::string isochrone_name;
    int ensemble_size = 0;
    std::vector<double> initial_Rho;
    std::vector<double> initial_Phi;
    std::vector<double> MFPT; // mean first passage times
    std::vector<double> VarFPT; // variance in the first passage times
    double T_tot = 0.0; // total simulated time to calculate mean periods Tbar
    std::vector<double> Tbar; // mean periods
    std::vector<double> VarT;
public:
    MFPTs(std::string);

    std::string get_isochrone_name();
    int get_sample_size();
    const double* get_initial_Rhos_buf_ptr();
    const double* get_initial_Phis_buf_ptr();
    const double* get_MFPTs_buf_ptr();
    const double* get_VarFPTs_buf_ptr();
    const double* get_Tbars_buf_ptr();
    const double* get_VarT_buf_ptr();
    void add(double, double, double, double, double);
    void add(std::vector<double>&,
            std::vector<double>&,
                    std::vector<double>&,
                            std::vector<double>&,
                                    std::vector<double>&,
                                            std::vector<double>&);

};


#endif //NEWBYSCHWEMMER_MFPTS_H