//
// Created by konstantin on 4/9/20.
//

#ifndef NEWBYSCHWEMMER_TRAJECTORY_H
#define NEWBYSCHWEMMER_TRAJECTORY_H

#include <vector>
#include <array>

class Trajectory {
private:
    std::vector<double>* Rho_ptr;
    std::vector<double>* Phi_ptr;
    std::vector<double>* t_ptr;
public:
    Trajectory();
    ~Trajectory();
    void push_back(double&, double&, double&);
    void push_back(std::array<double, 3>&);
    std::vector<double>& get_Rhos();
    double* get_Rho_buf_ptr();
    std::vector<double>& get_Phis();
    double* get_Phi_buf_ptr();
    std::vector<double>& get_ts();
    double* get_t_buf_ptr();
};


#endif //NEWBYSCHWEMMER_TRAJECTORY_H
