//
// Created by konstantin on 4/9/20.
//

#include "Trajectory.h"

Trajectory::Trajectory() {
    Rho_ptr = new std::vector<double>;
    Phi_ptr = new std::vector<double>;
    t_ptr = new std::vector<double>;
}

Trajectory::~Trajectory() {
    if (Rho_ptr != NULL)
        delete Rho_ptr;
    if (Phi_ptr != NULL)
        delete Phi_ptr;
    if (t_ptr != NULL)
        delete t_ptr;
}

void Trajectory::push_back(double& rho, double& phi, double& t) {
    Rho_ptr->push_back(rho);
    Phi_ptr->push_back(phi);
    t_ptr->push_back(t);
}

void Trajectory::push_back(std::array<double, 3>& state) {
    // state = [rho, phi, t]
    this->push_back(state[0], state[1], state[2]);
}

std::vector<double>& Trajectory::get_Rhos() {
    return *Rho_ptr;
}

double* Trajectory::get_Rho_buf_ptr() {
    return Rho_ptr->data();
}

std::vector<double>& Trajectory::get_Phis() {
    return *Phi_ptr;
}

double* Trajectory::get_Phi_buf_ptr() {
    return Phi_ptr->data();
}

std::vector<double>& Trajectory::get_ts() {
    return *t_ptr;
}

double* Trajectory::get_t_buf_ptr() {
    return t_ptr->data();
}