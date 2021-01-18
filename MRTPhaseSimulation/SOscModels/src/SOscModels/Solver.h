//
// Created by konstantin on 5/16/20.
//

#ifndef MRTPHASESIMULATION_SOLVER_H
#define MRTPHASESIMULATION_SOLVER_H

#include <memory>
#include <array>
#include "IsoPlanarSOsc.h"

class Solver {
private:
    IsoPlanarSOsc* model_ptr;
    Domain* domain_ptr;
    IsoPlanarSOsc::config_t* config;
    std::array<double, 3> x;
public:

    struct config_t {
        double dt;
        double T;
        double x0[3];
        MPI_Datatype mpiType();
    };

    virtual std::array<double, 3> evolve() = 0;
    std::array<double, 3> run();
    bool in_time();
};


#endif //MRTPHASESIMULATION_SOLVER_H
