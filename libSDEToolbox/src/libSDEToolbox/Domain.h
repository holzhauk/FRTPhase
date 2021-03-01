//
// Created by konstantin on 1/8/21.
//

#ifndef MRTPHASESIMULATION_DOMAIN_H
#define MRTPHASESIMULATION_DOMAIN_H

#include <array>

using namespace std;

class Domain {
public:
    virtual ~Domain() {};
    virtual array<double, 2> apply_boundary_conditions(array<double, 2>& x_init, array<double, 2>& x_final) = 0;
};

#endif //MRTPHASESIMULATION_DOMAIN_H
