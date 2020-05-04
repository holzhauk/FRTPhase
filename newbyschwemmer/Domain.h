//
// Created by konstantin on 4/7/20.
//

#ifndef NEWBYSCHWEMMER_DOMAIN_H
#define NEWBYSCHWEMMER_DOMAIN_H

#include <array>

typedef struct boundary_vals_t{
    double rho_min;
    double rho_max;
} boundary_vals_t;

class Domain {
private:
    double rho_LB = 0.5; // lower boundary
    double rho_UB = 2.0; // upper boundary
    double width = 1.5;
public:
    Domain() = default;
    Domain(double, double);
    Domain(boundary_vals_t&);

    bool hit_UpperBoundary(std::array<double, 2>&);
    bool hit_LowerBoundary(std::array<double, 2>&);
    std::array<double, 2> reflect_at_UpperBoundary(std::array<double, 2>&);
    std::array<double, 2> reflect_at_LowerBoundary(std::array<double, 2>&);
    std::array<double, 2> apply(std::array<double, 2>&);
    double get_rho_min();
    double get_rho_max();
};


#endif //NEWBYSCHWEMMER_DOMAIN_H
