//
// Created by konstantin on 4/7/20.
//
#include <cmath>

#include "Domain.h"

/*
 * Domain
 */
Domain::Domain(double rho_min, double rho_max) {
    if (rho_max >= rho_min) {
        rho_UB = rho_max;
        rho_LB = rho_min;
    } else
        {
        rho_LB = rho_max;
        rho_UB = rho_min;
    }
    width = rho_UB - rho_LB;
}

bool Domain::hit_LowerBoundary(std::array<double, 2>& pos) {
    return (pos[0] < rho_LB);
}

bool Domain::hit_UpperBoundary(std::array<double, 2>& pos) {
    return (pos[0] > rho_UB);
}

std::array<double, 2> Domain::reflect_at_LowerBoundary(std::array<double, 2>& pos) {
    // pos = [rho, phi]
    double Drho = abs(pos[0] - rho_LB);
    int m = floor(Drho / width);
    if (m % 2 == 0){
        return {rho_LB + (Drho - m*width), pos[1]};
    } else{
        return {rho_UB - (Drho - m*width), pos[1]};
    }
}

std::array<double, 2> Domain::reflect_at_UpperBoundary(std::array<double, 2>& pos) {
    // pos = [rho, phi]
    double Drho = abs(pos[0] - rho_UB);
    int m = floor(Drho / width);
    if (m % 2 == 0){
        return {rho_UB - (Drho - m*width), pos[1]};
    } else{
        return {rho_LB + (Drho - m*width), pos[1]};
    }
}

std::array<double, 2> Domain::apply(std::array<double, 2>& pos) {
    if (this->hit_LowerBoundary(pos)){
        return this->reflect_at_LowerBoundary(pos);
    } else{
        return this->reflect_at_UpperBoundary(pos);
    }
}