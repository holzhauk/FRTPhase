//
// Created by konstantin on 1/10/21.
//

#include <cmath>
#include "ReflectiveAnnulus.h"

ReflectiveAnnulus::ReflectiveAnnulus(double rho_min, double rho_max){
    if (rho_min > rho_max) {
        this->rho_min = rho_max;
        this->rho_max = rho_min;
    } else {
        this->rho_min = rho_min;
        this->rho_max = rho_max;
    }
    this->width = this->rho_max - this->rho_min;
}

array<double, 2> ReflectiveAnnulus::apply_boundary_conditions(array<double, 2>& x_init, array<double, 2>& x_final) {
    if (x_final[0] < rho_min){
        // reflect at lower boundary
        double Drho = abs(x_final[0] - rho_min);
        int m = floor(Drho / width);
        if (m % 2 == 0){
            return {rho_min + (Drho - m*width), x_final[1]};
        } else{
            return {rho_max - (Drho - m*width), x_final[1]};
        }
    } else {
        if (x_final[0] > rho_max){
            // reflect at upper boundary
            double Drho = abs(x_final[0] - rho_max);
            int m = floor(Drho / width);
            if (m % 2 == 0){
                return {rho_max - (Drho - m*width), x_final[1]};
            } else{
                return {rho_min + (Drho - m*width), x_final[1]};
            }
        }
        else{
            return x_final;
        }
    }
}