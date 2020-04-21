//
// Created by konstantin on 4/6/20.
//
#include "Isochrone.h"
#include <iostream>

Isochrone::Isochrone(std::string& isochrone_group_name) {
    name = isochrone_group_name;
}

void Isochrone::push_parameter(std::string key, double& val) {
    parameters[key] = val;
}

void Isochrone::set_curve(std::vector<double>& Rhos, std::vector<double>& Phis) {
    Rho = std::unique_ptr<std::vector<double>>(new std::vector<double>(Rhos.size()));
    Phi = std::unique_ptr<std::vector<double>>(new std::vector<double>(Phis.size()));
    *Rho = Rhos;
    *Phi = Phis;
}

void Isochrone::print() {
    std::cout << "Parameters:" << std::endl;
    for (auto p_i : parameters){
        std::cout << "\t" << p_i.first << ": " <<
                  p_i.second << std::endl;
    }
    std::cout << "Curve:" << std::endl;
    std::cout << "\t Rho: [" ;
    for (auto rho_i: *Rho){
        std::cout << rho_i << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << "\t Phi: [" ;
    for (auto phi_i: *Phi){
        std::cout << phi_i << ", ";
    }
    std::cout << "]" << std::endl;
}

std::vector<double>& Isochrone::get_Rho() {
    return *Rho;
}

std::vector<double>& Isochrone::get_Phi() {
    return *Phi;
}

double Isochrone::get_Parameter(const std::string Pname) {
    return parameters[Pname];
}

std::string Isochrone::get_name() {
    return name;
}