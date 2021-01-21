//
// Created by konstantin on 1/8/21.
//

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE libSDEToolbox test

#include <iostream>
#include <cmath>
#include <memory>
#include <array>
#include <boost/test/unit_test.hpp>
#include <libSDEToolbox/libSDEToolbox.h>

using namespace std;

BOOST_AUTO_TEST_CASE(SDEIntegrator_test){

   ParameterSet antirotatingPSet = {{"D", 0.2},
                                    {"omega", 1.0},
                                    {"gamma", 15.0},
                                    {"c", -15.0}};
   class NewbySchwemmer: public IsotropicPlanarSSDEAdditiveNoise {
   public:
       NewbySchwemmer(ParameterSet& pSet): IsotropicPlanarSSDEAdditiveNoise(pSet["D"], pSet["D"]) {
           this->pSet = pSet;
       };
       double g(double& rho) override {
           return -pSet["gamma"]*rho*(pow(rho, 2) - 1.0) + pSet["D"] / rho;
       };
       double f(double& rho) override {
           return pSet["omega"]*(1.0 + pSet["gamma"]*pSet["c"]*pow((rho - 1.0), 2.0));
       };
   };

    class NoBDomain: public Domain{
    public:
        array<double, 2> apply_boundary_conditions(array<double, 2>& x_init, array<double, 2>& x_final){
            return array<double, 2>(x_final);
        };
    };

   double rho_min = 0.5;
   double rho_max = 1.5;
   SDEIntegrator::config_t config;
   config.x0 = {0.6, 0.0};
   config.t0 = 0.0;
   config.dt = 0.001;
   config.T = 0.1;
   NewbySchwemmer model(antirotatingPSet);
   //unique_ptr<Domain> domain_ptr(new NoBDomain());
   ReflectiveAnnulus domain(rho_min, rho_max);
   ItoEulerIntegrator integrator(config, &domain, &model);
   for (int i = 0; i < 10; i++){
       integrator.configure(config);
       auto [x, t] = integrator.integrate();
       auto [rho, phi] = x;
       cout << i << "th realization: x = (" << rho << ", " << phi << ")" << endl;
   }

}