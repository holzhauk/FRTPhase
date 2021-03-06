//
// Created by konstantin on 1/8/21.
//

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE libSDEToolbox test

#include <iostream>
#include <cmath>
#include <array>
#include <memory>
#include <boost/test/unit_test.hpp>
#include <libSDEToolbox/libSDEToolbox.h>

using namespace std;

/*
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
*/

BOOST_AUTO_TEST_CASE(ReflectiveAnnulus_test){
    double rho_min, rho_max;
    rho_min = 0.5;
    rho_max = 1.5;
    ReflectiveAnnulus domain(rho_min, rho_max);
    Pos_t x0 = {0.6, 0.0};
    Pos_t x = {0.8, 0.02};
    Pos_t x_p = domain.apply_boundary_conditions(x0, x);
    BOOST_CHECK(x_p == x);
    x = {0.4, 0.02};
    x_p = domain.apply_boundary_conditions(x0, x);
    Pos_t x_c = {0.6, 0.02};
    BOOST_CHECK(x_p == x_c);
    x = {1.7, 0.6};
    x_p = domain.apply_boundary_conditions(x0, x);
    x_c = {1.3, 0.6};
    BOOST_CHECK(x_p == x_c);
    x = {3.8, 0.6};
    x_p = domain.apply_boundary_conditions(x0, x);
    BOOST_CHECK(x_p[0] - 1.2 < 0.00001);
    BOOST_CHECK_EQUAL(x_p[1], 0.6);
    //std::cout << "x = {" << x_p[0] << ", " << x_p[1] << "}" << std::endl;
    x = {-3.8, 0.6};
    x_p = domain.apply_boundary_conditions(x0, x);
    BOOST_CHECK(x_p[0] - 0.8 < 0.00001);
    BOOST_CHECK_EQUAL(x_p[1], 0.6);
}

BOOST_AUTO_TEST_CASE(ModelFactory_test, * boost::unit_test::tolerance(0.01)){
    ParameterSet NewbyS_pSet = { {"D", 0.5},
                                 {"omega", 1.5},
                                 {"gamma", 15.0},
                                 {"c", -15.0} };
    ParameterSet SchwabedalP_pSet = { {"sigma", 0.2},
                                      {"c", 1.8},
                                      {"omega", 1.0},
                                      {"delta", 1.5} };

    double rho = 0.278;
    ModelFactory theFactory;
    unique_ptr<IsotropicPlanarSSDE> theModel;
    string name;

    try {
        name = string("NewbySchwemmer");
        theModel = theFactory.createModel(name, NewbyS_pSet);

        // test if the selected model is indeed the NewbySchwemmer one
        BOOST_TEST(theModel->g(rho) == 5.64628);
        BOOST_TEST(theModel->f(rho) == -174.43335);
        BOOST_TEST(theModel->q_rho(rho) == 1.0);
        BOOST_TEST(theModel->q_phi(rho) == 3.5971);
    }
    catch (ModelNotDefined) {
        BOOST_CHECK(false);
    }

    try {
        name = string("SchwabedalPikovsky");
        theModel = theFactory.createModel(name, SchwabedalP_pSet);

        BOOST_TEST(theModel->g(rho) == 0.831543);
        BOOST_TEST(theModel->f(rho) == -3.548283);
        BOOST_TEST(theModel->q_rho(rho) == 0.05560);
        BOOST_TEST(theModel->q_phi(rho) == 0.0);
    }
    catch (ModelNotDefined) {
        BOOST_CHECK(false);
    }

    // check if the exception is catched
    bool is_catched = false;
    try {
        name = string("NotInZoo");
        theModel = theFactory.createModel(name, SchwabedalP_pSet);
    }
    catch (ModelNotDefined) {
        is_catched = true;
    }
    BOOST_TEST(is_catched);
}

