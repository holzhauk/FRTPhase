//
// Created by konstantin on 1/8/21.
//

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE libSDEToolbox test

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <libSDEToolbox/libSDEToolbox.h>
#include "../src/libSDEToolbox/SDEIntegrator.h"

using namespace std;

BOOST_AUTO_TEST_CASE(SDEIntegrator_test){

   class TestInt : public SDEIntegrator {
   public:
       TestInt(config_t& config): SDEIntegrator(config){};
       state_t evolve() {
           x[0] = x[0] + dt;
           x[1] = x[1] + dt;
           t += dt;
           return tie(x, t);
       }
   };

   TestInt::config_t config = {{0.0, 0.0}, 0.0, 1.0, 10.0};
   TestInt testInt(config);
   auto [x_1, t_1] = testInt.evolve();
   auto [x_f, t_f] = testInt.integrate();

   cout << "evolve: " << x_1[0] << "\t integrate: " << x_f[0] << endl;

}