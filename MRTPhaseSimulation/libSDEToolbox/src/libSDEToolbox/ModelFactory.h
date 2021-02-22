//
// Created by konstantin on 2/22/21.
//

#ifndef MRTPHASESIMULATION_MODELFACTORY_H
#define MRTPHASESIMULATION_MODELFACTORY_H

#include <memory>
#include <string>
#include <exception>

#include "ModelZoo.h"

using namespace std;

struct ModelNotDefined : public exception {
    const char* what () const throw() {
        return "ModelFactory: the model associated with the specified name is not part of the zoo";
    }
};

class ModelFactory {
public:
    unique_ptr<IsotropicPlanarSSDE> createModel(const string& modelName, ParameterSet pSet);
};


#endif //MRTPHASESIMULATION_MODELFACTORY_H
