//
// Created by konstantin on 2/22/21.
//

#include "ModelFactory.h"

unique_ptr<IsotropicPlanarSSDE> ModelFactory::createModel(const string& modelName, ParameterSet pSet) {
    unique_ptr<IsotropicPlanarSSDE> modelPtr;
    if (modelName == "NewbySchwemmer") {
        modelPtr = unique_ptr<IsotropicPlanarSSDE> (new ModelZoo::NewbySchwemmer(pSet));
        return move(modelPtr);
    }
    if (modelName == "SchwabedalPikovsky") {
        modelPtr = unique_ptr<IsotropicPlanarSSDE> (new ModelZoo::SchwabedalPikovsky(pSet));
        return move(modelPtr);
    }
    if (modelName == "SimpleModel") {
        modelPtr = unique_ptr<IsotropicPlanarSSDE> (new ModelZoo::SimpleModel(pSet));
        return move(modelPtr);
    }
    throw ModelNotDefined();
}
