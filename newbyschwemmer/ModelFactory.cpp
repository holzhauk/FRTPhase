//
// Created by konstantin on 5/2/20.
//

#include "ModelFactory.h"

std::shared_ptr<IsoPlanarSOsc> ModelFactory::create_Model(std::string name) {
    switch (modelMap[name]) {
        case NewbySchwemmer:
            return std::shared_ptr<IsoPlanarSOsc>(new NewSOsc());
        default:
            return nullptr;
    }
}

std::shared_ptr<IsoPlanarSOsc::config_t> ModelFactory::create_config(std::string name) {
    switch (modelMap[name]) {
        case NewbySchwemmer:
            return std::shared_ptr<IsoPlanarSOsc::config_t>(new NewSOsc::config_t);
        default:
            return nullptr;
    }
}


std::shared_ptr<IsoPlanarSOsc::pSet_t> ModelFactory::create_pSet(std::string name) {
    switch (modelMap[name]) {
        case NewbySchwemmer:
            return std::shared_ptr<IsoPlanarSOsc::pSet_t>(new NewSOsc::pSet_t);
        default:
            return nullptr;
    }
}

