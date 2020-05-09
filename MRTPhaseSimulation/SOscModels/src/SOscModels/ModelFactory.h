//
// Created by konstantin on 5/2/20.
//

#ifndef NEWBYSCHWEMMER_MODELFACTORY_H
#define NEWBYSCHWEMMER_MODELFACTORY_H

#include <tuple>
#include <memory>
#include <map>

#include "IsoPlanarSOsc.h"
#include "NewSOsc.h"


class ModelFactory {
private:
    enum models {NewbySchwemmer};
    std::map<std::string, models> modelMap;
public:
    ModelFactory(){
        modelMap["NewbySchwemmer"] = NewbySchwemmer;
    }
    std::shared_ptr<IsoPlanarSOsc> create_Model(std::string name);
    std::shared_ptr<IsoPlanarSOsc::config_t> create_config(std::string name);
    std::shared_ptr<IsoPlanarSOsc::pSet_t> create_pSet(std::string name);

};


#endif //NEWBYSCHWEMMER_MODELFACTORY_H
