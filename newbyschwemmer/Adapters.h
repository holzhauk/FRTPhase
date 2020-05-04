//
// Created by konstantin on 5/4/20.
//

#ifndef NEWBYSCHWEMMER_ADAPTERS_H
#define NEWBYSCHWEMMER_ADAPTERS_H

#include "IsoPlanarSOsc.h"
#include "SimConfig_t.h"


namespace Adapters {

    struct Config2SimConfig : IsoPlanarSOsc::config_t {
        Config2SimConfig(SimConfig_t& simConfig);
    };

}


#endif //NEWBYSCHWEMMER_ADAPTERS_H
