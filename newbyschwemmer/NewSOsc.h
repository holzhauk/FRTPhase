//
// Created by konstantin on 4/7/20.
//

#ifndef NEWBYSCHWEMMER_NEWSOSC_H
#define NEWBYSCHWEMMER_NEWSOSC_H

#include "IsoPlanarSOsc.h"

typedef struct pset_NewS_t{
    double D;
    double omega;
    double gamma;
    double c;
} pset_NewS_t;

class NewSOsc : public IsoPlanarSOsc {
protected:
    double omega = 0.0;
    double gamma = 0.0;
    double c = 0.0;

    /*
     * define previously pure virtual functions
     */
    double g(double); // deterministic rho-dynamics
    double f(double); // deterministic phi-dynamics

public:

    NewSOsc() = default;
    NewSOsc(Domain&);
    NewSOsc(Domain&, sim_config_t&);
    NewSOsc(Domain&, pset_NewS_t&, sim_config_t&);
    NewSOsc(Domain&, sim_config_C_t&);
    NewSOsc(Domain&, pset_NewS_t&, sim_config_C_t&);

    void configure(Domain&, pset_NewS_t&, sim_config_t&);
    void configure(Domain&, pset_NewS_t&, sim_config_C_t&);
};


#endif //NEWBYSCHWEMMER_NEWSOSC_H
