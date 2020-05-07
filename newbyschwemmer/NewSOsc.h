//
// Created by konstantin on 4/7/20.
//

#ifndef NEWBYSCHWEMMER_NEWSOSC_H
#define NEWBYSCHWEMMER_NEWSOSC_H

#include "IsoPlanarSOsc.h"

#include <iostream>
#include <memory>

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

    struct pSet_t : IsoPlanarSOsc::pSet_t {
        double omega;
        double gamma;
        double c;
        void load(std::map<std::string, double>& pMap) override;
        MPI_Datatype mpiType() override;
    };

    NewSOsc() = default;
    NewSOsc(Domain&);
    NewSOsc(Domain&, NewSOsc::config_t&, IsoPlanarSOsc::pSet_t*);

    void configure(Domain&, NewSOsc::config_t&, IsoPlanarSOsc::pSet_t*);
};


#endif //NEWBYSCHWEMMER_NEWSOSC_H
