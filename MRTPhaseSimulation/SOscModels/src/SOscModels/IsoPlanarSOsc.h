//
// Created by konstantin on 4/6/20.
//

#ifndef NEWBYSCHWEMMER_ISOPLANARSOSC_H
#define NEWBYSCHWEMMER_ISOPLANARSOSC_H

#include <array>
#include <map>
#include <string>
#include <random>
#include <memory>

#include <cmath>

#include "mpi.h"
#include "Domain.h"

/*
 * Abstract base class for isotropic, planar
 * stochastic oscillator models
 */
class IsoPlanarSOsc {
private:
    /*
     * Domain
     */
    Domain Disk;
    /*
     * Random Number Generation
     * Normal distributions N(0, 1)
     */
    // Rho
    std::random_device rd_Rho{};
    std::mt19937 gen_Rho{rd_Rho()};
    std::normal_distribution<double> ND_Rho{0.0, 1.0};
    // Phi
    std::random_device rd_Phi{};
    std::mt19937 gen_Phi{rd_Phi()};
    std::normal_distribution<double> ND_Phi{0.0, 1.0};
    /*
     * Parameters
     */
    double D = 0.0; // noise intensity
    /*
     * Initial Condition
     * Initial flow vector [Rho_0, Phi_0, t0]
     */
    std::array<double, 3> x0;
    /*
     * Simulation Parameters
     */
    double dt = 0.0; // time increment
    double T = 0.0; // Simulation time interval: [t0, T]

protected:
    /*
     * Pure virtual functions
     * g - deterministic rho dynamics
     * f - deterministic phase angle dynamics
     * This model is isotropic
     *      -> f and g only depend on the radial coordinate
     */
    virtual double g(double) = 0;
    virtual double f(double) = 0;

public:
    /*
     * Configuration Interface
     */
    struct pSet_t {
        double D;
        virtual ~pSet_t(){}
        virtual void load(std::map<std::string, double>& pMap);
        virtual MPI_Datatype mpiType() = 0;
    };

    struct config_t {
        double dt;
        double T;
        double x0[3];
        MPI_Datatype mpiType();
        void print(){
            std::cout << "IsoPlanarSOsc - config_t: ";
            std::cout << "dt: " << dt;
            std::cout << "| T: " << T;
            std::cout << "| x0: " << x0[0] << ", " << x0[1] << ", " << x0[2] << std::endl;
        }
    };

    /*
    * Iterator class for time evolution
    */
    class IsoPlanarSOscIt {
    private:
        IsoPlanarSOsc& model_ref;
        double t;
    public:
        IsoPlanarSOscIt(IsoPlanarSOsc& model, double t0)
                : model_ref(model), t(t0){}
        bool operator!=(IsoPlanarSOscIt&);
        IsoPlanarSOsc& operator++();
        IsoPlanarSOsc& operator*() const;
    };

    IsoPlanarSOsc() = default;
    IsoPlanarSOsc(Domain&, IsoPlanarSOsc::config_t&, IsoPlanarSOsc::pSet_t*);
    virtual ~IsoPlanarSOsc(){}

    virtual void configure(Domain&, IsoPlanarSOsc::config_t&, IsoPlanarSOsc::pSet_t*);
    std::array<double, 3> evolve();
    std::array<double, 3> get_state() const;
    bool in_time() const;

    IsoPlanarSOscIt begin() {
        return IsoPlanarSOscIt(*this, x0[2]);
    }

    IsoPlanarSOscIt end() {
        return IsoPlanarSOscIt(*this, T);
    }

};

#endif //NEWBYSCHWEMMER_ISOPLANARSOSC_H
