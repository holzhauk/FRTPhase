//
// Created by konstantin on 4/6/20.
//

#ifndef NEWBYSCHWEMMER_ISOPLANARSOSC_H
#define NEWBYSCHWEMMER_ISOPLANARSOSC_H

#include <array>
#include <random>

#include <cmath>

#include "Domain.h"

// Simulation configuration
typedef struct sim_config_t {
    double dt; // time increment
    double T; // total time
    std::array<double, 3> x0; // initial condition: [Rho_0, Phi_0, t0]
} sim_config_t;

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
    IsoPlanarSOsc(Domain&);
    IsoPlanarSOsc(Domain&, sim_config_t&);
    IsoPlanarSOsc(Domain&, sim_config_t&, double);

    void configure(Domain&, sim_config_t&, double);
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
