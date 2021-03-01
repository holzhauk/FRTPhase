//
// Created by konstantin on 1/16/21.
//

#ifndef MRTPHASESIMULATION_SAMPLER_H
#define MRTPHASESIMULATION_SAMPLER_H

#include <vector>
#include <array>
#include <libSPhaseFile/libSPhaseFile.h>

using namespace std;

class Sampler {
public:
    virtual vector<array<double, 2>> get_samples(size_t sample_size, InterpolatedCurve&) = 0;
};

#endif //MRTPHASESIMULATION_SAMPLER_H
