//
// Created by konstantin on 1/16/21.
//

#ifndef MRTPHASESIMULATION_EQUIDISTANTSAMPLER_H
#define MRTPHASESIMULATION_EQUIDISTANTSAMPLER_H

#include <libSPhaseFile/libSPhaseFile.h>
#include "Sampler.h"

class EquidistantSampler: public Sampler {
private:

public:
    vector<array<double, 2>> get_samples(size_t sample_size, InterpolatedCurve& curve);
};


#endif //MRTPHASESIMULATION_EQUIDISTANTSAMPLER_H
