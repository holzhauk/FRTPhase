//
// Created by konstantin on 1/16/21.
//

#include "EquidistantSampler.h"

vector<array<double,2>> EquidistantSampler::get_samples(size_t sample_size, InterpolatedCurve &curve) {
    auto [rhos, phis] = curve.get_nodes();
    vector<array<double, 2>> samples;
    size_t n = (size_t) rhos.size() / sample_size;
    for (int i = 0; i < rhos.size(); i += n){
        samples.push_back({rhos[i], phis[i]});
    }
    return samples;
}
