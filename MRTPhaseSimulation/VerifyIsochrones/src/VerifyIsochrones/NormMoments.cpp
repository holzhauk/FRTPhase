//
// Created by konstantin on 4/9/20.
//
#include <numeric>
#include <cmath>

#include "NormMoments.h"

NormMoments::NormMoments(std::vector<double>& samples) {
    mean = 0.0;
    for (auto s : samples){
        mean += s;
    }
    mean /= samples.size();

    variance = 0.0;
    for (auto s : samples){
        variance += pow(s - mean, 2.0);
    }
    variance /= samples.size();
}

double NormMoments::get_mean() {
    return mean;
}

double NormMoments::get_variance() {
    return variance;
}