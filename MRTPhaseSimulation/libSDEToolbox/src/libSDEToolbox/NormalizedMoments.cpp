//
// Created by konstantin on 1/11/21.
//

#include <cmath>
#include "NormalizedMoments.h"

double NormalizedMoments::get_mean() {
    mean = 0.0;
    for (auto s : data) {
        mean += s;
    }
    mean /= data.size();
    return mean;
}

double NormalizedMoments::get_variance() {
    mean = this->get_mean();
    variance = 0.0;
    for (auto s : data){
        variance += pow(s - mean, 2.0);
    }
    variance /= data.size();
}