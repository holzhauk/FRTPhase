//
// Created by konstantin on 1/11/21.
//

#ifndef MRTPHASESIMULATION_NORMALIZEDMOMENTS_H
#define MRTPHASESIMULATION_NORMALIZEDMOMENTS_H

#include<vector>
#include<utility>

using namespace std;

struct NormalizedMoments {
public:
    double mean = 0.0;
    double variance = 0.0;
    vector<double> data;
    double get_mean();
    double get_variance();
    tuple<double, double> get_moments();
};


#endif //MRTPHASESIMULATION_NORMALIZEDMOMENTS_H
