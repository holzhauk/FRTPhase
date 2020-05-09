//
// Created by konstantin on 4/9/20.
//

#ifndef NEWBYSCHWEMMER_NORMMOMENTS_H
#define NEWBYSCHWEMMER_NORMMOMENTS_H

#include <vector>

class NormMoments {
private:
    double mean;
    double variance;
public:
    NormMoments(std::vector<double>&);
    double get_mean();
    double get_variance();
};


#endif //NEWBYSCHWEMMER_NORMMOMENTS_H
