//
// Created by konstantin on 4/6/20.
//

#ifndef NEWBYSCHWEMMER_ISOCHRONE_H
#define NEWBYSCHWEMMER_ISOCHRONE_H

#include <string>
#include <map>
#include <vector>
#include <memory>

class Isochrone {
private:
    std::string name;
    std::map<std::string, double> parameters;
    std::unique_ptr<std::vector<double>> Rho;
    std::unique_ptr<std::vector<double>> Phi;

public:
    Isochrone(std::string&);

    void push_parameter(std::string, double&);
    void set_curve(std::vector<double>&, std::vector<double>&);
    std::vector<double>& get_Rho();
    std::vector<double>& get_Phi();
    std::map<std::string, double>& get_parameterMap();
    std::string get_name();
    void print();

};

#endif //NEWBYSCHWEMMER_ISOCHRONE_H
