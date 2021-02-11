//
// Created by konstantin on 1/3/21.
//

#ifndef MRTPHASESIMULATION_ISOSURFACE_H
#define MRTPHASESIMULATION_ISOSURFACE_H

#define SPhaseFile_IsoSurfaceFile_CLASS_ID "IsoSurfaceFile"

#include <list>
#include <map>
#include <utility>
#include <vector>

#include "SPhaseFile.h"

using namespace std;
namespace fs = filesystem;
using ParameterSet = map<string, double>;

class IsoSurfaceFile;

class InterpolatedCurve {
private:
    friend class IsoSurfaceFile;
    string name;
    ParameterSet parameterSet;
    void sortNodes(); // implicitely called by methods changing nodes
    double omegaBar = 0.0;
    array< vector<double>, 2 > nodes; // sorted in non-ascending order for rho
public:
    InterpolatedCurve() = default;
    InterpolatedCurve(string name): name(name) {};
    string get_name();
    double get_omegaBar();
    ParameterSet get_parameterSet();
    tuple<double, double> get_extensions();
    void add_parameter(string pName, double pVal);
    void set_omegaBar(double val);
    void add_node(array<double, 2>& x);
    void set_nodes(vector<double>& rho, vector<double>& phi);
    tuple<vector<double>, vector<double>> get_nodes();
    tuple<vector<string>, vector<double>> get_parameters();
    bool operator == (const InterpolatedCurve& other) const;
    InterpolatedCurve& operator = (const InterpolatedCurve& other);
    array<double, 2> get_random_point();
    bool is_first_return_event(array<double, 2>& x, bool pos_sense_of_rotation = true);
};

class IsoSurfaceFile: public SPhaseFile {
private:
    string modelName;
    list<InterpolatedCurve> curveList;
protected:
    void read_body(H5::H5File& file);
    void write_body(H5::H5File& file);
public:
    IsoSurfaceFile(): SPhaseFile(SPhaseFile_IsoSurfaceFile_CLASS_ID){};
    IsoSurfaceFile(string modelName):
        modelName(modelName),
        SPhaseFile(SPhaseFile_IsoSurfaceFile_CLASS_ID){};
    InterpolatedCurve& createInterpolatedCurve(string name);
    bool operator == (const IsoSurfaceFile& other) const;
    IsoSurfaceFile& operator = (const IsoSurfaceFile& other);
    list<InterpolatedCurve>::iterator begin();
    list<InterpolatedCurve>::iterator end();
    int get_NoSurfaces();
    string get_modelName();
};

#endif //MRTPHASESIMULATION_ISOSURFACE_H
