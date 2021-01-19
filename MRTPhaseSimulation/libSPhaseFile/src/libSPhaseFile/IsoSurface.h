//
// Created by konstantin on 1/3/21.
//

#ifndef MRTPHASESIMULATION_ISOSURFACE_H
#define MRTPHASESIMULATION_ISOSURFACE_H

#define SPhaseFile_IsoSurfaceFile_CLASS_ID "IsoSurfaceFile"

#include <list>
#include <map>
#include <memory>
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
    array< vector<double>, 2 > nodes; // sorted in non-ascending order for rho
public:
    InterpolatedCurve() = default;
    InterpolatedCurve(string name): name(name) {};
    string get_name();
    ParameterSet get_parameterSet();
    tuple<double, double> get_extensions();
    void add_parameter(string pName, double pVal);
    void add_node(array<double, 2>& x);
    void set_nodes(vector<double>& rho, vector<double>& phi);
    tuple<vector<double>, vector<double>> get_nodes();
    tuple<vector<string>, vector<double>> get_parameters();
    bool operator == (const InterpolatedCurve& other) const;
    InterpolatedCurve& operator = (const InterpolatedCurve& other);
    bool is_first_return_event(array<double, 2>& x, bool pos_sense_of_rotation = true);
};

class IsoSurfaceFile: public SPhaseFile {
public:
    class IsoSurfaceFileIt {
        using PtrListIt = list<unique_ptr<InterpolatedCurve>>::iterator;
    private:
        PtrListIt ptrListIt;
    public:
        IsoSurfaceFileIt(PtrListIt ptrListIt): ptrListIt(ptrListIt) {};
        bool operator != (IsoSurfaceFileIt&);
        InterpolatedCurve& operator++();
        InterpolatedCurve& operator*() const;
    };
private:
    string modelName;
    list<unique_ptr<InterpolatedCurve>> curve_ptr_Set;
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
    IsoSurfaceFileIt begin();
    IsoSurfaceFileIt end();
    int get_NoSurfaces();
    string get_modelName();
};

#endif //MRTPHASESIMULATION_ISOSURFACE_H
