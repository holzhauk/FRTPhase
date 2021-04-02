//
// Created by konstantin on 3/30/21.
//

#ifndef MRTPHASESIMULATION_TBARDATA_H
#define MRTPHASESIMULATION_TBARDATA_H

#define SPhaseFile_TbarDataFile_CLASS_ID "TbarDataFile"

#include <vector>
#include <list>

#include "SPhaseFile.h"

class TbarDataFile;

class TbarData {
private:
    friend class TbarDataFile;
    string isoSurfaceName;
public:
    double Tbar;
    TbarData(string isoSurfaceName): isoSurfaceName(isoSurfaceName) {};
    bool operator == (const TbarData& other) const;
    TbarData& operator = (const TbarData& other);
};

class TbarDataFile: public SPhaseFile{
private:
    string modelName;
    fs::path isoSurfaceFilePath;
    fs::path configFilePath;
    list<TbarData> dataList;
protected:
    void read_body(H5::H5File& file);
    void write_body(H5::H5File& file);
public:
    TbarDataFile(string modelName): modelName(modelName),
        SPhaseFile(SPhaseFile_TbarDataFile_CLASS_ID){};
    TbarDataFile(string modelName, fs::path& configFilePath, fs::path& isoSurfaceFilePath):
        modelName(modelName),
        isoSurfaceFilePath(isoSurfaceFilePath),
        configFilePath(configFilePath),
        SPhaseFile(SPhaseFile_TbarDataFile_CLASS_ID){};
    TbarData& createDataSet(string isoSurfaceName);
    bool operator == (const TbarDataFile& other) const;
    TbarDataFile& operator = (const TbarDataFile& other);
};


#endif //MRTPHASESIMULATION_TBARDATA_H
