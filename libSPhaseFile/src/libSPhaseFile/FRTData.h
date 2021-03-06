//
// Created by konstantin on 1/6/21.
//

#ifndef MRTPHASESIMULATION_FRTDATA_H
#define MRTPHASESIMULATION_FRTDATA_H

#define SPhaseFile_FRTDataFile_CLASS_ID "FRTDataFile"

#include <vector>
#include <list>

#include "SPhaseFile.h"

class FRTDataFile;

// define dataset container holding the data
struct FRTData {
private:
    friend class FRTDataFile;
    string isoSurfaceName;
public:
    // initial positions
    array<vector<double>, 2> x0;
    // first passage times FPT
    vector<double> mFRT;
    vector<double> varFRT;

    FRTData(string isoSurfaceName): isoSurfaceName(isoSurfaceName) {};
    bool operator == (const FRTData& other) const;
    FRTData& operator = (const FRTData& other);
};

class FRTDataFile: public SPhaseFile{
private:
    string modelName;
    fs::path isoSurfaceFilePath;
    fs::path configFilePath;
    list<FRTData> dataList;
protected:
    void read_body(H5::H5File& file);
    void write_body(H5::H5File& file);
public:
    FRTDataFile(string modelName): modelName(modelName),
        SPhaseFile(SPhaseFile_FRTDataFile_CLASS_ID){};
    FRTDataFile(string modelName, fs::path& configFilePath, fs::path& isoSurfaceFilePath):
        modelName(modelName),
        isoSurfaceFilePath(isoSurfaceFilePath),
        configFilePath(configFilePath),
        SPhaseFile(SPhaseFile_FRTDataFile_CLASS_ID){};
    FRTData& createDataSet(string isoSurfaceName);
    bool operator == (const FRTDataFile& other) const;
    FRTDataFile& operator = (const FRTDataFile& other);
};


#endif //MRTPHASESIMULATION_FRTDATA_H
