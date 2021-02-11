//
// Created by konstantin on 2/7/21.
//

#ifndef MRTPHASESIMULATION_SERIALCORRELATION_H
#define MRTPHASESIMULATION_SERIALCORRELATION_H

#define SPhaseFile_SerialCorrFile_CLASS_ID "SerialCorrelationFile"

#include <string>
#include <list>
#include <filesystem>

#include <H5Cpp.h>

#include "Exceptions.h"
#include "libSPhaseDefinitionsAndConstants.h"
#include "SPhaseFile.h"

namespace fs = std::filesystem;

class CorrelationFile;

struct IsoSurfaceCorr {
public:
    size_t N = 0; // total number of intervals
    size_t offset = 0; // offset
    vector<double> rho_k; // serial correlations up to lag k
private:
    friend class CorrelationFile;
    string key;
public:
    IsoSurfaceCorr(const string& isoSurfaceName): key(isoSurfaceName) {};
    IsoSurfaceCorr& operator = (const IsoSurfaceCorr& other);
    bool operator== (const IsoSurfaceCorr& other) const;
    string get_key();
};

class SerialCorrFile: public SPhaseFile{
private:
    string modelName;
    fs::path isoSurfaceFilePath;
    fs::path configFilePath;
    list<IsoSurfaceCorr> corrList;
protected:
    void write_body(H5::H5File& file);
    void read_body(H5::H5File& file);
public:
    SerialCorrFile(string modelName): modelName(modelName),
        SPhaseFile(SPhaseFile_SerialCorrFile_CLASS_ID) {};
    SerialCorrFile(string modelName,
                    const fs::path& isoSurfaceFilePath,
                    const fs::path& configFilePath):
            modelName(modelName),
            isoSurfaceFilePath(isoSurfaceFilePath),
            configFilePath(configFilePath),
            SPhaseFile(SPhaseFile_SerialCorrFile_CLASS_ID){};
    IsoSurfaceCorr& create_isoSurfaceCorr(const string& isoSurfaceName);
    SerialCorrFile& operator = (const SerialCorrFile& other);
    bool operator== (const SerialCorrFile& other) const;
};


#endif //MRTPHASESIMULATION_SERIALCORRELATION_H
