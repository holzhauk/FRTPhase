//
// Created by konstantin on 1/6/21.
//

#ifndef MRTPHASESIMULATION_SPHASEFILE_H
#define MRTPHASESIMULATION_SPHASEFILE_H

#include <string>
#include <filesystem>
#include <H5Cpp.h>

#include "Exceptions.h"
#include "libSPhaseDefinitionsAndConstants.h"

namespace fs = filesystem;

class SPhaseFile {
private:
    const string SPHASEFILE_CLASS_ID;
protected:
    virtual void read_body(H5::H5File& file) = 0;
    virtual void write_body(H5::H5File& file) = 0;
public:
    SPhaseFile(string class_id): SPHASEFILE_CLASS_ID(class_id) {};
    void read(const fs::path& filePath);
    void write(const fs::path& filePath);
};


#endif //MRTPHASESIMULATION_SPHASEFILE_H
