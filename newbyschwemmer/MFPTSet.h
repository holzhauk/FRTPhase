//
// Created by konstantin on 4/8/20.
//

#ifndef NEWBYSCHWEMMER_MFPTSET_H
#define NEWBYSCHWEMMER_MFPTSET_H

#include <filesystem>
#include <list>
#include <string>

#include "H5Cpp.h"

#include "MFPTs.h"

namespace fs = std::filesystem;
using namespace H5;

const H5std_string FORMAT_TAG_MFPTS = "SOscFile";
const H5std_string VERSION_TAG_MFPTS = "0.0.1";
const H5std_string CLASS_TAG_MFPTS = "MeanFirstPassageTimesSet";

class MFPTSet {
private:
    std::string model_name;
    fs::path IsochroneFile_name;
    std::list<MFPTs*> MFPTs_ptr_list;
public:
    MFPTSet(std::string&, fs::path&);
    ~MFPTSet();
    void push_back(MFPTs*);
    void write_to_file(fs::path);
};


#endif //NEWBYSCHWEMMER_MFPTSET_H
