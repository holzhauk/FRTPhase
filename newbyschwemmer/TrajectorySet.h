//
// Created by konstantin on 4/9/20.
//

#ifndef NEWBYSCHWEMMER_TRAJECTORYSET_H
#define NEWBYSCHWEMMER_TRAJECTORYSET_H

#include <filesystem>
#include <map>
#include <string>
#include <list>

#include "H5Cpp.h"

#include "Trajectory.h"

namespace fs = std::filesystem;
using namespace H5;

const H5std_string FORMAT_TAG_TRAJECTORY_SET = "SOscFile";
const H5std_string VERSION_TAG_TRAJECTORY_SET = "0.0.1";
const H5std_string CLASS_TAG_TRAJECTORY_SET = "TrajectorySet";

class TrajectorySet {
private:
    std::map<std::string, double> parameters;
    std::string model_name;
    std::list<Trajectory*> Trajectory_ptr_list;

public:
    TrajectorySet(std::map<std::string, double>&, std::string&);
    ~TrajectorySet();

    void push_back(Trajectory*);
    void write_to_file(fs::path);
};


#endif //NEWBYSCHWEMMER_TRAJECTORYSET_H
