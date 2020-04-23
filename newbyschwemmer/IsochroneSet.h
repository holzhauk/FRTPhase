//
// Created by konstantin on 4/6/20.
//

#ifndef NEWBYSCHWEMMER_ISOCHRONESET_H
#define NEWBYSCHWEMMER_ISOCHRONESET_H

#include <filesystem>
#include <list>

#include "H5Cpp.h"

#include "Isochrone.h"

/*
 * C- operator functions to iterate through the HDF5 File
 *  - not covered by the H5Cpp API yet
 *  -> using the C-Api for this part
 */
extern "C" herr_t isochrone_handler (hid_t, const char*,
                                     const H5L_info_t* , void* );

extern "C" herr_t parameter_handler (hid_t, const char*,
                                     const H5L_info_t*, void* );

namespace fs = std::filesystem;
using namespace H5;

class IsochroneSet {
private:
    std::list<Isochrone*> Isochrone_list_ptr;

public:
    ~IsochroneSet();

    void load(fs::path);
    std::list<Isochrone*>::iterator begin();
    std::list<Isochrone*>::iterator end();
    int get_NoIsochrones();

};

#endif //NEWBYSCHWEMMER_ISOCHRONESET_H
