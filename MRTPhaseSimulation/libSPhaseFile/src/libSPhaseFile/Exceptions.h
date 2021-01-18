//
// Created by konstantin on 1/4/21.
//

#ifndef MRTPHASESIMULATION_EXCEPTIONS_H
#define MRTPHASESIMULATION_EXCEPTIONS_H

#include <exception>
#include <string>

using namespace std;

struct SPhaseFileWrongFormat : public exception {
    const char* what () const throw() {
        return "SPhaseFile - wrong file format";
    }
};

struct SPhaseFileVersionConflict : public exception {
    const char* what () const throw() {
        return "SPhaseFile - file was created with a different version of SPhaseFile";
    }
};

struct SPhaseFileModelConflict : public exception {
    const char* what () const throw() {
        return "SPhaseFile - file contains information about a different model than specified";
    }
};

struct SPhaseFileClassConflict : public exception {
    const char* what () const throw() {
        return "SPhaseFile - the file contains information about a different class";
    }
};

struct SPhaseFileDimError : public exception {
    const char* whar () const throw() {
        return "SPhaseFile - the extensions of one object are faulty";
    }
};

struct SPhaseFileOutOfDomain : public exception {
    const char* whar () const throw() {
        return "SPhaseFile - state does not lie within the domain of the isosurface curve ";
    }
};

struct SPhaseFileHDF5APIError : public exception {
    string ErrorMsg;

    SPhaseFileHDF5APIError(string msg) : exception() {
        ErrorMsg = msg;
    }

    const char* what () const throw() {
        return ErrorMsg.c_str();
    }
};

#endif //MRTPHASESIMULATION_EXCEPTIONS_H
