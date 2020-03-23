// Unix standard libraries
#include <sys/types.h>
#include <unistd.h>

// CXX standard libraries
#include <iostream>
#include <ctime>

// GNU Scientific Library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// HDF5 libraries
#include "hdf5.h"
#define FILENAME "file.h5"

unsigned long gen_seed();

int main() {
    std::cout << "Hello, World!" << std::endl;

    double dt = 0.1;

    // initialize random number generator "Mersenne Twister"
    gsl_rng* rng_ptr = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng_ptr, gen_seed());
    double dW = gsl_ran_gaussian(rng_ptr, dt);

    std::cout << "dW = " << dW << std::endl;

    gsl_rng_free(rng_ptr);

    hid_t file_id; // file handler / identifier
    herr_t status; // error status variable

    // create a new file using default properties
    file_id = H5Fcreate(FILENAME,H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // close file
    status = H5Fclose(file_id);

    return EXIT_SUCCESS;
}

unsigned long gen_seed(){
    pid_t pid = getpid(); // get id of current process
    unsigned long t = time(NULL);
    return t*pid;
}

