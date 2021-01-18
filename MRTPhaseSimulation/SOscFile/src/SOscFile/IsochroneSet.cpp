//
// Created by konstantin on 4/6/20.
//

#include "IsochroneSet.h"

void IsochroneSet::load(fs::path filepath) {
    try {
        Exception::dontPrint();
        H5std_string FILENAME(filepath);
        H5File* file = new H5File(FILENAME, H5F_ACC_RDONLY);
        herr_t idx = H5Literate(file->getId(), H5_INDEX_NAME, H5_ITER_INC,
                                NULL, isochrone_handler, &Isochron_list_ptr);
        delete file;
    }
        // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printErrorStack();
    }
        // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printErrorStack();
    }
        // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printErrorStack();
    }
        // catch failure caused by the Attribute operations
    catch( AttributeIException error )
    {
        error.printErrorStack();
    }
}

std::list<std::shared_ptr<Isochron>>::iterator IsochroneSet::begin() {
    return Isochron_list_ptr.begin();
}

std::list<std::shared_ptr<Isochron>>::iterator IsochroneSet::end() {
    return Isochron_list_ptr.end();
}

/***************************************
 * C - operator function definitions
 ***************************************/

/*
 * Operator function to iterate over the Isochrone groups in the hdf5 file
 */
herr_t isochrone_handler(hid_t loc_id, const char* name, const H5L_info_t* linfo, void* opdata){

    auto Isochrone_list_ptr = reinterpret_cast<std::list<std::shared_ptr<Isochron>>*>(opdata);

    std::string name_s(name); // convert C-String to std::string object
    std::shared_ptr<Isochron> Iso = std::shared_ptr<Isochron>(new Isochron(name_s));

    hid_t iso_g = H5Gopen2(loc_id, name, H5P_DEFAULT);

    /*
     * get the parameters in the Parameters subgroup
     */
    hid_t pset_sg = H5Gopen1(iso_g, "Parameters");
    herr_t idx = H5Literate(pset_sg, H5_INDEX_NAME, H5_ITER_INC,
                            NULL, parameter_handler, Iso.get());
    H5Gclose(pset_sg);

    /*
     * get the coordinates of the Isochrone stored in the Curve subgroup
     */
    std::vector<double>* Rhos;
    std::vector<double>* Phis;
    hid_t curve_sg = H5Gopen1(iso_g, "Curve");

    // get radial coordinates
    hid_t dset_rho = H5Dopen(curve_sg, "Rho", H5P_DEFAULT);
    hid_t space_rho = H5Dget_space(dset_rho);
    hsize_t size_dim_rho[2];
    int ndims_rho = H5Sget_simple_extent_dims(space_rho, size_dim_rho, NULL);
    Rhos = new std::vector<double>(size_dim_rho[0]);
    herr_t status_rho = H5Dread(dset_rho, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Rhos->data());
    H5Sclose(space_rho);
    H5Dclose(dset_rho);

    // get angular coordinates
    hid_t dset_phi = H5Dopen(curve_sg, "Phi", H5P_DEFAULT);
    hid_t space_phi = H5Dget_space(dset_phi);
    hsize_t size_dim0_phi[2];
    int ndims_phi = H5Sget_simple_extent_dims(space_phi, size_dim0_phi, NULL);
    Phis = new std::vector<double>(size_dim0_phi[0]);
    herr_t status_phi = H5Dread(dset_phi, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Phis->data());
    H5Sclose(space_phi);
    H5Dclose(dset_phi);

    Iso->set_curve(*Rhos, *Phis);

    delete Rhos;
    delete Phis;
    H5Gclose(curve_sg);

    H5Gclose(iso_g);
    Isochrone_list_ptr->push_back(Iso);
    return 0;

}

/*
 * Operator function to iterate over the datasets in the Parameters group in the Isochrone group
 */
herr_t parameter_handler(hid_t loc_id, const char* name, const H5L_info_t* linfo, void* opdata) {

    auto Iso = reinterpret_cast<Isochron*>(opdata);

    hid_t dset = H5Dopen(loc_id, name, H5P_DEFAULT);
    hid_t space = H5Dget_space(dset);

    hsize_t size_dim0[2];
    int ndims = H5Sget_simple_extent_dims(space, size_dim0, NULL);
    double val;
    herr_t idx = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
    Iso->push_parameter(name, val);

    H5Sclose(space);
    H5Dclose(dset);

    return 0;

}
