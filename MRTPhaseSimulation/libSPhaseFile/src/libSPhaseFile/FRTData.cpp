//
// Created by konstantin on 1/6/21.
//

#include "FRTData.h"
#include <iostream>

extern "C" herr_t surface_group_handler_frtdata(hid_t, const char*, const H5L_info_t*, void*);

bool FRTData::operator == (const FRTData& other) const {
    bool is_equal = true;
    is_equal = is_equal && (this->isoSurfaceName == other.isoSurfaceName);
    is_equal = is_equal && (this->x0 == other.x0);
    is_equal = is_equal && (this->mPhiT == other.mPhiT);
    is_equal = is_equal && (this->varPhiT == other.varPhiT);
    is_equal = is_equal && (this->mT == other.mT);
    is_equal = is_equal && (this->varT == other.varT);
    is_equal = is_equal && (this->mFRT == other.mFRT);
    is_equal = is_equal && (this->varFRT == other.varFRT);
    return is_equal;
}



FRTData& FRTDataFile::createDataSet(string isoSurfaceName) {
    unique_ptr<FRTData> dataSet_ptr(new FRTData(isoSurfaceName));
    data_ptr_Set.push_back(move(dataSet_ptr));
    return *data_ptr_Set.back();
}

void FRTDataFile::write_body(H5::H5File& file) {

    // define attribute types
    H5::StrType str_type(0, H5T_VARIABLE);
    str_type.setCset(H5T_CSET_UTF8);
    str_type.setStrpad(H5T_STR_NULLTERM);
    H5::DataSpace dspace(H5S_SCALAR);

    // write model name attribute
    H5::Attribute attr = file.createAttribute("model", str_type, dspace);
    attr.write(str_type, modelName);
    attr.close();

    // write configuration file path
    H5::DataSet dset = file.createDataSet("configuration_file", str_type, dspace);
    dset.write(configFilePath, str_type);
    dset.close();

    // write isosurface file path
    dset = file.createDataSet("isosurface_file", str_type, dspace);
    dset.write(isoSurfaceFilePath, str_type);
    dset.close();

    dspace.close();
    str_type.close();

    // write the individual data sets to the file
    for (auto& data_ptr: data_ptr_Set){

        hsize_t size = data_ptr->x0[0].size();
        if ((data_ptr->x0[0].size() != size) ||
            (data_ptr->mPhiT.size() != size) ||
            (data_ptr->varPhiT.size() != size) ||
            (data_ptr->mT.size() != size) ||
            (data_ptr->varT.size() != size) ||
            (data_ptr->mFRT.size() != size) ||
            (data_ptr->varFRT.size() != size)) throw SPhaseFileDimError();

        H5::DataSpace dspace;

        H5::Group surface = file.createGroup(data_ptr->isoSurfaceName);

        // INITIAL POSITIONS
        H5::Group xinit_g = surface.createGroup("initial_position");
        hsize_t dim_size[] = {data_ptr->x0[0].size()};
        dspace = H5::DataSpace(1, dim_size);
        H5::DataSet rhos = xinit_g.createDataSet("rho",
                                               H5::PredType::NATIVE_DOUBLE,
                                               dspace);
        rhos.write(data_ptr->x0[0].data(), H5::PredType::NATIVE_DOUBLE);
        rhos.close();
        H5::DataSet phis = xinit_g.createDataSet("phi",
                                               H5::PredType::NATIVE_DOUBLE,
                                               dspace);
        phis.write(data_ptr->x0[1].data(), H5::PredType::NATIVE_DOUBLE);
        phis.close();
        xinit_g.close();

        // ASYMPTOTIC ANGLE PHIT
        H5::Group aangle_g = surface.createGroup("asymptotic_angle");
        H5::DataSet mPhiTs = aangle_g.createDataSet("mPhiT",
                                                    H5::PredType::NATIVE_DOUBLE,
                                                    dspace);
        mPhiTs.write(data_ptr->mPhiT.data(), H5::PredType::NATIVE_DOUBLE);
        mPhiTs.close();
        H5::DataSet varPhiTs = aangle_g.createDataSet("varPhiT",
                                                      H5::PredType::NATIVE_DOUBLE,
                                                      dspace);
        varPhiTs.write(data_ptr->varPhiT.data(), H5::PredType::NATIVE_DOUBLE);
        varPhiTs.close();
        aangle_g.close();

        // STATIONARY PERIOD T
        H5::Group period_g = surface.createGroup("stationary_period");
        H5::DataSet mTs = period_g.createDataSet("mT",
                                                 H5::PredType::NATIVE_DOUBLE,
                                                 dspace);
        mTs.write(data_ptr->mT.data(), H5::PredType::NATIVE_DOUBLE);
        mTs.close();
        H5::DataSet varTs = period_g.createDataSet("varT",
                                                   H5::PredType::NATIVE_DOUBLE,
                                                   dspace);
        varTs.write(data_ptr->varT.data(), H5::PredType::NATIVE_DOUBLE);
        varTs.close();
        period_g.close();

        // FIRST RETURN TIMES FRT
        H5::Group frt_g = surface.createGroup("first_return_time");
        H5::DataSet mFRTs = frt_g.createDataSet("mFRT",
                                                H5::PredType::NATIVE_DOUBLE,
                                                dspace);
        mFRTs.write(data_ptr->mFRT.data(), H5::PredType::NATIVE_DOUBLE);
        mFRTs.close();
        H5::DataSet varFRTs = frt_g.createDataSet("varFRT",
                                                  H5::PredType::NATIVE_DOUBLE,
                                                  dspace);
        varFRTs.write(data_ptr->varFRT.data(), H5::PredType::NATIVE_DOUBLE);
        varFRTs.close();

        dspace.close();
        surface.close();

    }
}

void FRTDataFile::read_body(H5::H5File &file) {

    H5std_string str_buf;

    H5::Attribute attr = file.openAttribute("model");
    H5::DataType dtype = attr.getDataType();
    attr.read(dtype, str_buf);
    if (str_buf != modelName)
        throw SPhaseFileModelConflict();
    dtype.close();
    attr.close();

    H5::DataSet dset = file.openDataSet("configuration_file");
    dtype = dset.getDataType();
    dset.read(str_buf, dtype);
    configFilePath = string(str_buf);
    dtype.close();
    dset.close();

    dset = file.openDataSet("isosurface_file");
    dtype = dset.getDataType();
    dset.read(str_buf, dtype);
    isoSurfaceFilePath = string(str_buf);
    dtype.close();
    dset.close();

    herr_t idx = H5Literate(file.getId(),
                            H5_INDEX_NAME,
                            H5_ITER_INC,
                            NULL,
                            surface_group_handler_frtdata,
                            this);
    if (idx < 0)
        throw SPhaseFileHDF5APIError("H5Literate encountered failure iterating through surface groups");
}

herr_t surface_group_handler_frtdata(hid_t group_id, const char* group_name,
                             const H5L_info_t* group_info, void* opdata){

    string g_name(group_name);
    if ((g_name == "isosurface_file") || (g_name == "configuration_file"))
        return 0;

    auto dataFile_ptr = reinterpret_cast<FRTDataFile*> (opdata);

    // using C - interface from now on
    // open the surface group that are iterated through
    hid_t surface_group_id = H5Gopen(group_id, group_name, H5P_DEFAULT);
    if (surface_group_id < 0)
        throw SPhaseFileHDF5APIError(string("failed to open group: ") + string(group_name));

    ////////////////////////////////////////////////////////////////////
    // open initial_position subgroup to retrieve the node's coordinates
    hid_t xinit_subgroup_id = H5Gopen(surface_group_id, "initial_position", H5P_DEFAULT);
    if (xinit_subgroup_id < 0)
        throw SPhaseFileHDF5APIError("failed to open 'initial_position' subgroup in surface group: "
                                     + string(group_name));

    // get radial coordinates
    hid_t dset_rho_id = H5Dopen(xinit_subgroup_id, "rho", H5P_DEFAULT);
    if (dset_rho_id < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset 'rho' in surface group: "
                                     + string(group_name));
    hid_t space_rho_id = H5Dget_space(dset_rho_id);
    if (space_rho_id < 0)
        throw SPhaseFileHDF5APIError("failed to get dataspace from 'rho' in surface group: "
                                     + string(group_name));
    hsize_t size_dim_rho[2];
    int ndims_rho = H5Sget_simple_extent_dims(space_rho_id, size_dim_rho, NULL);
    if (ndims_rho < 0)
        throw SPhaseFileHDF5APIError("failed to get the dspace's dims from 'rho' in surface group: "
                                     + string(group_name));
    vector<double> rho_dummy_buf(size_dim_rho[0]);
    herr_t status_rho = H5Dread(dset_rho_id,
                                H5T_NATIVE_DOUBLE,
                                H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                rho_dummy_buf.data());
    if (status_rho < 0)
        throw SPhaseFileHDF5APIError("failed to read dataset 'rho' in surface group: "
                                     + string(group_name));
    H5Sclose(space_rho_id);
    H5Dclose(dset_rho_id);

    // get angular coordinates
    hid_t dset_phi_id = H5Dopen(xinit_subgroup_id, "phi", H5P_DEFAULT);
    if (dset_rho_id < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset 'phi' in surface group: "
                                     + string(group_name));
    hid_t space_phi_id = H5Dget_space(dset_phi_id);
    if (space_rho_id < 0)
        throw SPhaseFileHDF5APIError("failed to get dataspace from 'phi' in surface group: "
                                     + string(group_name));
    hsize_t size_dim_phi[2];
    int ndims_phi = H5Sget_simple_extent_dims(space_phi_id, size_dim_phi, NULL);
    if (ndims_phi < 0)
        throw SPhaseFileHDF5APIError("failed to get the dspace's dims from 'phi' in surface group: "
                                     + string(group_name));
    vector<double> phi_dummy_buf(size_dim_phi[0]);
    herr_t status_phi = H5Dread(dset_phi_id,
                                H5T_NATIVE_DOUBLE,
                                H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                phi_dummy_buf.data());
    if (status_phi < 0)
        throw SPhaseFileHDF5APIError("failed to read dataset 'phi' in surface group: "
                                     + string(group_name));
    H5Sclose(space_phi_id);
    H5Dclose(dset_phi_id);

    herr_t xinit_subgroup_err = H5Gclose(xinit_subgroup_id);
    if (xinit_subgroup_err < 0)
        throw SPhaseFileHDF5APIError("failed to close 'initial_position' subgroup in surface group: "
                                     + string(group_name));

    //////////////////////////////////////////////////////////
    // open asymptotic angle subgroup to retrieve statistics
    hid_t phiT_subgroup_id = H5Gopen(surface_group_id, "asymptotic_angle", H5P_DEFAULT);
    if (phiT_subgroup_id < 0)
        throw SPhaseFileHDF5APIError("failed to open 'asymptotic_angle' subgroup in surface group: "
                                     + string(group_name));

    // get mean values
    hid_t dset_mPhiT_id = H5Dopen(phiT_subgroup_id, "mPhiT", H5P_DEFAULT);
    if (dset_mPhiT_id < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset 'mPhiT' in surface group: "
                                     + string(group_name));
    hid_t space_mPhiT_id = H5Dget_space(dset_mPhiT_id);
    if (space_mPhiT_id < 0)
        throw SPhaseFileHDF5APIError("failed to get dataspace from 'mPhiT' in surface group: "
                                     + string(group_name));
    hsize_t size_dim_mPhiT[2];
    int ndims_mPhiT = H5Sget_simple_extent_dims(space_mPhiT_id, size_dim_mPhiT, NULL);
    if (ndims_mPhiT < 0)
        throw SPhaseFileHDF5APIError("failed to get the dspace's dims from 'mPhiT' in surface group: "
                                     + string(group_name));
    vector<double> mPhiT_dummy_buf(size_dim_mPhiT[0]);
    herr_t status_mPhiT = H5Dread(dset_mPhiT_id,
                               H5T_NATIVE_DOUBLE,
                               H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               mPhiT_dummy_buf.data());
    if (status_mPhiT < 0)
        throw SPhaseFileHDF5APIError("failed to read dataset 'mPhiT' in surface group: "
                                     + string(group_name));
    H5Sclose(space_mPhiT_id);
    H5Dclose(dset_mPhiT_id);

    // get variances
    hid_t dset_varPhiT_id = H5Dopen(phiT_subgroup_id, "varPhiT", H5P_DEFAULT);
    if (dset_varPhiT_id < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset 'varPhiT' in surface group: "
                                     + string(group_name));
    hid_t space_varPhiT_id = H5Dget_space(dset_varPhiT_id);
    if (space_varPhiT_id < 0)
        throw SPhaseFileHDF5APIError("failed to get dataspace from 'varPhiT' in surface group: "
                                     + string(group_name));
    hsize_t size_dim_varPhiT[2];
    int ndims_varPhiT = H5Sget_simple_extent_dims(space_varPhiT_id, size_dim_varPhiT, NULL);
    if (ndims_varPhiT < 0)
        throw SPhaseFileHDF5APIError("failed to get the dspace's dims from 'varPhiT' in surface group: "
                                     + string(group_name));
    vector<double> varPhiT_dummy_buf(size_dim_varPhiT[0]);
    herr_t status_varPhiT = H5Dread(dset_varPhiT_id,
                                 H5T_NATIVE_DOUBLE,
                                 H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                 varPhiT_dummy_buf.data());
    if (status_varPhiT < 0)
        throw SPhaseFileHDF5APIError("failed to read dataset 'varPhiT' in surface group: "
                                     + string(group_name));
    H5Sclose(space_varPhiT_id);
    H5Dclose(dset_varPhiT_id);

    herr_t phiT_subgroup_err = H5Gclose(phiT_subgroup_id);
    if (phiT_subgroup_err < 0)
        throw SPhaseFileHDF5APIError("failed to close 'asymptotic_angle' subgroup in surface group: "
                                     + string(group_name));


    //////////////////////////////////////////////////////////
    // open stationary period subgroup to retrieve statistics
    hid_t tbar_subgroup_id = H5Gopen(surface_group_id, "stationary_period", H5P_DEFAULT);
    if (tbar_subgroup_id < 0)
        throw SPhaseFileHDF5APIError("failed to open 'stationary_period' subgroup in surface group: "
                                     + string(group_name));

    // get mean values
    hid_t dset_mT_id = H5Dopen(tbar_subgroup_id, "mT", H5P_DEFAULT);
    if (dset_mT_id < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset 'mT' in surface group: "
                                     + string(group_name));
    hid_t space_mT_id = H5Dget_space(dset_mT_id);
    if (space_mT_id < 0)
        throw SPhaseFileHDF5APIError("failed to get dataspace from 'mT' in surface group: "
                                     + string(group_name));
    hsize_t size_dim_mT[2];
    int ndims_mT = H5Sget_simple_extent_dims(space_mT_id, size_dim_mT, NULL);
    if (ndims_mT < 0)
        throw SPhaseFileHDF5APIError("failed to get the dspace's dims from 'mT' in surface group: "
                                     + string(group_name));
    vector<double> mT_dummy_buf(size_dim_mT[0]);
    herr_t status_mT = H5Dread(dset_mT_id,
                                H5T_NATIVE_DOUBLE,
                                H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                mT_dummy_buf.data());
    if (status_mT < 0)
        throw SPhaseFileHDF5APIError("failed to read dataset 'mT' in surface group: "
                                     + string(group_name));
    H5Sclose(space_mT_id);
    H5Dclose(dset_mT_id);

    // get variances
    hid_t dset_varT_id = H5Dopen(tbar_subgroup_id, "varT", H5P_DEFAULT);
    if (dset_varT_id < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset 'varT' in surface group: "
                                     + string(group_name));
    hid_t space_varT_id = H5Dget_space(dset_varT_id);
    if (space_varT_id < 0)
        throw SPhaseFileHDF5APIError("failed to get dataspace from 'varT' in surface group: "
                                     + string(group_name));
    hsize_t size_dim_varT[2];
    int ndims_varT = H5Sget_simple_extent_dims(space_varT_id, size_dim_varT, NULL);
    if (ndims_varT < 0)
        throw SPhaseFileHDF5APIError("failed to get the dspace's dims from 'varT' in surface group: "
                                     + string(group_name));
    vector<double> varT_dummy_buf(size_dim_varT[0]);
    herr_t status_varT = H5Dread(dset_varT_id,
                                H5T_NATIVE_DOUBLE,
                                H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                varT_dummy_buf.data());
    if (status_varT < 0)
        throw SPhaseFileHDF5APIError("failed to read dataset 'varT' in surface group: "
                                     + string(group_name));
    H5Sclose(space_varT_id);
    H5Dclose(dset_varT_id);

    herr_t tbar_subgroup_err = H5Gclose(tbar_subgroup_id);
    if (tbar_subgroup_err < 0)
        throw SPhaseFileHDF5APIError("failed to close 'stationary_period' subgroup in surface group: "
                                     + string(group_name));

    //////////////////////////////////////////////////////////
    // open first return time subgroup to retrieve statistics
    hid_t frt_subgroup_id = H5Gopen(surface_group_id, "first_return_time", H5P_DEFAULT);
    if (frt_subgroup_id < 0)
        throw SPhaseFileHDF5APIError("failed to open 'first_return_time' subgroup in surface group: "
                                     + string(group_name));

    // get mean values
    hid_t dset_mFRT_id = H5Dopen(frt_subgroup_id, "mFRT", H5P_DEFAULT);
    if (dset_mFRT_id < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset 'mFRT' in surface group: "
                                     + string(group_name));
    hid_t space_mFRT_id = H5Dget_space(dset_mFRT_id);
    if (space_mFRT_id < 0)
        throw SPhaseFileHDF5APIError("failed to get dataspace from 'mFRT' in surface group: "
                                     + string(group_name));
    hsize_t size_dim_mFRT[2];
    int ndims_mFRT = H5Sget_simple_extent_dims(space_mFRT_id, size_dim_mFRT, NULL);
    if (ndims_mFRT < 0)
        throw SPhaseFileHDF5APIError("failed to get the dspace's dims from 'mFRT' in surface group: "
                                     + string(group_name));
    vector<double> mFRT_dummy_buf(size_dim_mFRT[0]);
    herr_t status_mFRT = H5Dread(dset_mFRT_id,
                               H5T_NATIVE_DOUBLE,
                               H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               mFRT_dummy_buf.data());
    if (status_mFRT < 0)
        throw SPhaseFileHDF5APIError("failed to read dataset 'mFRT' in surface group: "
                                     + string(group_name));
    H5Sclose(space_mFRT_id);
    H5Dclose(dset_mFRT_id);

    // get variances
    hid_t dset_varFRT_id = H5Dopen(frt_subgroup_id, "varFRT", H5P_DEFAULT);
    if (dset_varFRT_id < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset 'varFRT' in surface group: "
                                     + string(group_name));
    hid_t space_varFRT_id = H5Dget_space(dset_varFRT_id);
    if (space_varFRT_id < 0)
        throw SPhaseFileHDF5APIError("failed to get dataspace from 'varFRT' in surface group: "
                                     + string(group_name));
    hsize_t size_dim_varFRT[2];
    int ndims_varFRT = H5Sget_simple_extent_dims(space_varFRT_id, size_dim_varFRT, NULL);
    if (ndims_varFRT < 0)
        throw SPhaseFileHDF5APIError("failed to get the dspace's dims from 'varFRT' in surface group: "
                                     + string(group_name));
    vector<double> varFRT_dummy_buf(size_dim_varFRT[0]);
    herr_t status_varFRT = H5Dread(dset_varFRT_id,
                                 H5T_NATIVE_DOUBLE,
                                 H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                 varFRT_dummy_buf.data());
    if (status_varFRT < 0)
        throw SPhaseFileHDF5APIError("failed to read dataset 'varFRT' in surface group: "
                                     + string(group_name));
    H5Sclose(space_varFRT_id);
    H5Dclose(dset_varFRT_id);

    herr_t frt_subgroup_err = H5Gclose(frt_subgroup_id);
    if (frt_subgroup_err < 0)
        throw SPhaseFileHDF5APIError("failed to close 'first_return_time' subgroup in surface group: "
                                     + string(group_name));
    size_t size = rho_dummy_buf.size();
    if ((phi_dummy_buf.size() != size) ||
        (mPhiT_dummy_buf.size() != size) ||
        (varPhiT_dummy_buf.size() != size) ||
        (mT_dummy_buf.size() != size) ||
        (varT_dummy_buf.size() != size) ||
        (mFRT_dummy_buf.size() != size) ||
        (varFRT_dummy_buf.size() != size))
        throw SPhaseFileHDF5APIError("the object: "
                                     + string(group_name) + string(" is faulty"));

    // store retrieved data in FRTData object
    FRTData& frtData = dataFile_ptr->createDataSet(g_name);
    frtData.x0[0] = rho_dummy_buf;
    frtData.x0[1] = phi_dummy_buf;
    frtData.mPhiT = mPhiT_dummy_buf;
    frtData.varPhiT = varPhiT_dummy_buf;
    frtData.mT = mT_dummy_buf;
    frtData.varT = varT_dummy_buf;
    frtData.mFRT = mFRT_dummy_buf;
    frtData.varFRT = varFRT_dummy_buf;

    herr_t surface_group_err = H5Gclose(surface_group_id);
    if (surface_group_err < 0)
        throw SPhaseFileHDF5APIError("failed to close group: " + string(group_name));

    return 0;

}

bool FRTDataFile::operator == (const FRTDataFile& other) const {
    bool is_equal = true;
    is_equal = is_equal && (this->configFilePath == other.configFilePath);
    is_equal = is_equal && (this->isoSurfaceFilePath == other.isoSurfaceFilePath);
    is_equal = is_equal && (this->modelName == other.modelName);
    for (auto t_it = this->data_ptr_Set.begin(), o_it = other.data_ptr_Set.begin();
            t_it != this->data_ptr_Set.end(), o_it != other.data_ptr_Set.end();
            ++t_it, ++o_it){
        is_equal = is_equal && (**t_it == **o_it);
    }
    return is_equal;
}