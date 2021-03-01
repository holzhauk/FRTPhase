//
// Created by konstantin on 1/6/21.
//

#include "FRTData.h"
#include <iostream>

extern "C" herr_t surface_group_handler_frtdata(hid_t, const char*, const H5L_info_t*, void*);

bool FRTData::operator == (const FRTData& other) const {
    bool is_equal;
    is_equal = (this->isoSurfaceName == other.isoSurfaceName);
    is_equal = is_equal && (this->x0 == other.x0);
    is_equal = is_equal && (this->mFRT == other.mFRT);
    is_equal = is_equal && (this->varFRT == other.varFRT);
    return is_equal;
}

FRTData& FRTData::operator = (const FRTData& other) {
    this->isoSurfaceName = other.isoSurfaceName;
    this->x0 = other.x0;
    this->mFRT = other.mFRT;
    this->varFRT = other.varFRT;
    return *this;
}

FRTData& FRTDataFile::createDataSet(string isoSurfaceName) {
    FRTData dataSet(isoSurfaceName);
    dataList.push_back(dataSet);
    return dataList.back();
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
    for (auto data: dataList){

        hsize_t size = data.x0[0].size();
        if ((data.x0[0].size() != size) ||
            (data.mFRT.size() != size) ||
            (data.varFRT.size() != size)) throw SPhaseFileDimError();

        H5::DataSpace dspace;

        H5::Group surface = file.createGroup(data.isoSurfaceName);

        // INITIAL POSITIONS
        H5::Group xinit_g = surface.createGroup("initial_position");
        hsize_t dim_size[] = {data.x0[0].size()};
        dspace = H5::DataSpace(1, dim_size);
        H5::DataSet rhos = xinit_g.createDataSet("rho",
                                               H5::PredType::NATIVE_DOUBLE,
                                               dspace);
        rhos.write(data.x0[0].data(), H5::PredType::NATIVE_DOUBLE);
        rhos.close();
        H5::DataSet phis = xinit_g.createDataSet("phi",
                                               H5::PredType::NATIVE_DOUBLE,
                                               dspace);
        phis.write(data.x0[1].data(), H5::PredType::NATIVE_DOUBLE);
        phis.close();
        xinit_g.close();

        // FIRST RETURN TIMES FRT
        H5::Group frt_g = surface.createGroup("first_return_time");
        H5::DataSet mFRTs = frt_g.createDataSet("mFRT",
                                                H5::PredType::NATIVE_DOUBLE,
                                                dspace);
        mFRTs.write(data.mFRT.data(), H5::PredType::NATIVE_DOUBLE);
        mFRTs.close();
        H5::DataSet varFRTs = frt_g.createDataSet("varFRT",
                                                  H5::PredType::NATIVE_DOUBLE,
                                                  dspace);
        varFRTs.write(data.varFRT.data(), H5::PredType::NATIVE_DOUBLE);
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
        (mFRT_dummy_buf.size() != size) ||
        (varFRT_dummy_buf.size() != size))
        throw SPhaseFileHDF5APIError("the object: "
                                     + string(group_name) + string(" is faulty"));

    // store retrieved data in FRTData object
    FRTData& frtData = dataFile_ptr->createDataSet(g_name);
    frtData.x0[0] = rho_dummy_buf;
    frtData.x0[1] = phi_dummy_buf;
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
    for (auto t_it = this->dataList.begin(), o_it = other.dataList.begin();
         (t_it != this->dataList.end()) && (o_it != other.dataList.end());
            ++t_it, ++o_it){
        is_equal = is_equal && (*t_it == *o_it);
    }
    return is_equal;
}

FRTDataFile& FRTDataFile::operator = (const FRTDataFile& other) {
    this->modelName = other.modelName;
    this->configFilePath = other.configFilePath;
    this->isoSurfaceFilePath = other.isoSurfaceFilePath;
    for (auto data: dataList){
        FRTData& new_data = this->createDataSet(data.isoSurfaceName);
        new_data = data;
    }
    return *this;
}