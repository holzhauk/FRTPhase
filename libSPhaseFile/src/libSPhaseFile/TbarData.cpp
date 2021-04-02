//
// Created by konstantin on 3/30/21.
//

#include "TbarData.h"

extern "C" herr_t surface_group_handler_tbardata(hid_t, const char*, const H5L_info_t*, void*);

bool TbarData::operator == (const TbarData& other) const {
    bool is_equal;
    is_equal = (this->isoSurfaceName == other.isoSurfaceName);
    is_equal = is_equal && (this->Tbar == other.Tbar);
    return is_equal;
}

TbarData& TbarData::operator = (const TbarData& other) {
    this->isoSurfaceName = other.isoSurfaceName;
    this->Tbar = other.Tbar;
    return *this;
}

TbarData& TbarDataFile::createDataSet(string isoSurfaceName) {
    TbarData dataSet(isoSurfaceName);
    dataList.push_back(dataSet);
    return dataList.back();
}

void TbarDataFile::write_body(H5::H5File& file) {

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

        H5::DataSpace dspace;

        H5::Group surface = file.createGroup(data.isoSurfaceName);

        // INITIAL POSITIONS
        H5::Group xinit_g = surface.createGroup("mean_period");
        hsize_t dim_size[] = {1};
        dspace = H5::DataSpace(1, dim_size);
        H5::DataSet tbar = xinit_g.createDataSet("Tbar",
                                                 H5::PredType::NATIVE_DOUBLE,
                                                 dspace);
        tbar.write(&data.Tbar, H5::PredType::NATIVE_DOUBLE);
        tbar.close();

        dspace.close();
        surface.close();

    }
}

void TbarDataFile::read_body(H5::H5File &file) {

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
                            surface_group_handler_tbardata,
                            this);
    if (idx < 0)
        throw SPhaseFileHDF5APIError("H5Literate encountered failure iterating through surface groups");
}

herr_t surface_group_handler_tbardata(hid_t group_id, const char* group_name,
                                     const H5L_info_t* group_info, void* opdata){

    string g_name(group_name);
    if ((g_name == "isosurface_file") || (g_name == "configuration_file"))
        return 0;

    auto dataFile_ptr = reinterpret_cast<TbarDataFile*> (opdata);

    // using C - interface from now on
    // open the surface group that are iterated through
    hid_t surface_group_id = H5Gopen(group_id, group_name, H5P_DEFAULT);
    if (surface_group_id < 0)
        throw SPhaseFileHDF5APIError(string("failed to open group: ") + string(group_name));

    ////////////////////////////////////////////////////////////////////
    // open mean_period subgroup to retrieve the mean period
    hid_t tbar_subgroup_id = H5Gopen(surface_group_id, "mean_period", H5P_DEFAULT);
    if (tbar_subgroup_id < 0)
        throw SPhaseFileHDF5APIError("failed to open 'mean_period' subgroup in surface group: "
                                     + string(group_name));

    // get mean period value
    hid_t dset_tbar_id = H5Dopen(tbar_subgroup_id, "Tbar", H5P_DEFAULT);
    if (dset_tbar_id < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset 'Tbar' in surface group: "
                                     + string(group_name));
    hid_t space_tbar_id = H5Dget_space(dset_tbar_id);
    if (space_tbar_id < 0)
        throw SPhaseFileHDF5APIError("failed to get dataspace from 'Tbar' in surface group: "
                                     + string(group_name));
    hsize_t size_dim_tbar[2];
    int ndims_tbar = H5Sget_simple_extent_dims(space_tbar_id, size_dim_tbar, NULL);
    if (ndims_tbar < 0)
        throw SPhaseFileHDF5APIError("failed to get the dspace's dims from 'Tbar' in surface group: "
                                     + string(group_name));

    if (size_dim_tbar[0] != 1)
        throw SPhaseFileHDF5APIError("Dataspace 'Tbar' in surface group: "
                                     + string(group_name) +
                                     " is not a scalar dataspace in contrast to expectations");

    double tbar_dummy_buf;
    herr_t status_tbar = H5Dread(dset_tbar_id,
                                H5T_NATIVE_DOUBLE,
                                H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                &tbar_dummy_buf);
    if (status_tbar < 0)
        throw SPhaseFileHDF5APIError("failed to read dataset 'Tbar' in surface group: "
                                     + string(group_name));
    H5Sclose(space_tbar_id);
    H5Dclose(dset_tbar_id);

    // store retrieved data in FRTData object
    TbarData& tbarData = dataFile_ptr->createDataSet(g_name);
    tbarData.Tbar = tbar_dummy_buf;

    herr_t surface_group_err = H5Gclose(surface_group_id);
    if (surface_group_err < 0)
        throw SPhaseFileHDF5APIError("failed to close group: " + string(group_name));

    return 0;

}

bool TbarDataFile::operator == (const TbarDataFile& other) const {
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

TbarDataFile& TbarDataFile::operator = (const TbarDataFile& other) {
    this->modelName = other.modelName;
    this->configFilePath = other.configFilePath;
    this->isoSurfaceFilePath = other.isoSurfaceFilePath;
    for (auto data: dataList){
        TbarData& new_data = this->createDataSet(data.isoSurfaceName);
        new_data = data;
    }
    return *this;
}