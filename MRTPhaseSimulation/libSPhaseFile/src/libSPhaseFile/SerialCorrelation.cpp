//
// Created by konstantin on 2/7/21.
//

#include <vector>

#include "SerialCorrelation.h"

using namespace std;

extern "C" herr_t surface_group_handler_corr(hid_t, const char*, const H5L_info_t*, void*);

string IsoSurfaceCorr::get_key() {
    return string(key);
}

IsoSurfaceCorr& IsoSurfaceCorr::operator=(const IsoSurfaceCorr& other) {
    this->key = other.key;
    this->N = other.N;
    this->offset = other.offset;
    this->rho_k = other.rho_k;
    return *this;
}

bool IsoSurfaceCorr::operator== (const IsoSurfaceCorr& other) const {
    bool is_equal;
    is_equal = (this->key == other.key);
    is_equal = is_equal && (this->N == other.N);
    is_equal = is_equal && (this->offset == other.offset);
    is_equal = is_equal && (this->rho_k == other.rho_k);
    return is_equal;
}

void SerialCorrFile::read_body(H5::H5File &file) {

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
                            surface_group_handler_corr,
                            this);
    if (idx < 0)
        throw SPhaseFileHDF5APIError("H5Literate encountered failure iterating through surface groups");

}

herr_t surface_group_handler_corr(hid_t group_id, const char* group_name,
                                        const H5L_info_t* group_info, void* opdata){

    string g_name(group_name);
    if ((g_name == "isosurface_file") || (g_name == "configuration_file"))
        return 0;

    auto SerialCorrFile_ptr = reinterpret_cast<SerialCorrFile*>(opdata);
    IsoSurfaceCorr& isoSurfaceCorr = SerialCorrFile_ptr->create_isoSurfaceCorr(group_name);

    // using C - interface from now on
    // open the surface group that are iterated through
    hid_t surface_group_id = H5Gopen(group_id, group_name, H5P_DEFAULT);
    if (surface_group_id < 0)
        throw SPhaseFileHDF5APIError(string("failed to open group: ") + string(group_name));

    // read interval count
    hid_t dSet_N = H5Dopen(surface_group_id, "total_number_of_intervals", H5P_DEFAULT);
    if (dSet_N < 0)
        throw SPhaseFileHDF5APIError("failed to open 'total_number_of_intervals' in surface group:"
                                    + string(group_name));

    size_t scalar_dummy;
    herr_t idx = H5Dread(dSet_N, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &scalar_dummy);
    if (idx < 0)
        throw SPhaseFileHDF5APIError("failed to read 'total_number_of_intervals' in surface group: "
                                    + string(group_name));
    isoSurfaceCorr.N = scalar_dummy;

    herr_t dSet_N_err = H5Dclose(dSet_N);
    if (dSet_N_err < 0)
        throw SPhaseFileHDF5APIError("failed to close dataset 'total_number_of_intervals' in surface group: "
                                    + string(group_name));

    // read offset
    hid_t dSet_offset = H5Dopen(surface_group_id, "stationary_offset", H5P_DEFAULT);
    if (dSet_offset < 0)
        throw SPhaseFileHDF5APIError("failed to open 'stationary_offset' in surface group:"
                                     + string(group_name));

    herr_t idy = H5Dread(dSet_offset, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &scalar_dummy);
    if (idy < 0)
        throw SPhaseFileHDF5APIError("failed to read 'stationary_offset' in surface group: "
                                     + string(group_name));
    isoSurfaceCorr.offset = scalar_dummy;

    herr_t dSet_offset_err = H5Dclose(dSet_offset);
    if (dSet_offset_err < 0)
        throw SPhaseFileHDF5APIError("failed to close dataset 'stationary_offset' in surface group: "
                                     + string(group_name));

    // read serial correlation coefficients up to lag k
    hid_t dSet_rho_k = H5Dopen(surface_group_id, "rho_k", H5P_DEFAULT);
    if (dSet_rho_k < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset 'rho_k' in surface group: "
                                     + string(group_name));

    hid_t dSpace_rho_k = H5Dget_space(dSet_rho_k);
    if (dSpace_rho_k < 0)
        throw SPhaseFileHDF5APIError("failed to retrieve data space info about 'rho_k' in surface group: "
                                     + string(group_name));

    hsize_t dims[2];
    int ndims = H5Sget_simple_extent_dims(dSpace_rho_k, dims, NULL);
    if (ndims < 0)
        throw SPhaseFileHDF5APIError("faild to retrieve dspace's dims of 'rho_k' in surface group: "
                                     + string(group_name));

    isoSurfaceCorr.rho_k = vector<double> (dims[0]);
    herr_t idz = H5Dread(dSet_rho_k,
                         H5T_NATIVE_DOUBLE,
                         H5S_ALL,
                         H5S_ALL,
                         H5P_DEFAULT,
                         isoSurfaceCorr.rho_k.data());
    if (idz < 0)
        throw SPhaseFileHDF5APIError("failed to read data from 'rho_k' in surface group: "
                                    + string(group_name));

    herr_t dSpace_rho_k_err = H5Sclose(dSpace_rho_k);
    if (dSpace_rho_k_err < 0)
        throw SPhaseFileHDF5APIError("failed to close dataspace of 'rho_k' in surface group: "
                                    + string(group_name));

    herr_t dSet_rho_k_err = H5Dclose(dSet_rho_k);
    if (dSet_rho_k_err < 0)
        throw SPhaseFileHDF5APIError("failed to close dataset 'rho_k' in surface group: "
                                     + string(group_name));

    herr_t surface_group_err = H5Gclose(surface_group_id);
    if (surface_group_err < 0)
        throw SPhaseFileHDF5APIError("failed to close group: " + string(group_name));

    return 0;
}

void SerialCorrFile::write_body(H5::H5File &file) {

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

    for (auto isoSurfaceCorr: this->corrList) {

        H5::Group surface = file.createGroup(isoSurfaceCorr.get_key());
        H5::DataSet dset;

        const int rank = 1;
        // write down interval count and offset
        hsize_t scalar_dim[] = {1};
        H5::DataSpace scalar_dspace (rank, scalar_dim);
        size_t dummy;

        dset = surface.createDataSet("total_number_of_intervals", H5::PredType::NATIVE_UINT, scalar_dspace);
        dummy = isoSurfaceCorr.N;
        dset.write(&dummy, H5::PredType::NATIVE_HSIZE);
        dset.close();

        dset = surface.createDataSet("stationary_offset", H5::PredType::NATIVE_UINT, scalar_dspace);
        dummy = isoSurfaceCorr.offset;
        dset.write(&dummy, H5::PredType::NATIVE_HSIZE);
        dset.close();

        scalar_dspace.close();

        // write down serial correlation coefficients up to lag k
        hsize_t dims[] = {isoSurfaceCorr.rho_k.size()};
        H5::DataSpace dspace (rank, dims);
        dset = surface.createDataSet("rho_k", H5::PredType::NATIVE_DOUBLE, dspace);
        dset.write(isoSurfaceCorr.rho_k.data(), H5::PredType::NATIVE_DOUBLE);
        dset.close();
        dspace.close();

        surface.close();

    }

}

IsoSurfaceCorr& SerialCorrFile::create_isoSurfaceCorr(const string &isoSurfaceName) {
    IsoSurfaceCorr isoSurfaceCorr (isoSurfaceName);
    corrList.push_back(isoSurfaceCorr);
    return corrList.back();
}

SerialCorrFile& SerialCorrFile::operator = (const SerialCorrFile& other) {
    this->modelName = other.modelName;
    this->isoSurfaceFilePath = other.isoSurfaceFilePath;
    this->configFilePath = other.configFilePath;
    for (auto corr: other.corrList) {
        IsoSurfaceCorr corr_cpy (other.modelName);
        corr_cpy = corr;
        this->corrList.push_back(corr_cpy);
    }
    return *this;
}

bool SerialCorrFile::operator== (const SerialCorrFile& other) const{
    bool is_equal;
    is_equal = (this->modelName == other.modelName);
    is_equal = is_equal && (this->configFilePath == other.configFilePath);
    is_equal = is_equal && (this->isoSurfaceFilePath == other.isoSurfaceFilePath);

    auto LHScorrIt = this->corrList.begin();
    auto RHScorrIt = other.corrList.begin();
    for (LHScorrIt = this->corrList.begin(), RHScorrIt = other.corrList.begin();
         (LHScorrIt != this->corrList.end()) && (RHScorrIt != other.corrList.end());
        ++LHScorrIt, ++RHScorrIt){
        is_equal = is_equal && (*LHScorrIt == *RHScorrIt);
    }
    return is_equal;
}