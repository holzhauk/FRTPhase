//
// Created by konstantin on 1/3/21.
//
#include <cmath>
#include <utility>
#include <algorithm>
#include <H5Cpp.h>

#include "IsoSurface.h"

extern "C" {
    herr_t surface_group_handler_isosurface(hid_t, const char*, const H5L_info_t*, void*);
    herr_t parameter_dset_handler(hid_t, const char*, const H5L_info_t*, void*);
}

string InterpolatedCurve::get_name(){
    return name;
}

tuple<double, double> InterpolatedCurve::get_extensions() {
    if (nodes[0].size() == 0)
        return {0.0, 0.0};

    double rho_min, rho_max;
    rho_min = nodes[0][0];
    rho_max = nodes[0][nodes[0].size() - 1];
    return tie(rho_min, rho_max);
}

ParameterSet InterpolatedCurve::get_parameterSet() {
    return ParameterSet(parameterSet);
}

void InterpolatedCurve::add_parameter(string pName, double pVal){
    parameterSet[pName] = pVal;
}

void InterpolatedCurve::add_node(array<double, 2>& x) {
    double Rho = x[0];
    double Phi = x[1];
    nodes[0].push_back(Rho);
    nodes[1].push_back(Phi);
    this->sortNodes();
}

void InterpolatedCurve::set_nodes(vector<double>& rho, vector<double>& phi) {
    nodes[0] = vector<double>(rho);
    nodes[1] = vector<double>(phi);
    this->sortNodes();
}

tuple<vector<double>, vector<double>> InterpolatedCurve::get_nodes() {
    return {vector<double>(nodes[0]), vector<double>(nodes[1])};
}

tuple<vector<string>, vector<double>> InterpolatedCurve::get_parameters() {
    vector<string> keys;
    vector<double> vals;
    for (auto& p: parameterSet){
        keys.push_back(p.first);
        vals.push_back(p.second);
    }
    return {keys, vals};
}

bool InterpolatedCurve::operator==(const InterpolatedCurve &other) const {
    bool is_equal = true;
    is_equal = is_equal && (this->name == other.name);
    is_equal = is_equal && (this->parameterSet == other.parameterSet);
    is_equal = is_equal && (this->nodes == other.nodes);
    return is_equal;
}


bool IsoSurfaceFile::operator==(const IsoSurfaceFile& other) const {
    bool is_equal = true;
    if (this->curve_ptr_Set.size() != other.curve_ptr_Set.size())
        return false;
    is_equal = is_equal && (this->modelName == other.modelName);

    for (auto t_it = this->curve_ptr_Set.begin(), o_it = other.curve_ptr_Set.begin();
            t_it != this->curve_ptr_Set.end(), o_it != other.curve_ptr_Set.end();
            ++t_it, ++o_it){
        is_equal = is_equal && (**(t_it) == **(o_it));
    }
    return is_equal;
}

int IsoSurfaceFile::get_NoSurfaces() {
    return curve_ptr_Set.size();
}

void IsoSurfaceFile::read_body(H5::H5File& file) {

    H5std_string attr_content;

    H5::Attribute attr = file.openAttribute("model");
    H5::DataType dtype = attr.getDataType();
    attr.read(dtype, attr_content);
    if (attr_content != modelName)
        throw SPhaseFileModelConflict();
    dtype.close();
    attr.close();

    herr_t idx = H5Literate(file.getId(),
                            H5_INDEX_NAME,
                            H5_ITER_INC,
                            NULL,
                            surface_group_handler_isosurface,
                            this);
    if (idx < 0)
        throw SPhaseFileHDF5APIError("H5Literate encountered failure iterating through surface groups");
}

herr_t surface_group_handler_isosurface(hid_t group_id, const char* group_name,
                             const H5L_info_t* group_info, void* opdata){

    auto isoSurfaceFile_ptr = reinterpret_cast<IsoSurfaceFile*>(opdata);
    InterpolatedCurve& interpolatedCurve = isoSurfaceFile_ptr->createInterpolatedCurve(group_name);

    // using C - interface from now on
    // open the surface group that are iterated through
    hid_t surface_group_id = H5Gopen(group_id, group_name, H5P_DEFAULT);
    if (surface_group_id < 0)
        throw SPhaseFileHDF5APIError(string("failed to open group: ") + string(group_name));

    // open Parameter subgroup to retrieve the parameters of the mode
    // that lead to the curve representing the respective isosurface
    hid_t param_subgroup_id = H5Gopen(surface_group_id, "parameters", H5P_DEFAULT);
    if (param_subgroup_id < 0)
        throw SPhaseFileHDF5APIError("failed to open subgroup 'parameters' in: "
                                    + string(group_name));
    // iterate over the scalar parameter datasets
    herr_t idy = H5Literate(param_subgroup_id,
                            H5_INDEX_NAME,
                            H5_ITER_INC,
                            NULL,
                            parameter_dset_handler,
                            &interpolatedCurve);
    if (idy < 0)
        throw SPhaseFileHDF5APIError("H5Literate encountered failure iterating"
                                        + string("through parameter groups in surface group: ")
                                        + string(group_name));

    herr_t param_subgroup_err = H5Gclose(param_subgroup_id);
    if (param_subgroup_err < 0)
        throw SPhaseFileHDF5APIError("failed to close subgroup 'parameters' in: "
                                    + string(group_name));

    // open curve subgroup to retrieve the node's coordinates
    hid_t curve_subgroup_id = H5Gopen(surface_group_id, "curve", H5P_DEFAULT);
    if (curve_subgroup_id < 0)
        throw SPhaseFileHDF5APIError("failed to open 'curve' subgroup in surface group: "
            + string(group_name));

    // get radial coordinates
    hid_t dset_rho_id = H5Dopen(curve_subgroup_id, "rho", H5P_DEFAULT);
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
    hid_t dset_phi_id = H5Dopen(curve_subgroup_id, "phi", H5P_DEFAULT);
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

    if (rho_dummy_buf.size() != phi_dummy_buf.size())
        throw SPhaseFileHDF5APIError("rho and phi have different sizes in surface group: "
                                    + string(group_name));

    interpolatedCurve.set_nodes(rho_dummy_buf, phi_dummy_buf);

    herr_t curve_subgroup_err = H5Gclose(curve_subgroup_id);
    if (curve_subgroup_err < 0)
        throw SPhaseFileHDF5APIError("failed to close 'curve' subgroup in surface group: "
            + string(group_name));

    herr_t surface_group_err = H5Gclose(surface_group_id);
    if (surface_group_err < 0)
        throw SPhaseFileHDF5APIError("failed to close group: " + string(group_name));

    return 0;
}

herr_t parameter_dset_handler(hid_t p_id, const char* p_name, const H5L_info_t* p_info, void* opdata){

    InterpolatedCurve* curve_ptr = reinterpret_cast<InterpolatedCurve*> (opdata);

    hid_t p_dset = H5Dopen(p_id, p_name, H5P_DEFAULT);
    if (p_dset < 0)
        throw SPhaseFileHDF5APIError("failed to open dataset: " + string(p_name));
    hid_t p_dspace = H5Dget_space(p_dset);
    if (p_dspace < 0)
        throw SPhaseFileHDF5APIError("failed to get data space info of dataset: "
                                    + string(p_name));

    hsize_t size_dim0[2];
    int ndims = H5Sget_simple_extent_dims(p_dspace, size_dim0, NULL);
    if (ndims < 0)
        throw SPhaseFileHDF5APIError("failed to get extensions of dataspace of dataset: "
                                        + string(p_name));
    double double_buf;
    herr_t idx = H5Dread(p_dset,
                         H5T_NATIVE_DOUBLE,
                         H5S_ALL, H5S_ALL, H5P_DEFAULT,
                         &double_buf);
    if (idx < 0)
        throw SPhaseFileHDF5APIError("failed to read dataset: "
                                        + string(p_name));
    curve_ptr->add_parameter(string(p_name), double_buf);

    herr_t p_dspace_err = H5Sclose(p_dspace);
    if (p_dspace_err < 0)
       throw SPhaseFileHDF5APIError("failed to close data space of dataset: "
                                    + string(p_name));
    herr_t p_dset_err = H5Dclose(p_dset);
    if (p_dset_err < 0)
        throw SPhaseFileHDF5APIError("failed to close dataset: "
                                    + string(p_name));

    return 0;
}

void IsoSurfaceFile::write_body(H5::H5File& file){

    // define attribute types
    H5::StrType attr_str_type(0, H5T_VARIABLE);
    attr_str_type.setCset(H5T_CSET_UTF8);
    attr_str_type.setStrpad(H5T_STR_NULLTERM);
    H5::DataSpace attr_ds(H5S_SCALAR);
    // write model name attribute
    H5::Attribute attr = file.createAttribute("model", attr_str_type, attr_ds);
    attr.write(attr_str_type, modelName);
    attr.close();
    attr_ds.close();
    attr_str_type.close();

    // write the individual curves representing
    // the isosurfaces to the file
    for (auto& curve_ptr: curve_ptr_Set){

        if (curve_ptr->nodes[0].size() != curve_ptr->nodes[1].size())
            throw SPhaseFileDimError();

        H5::DataSpace dspace;

        H5::Group surface = file.createGroup(curve_ptr->name);

        H5::Group curve = surface.createGroup("curve");
        hsize_t dim_size[] = {curve_ptr->nodes[0].size()};
        dspace = H5::DataSpace(1, dim_size);
        H5::DataSet rhos = curve.createDataSet("rho",
                                               H5::PredType::NATIVE_DOUBLE,
                                               dspace);
        rhos.write(curve_ptr->nodes[0].data(), H5::PredType::NATIVE_DOUBLE);
        rhos.close();
        H5::DataSet phis = curve.createDataSet("phi",
                                               H5::PredType::NATIVE_DOUBLE,
                                               dspace);
        phis.write(curve_ptr->nodes[1].data(), H5::PredType::NATIVE_DOUBLE);
        phis.close();
        dspace.close();
        curve.close();

        H5::Group parameters = surface.createGroup("parameters");
        for (auto p: curve_ptr->parameterSet){
            hsize_t dim_size[] = {1};
            dspace = H5::DataSpace(1, dim_size);
            H5::DataSet param = parameters.createDataSet(p.first,
                                                         H5::PredType::NATIVE_DOUBLE,
                                                         dspace);
            double buf = p.second;
            param.write(&buf, H5::PredType::NATIVE_DOUBLE);
            param.close();
            dspace.close();
        }
        parameters.close();

        surface.close();

    }
}

InterpolatedCurve& IsoSurfaceFile::createInterpolatedCurve(string name) {
    unique_ptr<InterpolatedCurve> curve_ptr(new InterpolatedCurve(name));
    curve_ptr_Set.push_back(move(curve_ptr));
    return *curve_ptr_Set.back();
}

bool InterpolatedCurve::is_first_return_event(array<double, 2>& x, bool pos_sense_of_rotation) {

    // assuming rho_i of isosurface are equally spaced
    // and ordered along the rho-axis
    // -> we get the index as follows
    // this could also be solved more elegantly by using the STL
    auto& [rho, phi] = x;
    size_t dof = nodes[0].size();
    double& rho_min = nodes[0][0];
    double& rho_max = nodes[0][dof - 1];

    if ((rho > rho_max) || (rho < rho_min))
        throw SPhaseFileOutOfDomain();

    int k = floor(((double)dof - 1.0) * (rho - rho_min) / (rho_max - rho_min));

    if (pos_sense_of_rotation)
        return (phi >= nodes[1][k] + 2 * M_PI);
    else
        return (phi <= nodes[1][k] - 2 * M_PI);

}

void InterpolatedCurve::sortNodes(){
    if (nodes[0].size() != nodes[1].size())
        throw SPhaseFileDimError();

    // rewrite array of vectors into vector of pairs
    vector<pair<double, double>> dummy;
    for(auto rho_it = nodes[0].begin(), phi_it = nodes[1].begin();
            rho_it != nodes[0].end(), phi_it != nodes[1].end();
            ++rho_it, ++phi_it){
        dummy.push_back({*rho_it, *phi_it});
    }
    sort(dummy.begin(), dummy.end()); // sort for non-descending order in rho
    nodes[0] = vector<double>();
    nodes[1] = vector<double>();
    for (auto& x: dummy){
        nodes[0].push_back(x.first);
        nodes[1].push_back(x.second);
    }
}

IsoSurfaceFile::IsoSurfaceFileIt IsoSurfaceFile::begin() {
    return IsoSurfaceFileIt(curve_ptr_Set.begin());
}

IsoSurfaceFile::IsoSurfaceFileIt IsoSurfaceFile::end() {
    return IsoSurfaceFileIt(curve_ptr_Set.end());
}

bool IsoSurfaceFile::IsoSurfaceFileIt::operator!=(IsoSurfaceFileIt& other) {
    return (this->ptrListIt != other.ptrListIt);
}

InterpolatedCurve& IsoSurfaceFile::IsoSurfaceFileIt::operator*() const {
    return **ptrListIt;
}

InterpolatedCurve& IsoSurfaceFile::IsoSurfaceFileIt::operator++() {
    return **(ptrListIt++);
}

string IsoSurfaceFile::get_modelName() {
    return string(modelName);
}

InterpolatedCurve& InterpolatedCurve::operator = (const InterpolatedCurve& other) {
    this->name = other.name;
    this->parameterSet = other.parameterSet;
    this->nodes = other.nodes;
    return *this;
}

IsoSurfaceFile& IsoSurfaceFile::operator = (const IsoSurfaceFile& other) {
    this->modelName = other.modelName;
    for (auto& curve_ptr: other.curve_ptr_Set) {
        InterpolatedCurve& new_curve = this->createInterpolatedCurve(curve_ptr->name);
        new_curve = *curve_ptr;
    }
    return *this;
}