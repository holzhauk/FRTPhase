//
// Created by konstantin on 1/6/21.
//

#include "SPhaseFile.h"

void SPhaseFile::read(const fs::path& filePath){
    try{
        H5::Exception::dontPrint();

        // write the general header for all SPhaseFiles
        // containing attributes with information
        // about the file format and the objects it represents
        H5std_string FILENAME(filePath);
        H5::H5File file(FILENAME, H5F_ACC_RDONLY);

        // read in the file's attributes
        H5std_string attr_content;

        // check if the file specification meets the API's
        H5::Attribute attr = file.openAttribute("format");
        H5::DataType dtype = attr.getDataType();
        attr.read(dtype, attr_content);
        if (attr_content != SPhaseFile_FORMAT_ID)
            throw SPhaseFileWrongFormat();
        dtype.close();
        attr.close();

        attr = file.openAttribute("version");
        dtype = attr.getDataType();
        attr.read(dtype, attr_content);
        if (attr_content != SPhaseFile_VERSION_ID)
            throw SPhaseFileVersionConflict();
        dtype.close();
        attr.close();

        attr = file.openAttribute("class");
        dtype = attr.getDataType();
        attr.read(dtype, attr_content);
        if (attr_content != SPHASEFILE_CLASS_ID)
            throw SPhaseFileClassConflict();
        dtype.close();
        attr.close();

        // call pure virtual function that needs to be specified by
        // the children of this parent class containing detailed information
        // about the actual object that is represented in the file
        this->read_body(file);

        file.close();
    }
        // catch failure caused by the H5File operations
    catch( H5::FileIException error )
    {
        throw SPhaseFileHDF5APIError(error.getDetailMsg());
    }
        // catch failure caused by the DataSet operations
    catch( H5::DataSetIException error )
    {
        throw SPhaseFileHDF5APIError(error.getDetailMsg());
    }
        // catch failure caused by the DataSpace operations
    catch( H5::DataSpaceIException error )
    {
        throw SPhaseFileHDF5APIError(error.getDetailMsg());
    }
        // catch failure caused by the Attribute operations
    catch( H5::AttributeIException error )
    {
        throw SPhaseFileHDF5APIError(error.getDetailMsg());
    }
}

void SPhaseFile::write(const fs::path& filePath){
    try {
        H5::Exception::dontPrint();

        // write the general header for all SPhaseFiles
        // containing attributes with information
        // about the file format and the objects it represents
        H5std_string FILENAME(filePath);
        H5::H5File file(FILENAME, H5F_ACC_TRUNC);

        // define attribute types
        H5::StrType attr_str_type(0, H5T_VARIABLE);
        attr_str_type.setCset(H5T_CSET_UTF8);
        attr_str_type.setStrpad(H5T_STR_NULLTERM);
        H5::DataSpace attr_ds(H5S_SCALAR);

        // write file attributes
        H5::Attribute attr = file.createAttribute("class", attr_str_type, attr_ds);
        attr.write(attr_str_type, SPHASEFILE_CLASS_ID);
        attr.close();

        attr = file.createAttribute("format", attr_str_type, attr_ds);
        string str_dummy = string(SPhaseFile_FORMAT_ID);
        attr.write(attr_str_type, str_dummy);
        attr.close();

        attr = file.createAttribute("version", attr_str_type, attr_ds);
        str_dummy = string(SPhaseFile_VERSION_ID);
        attr.write(attr_str_type, str_dummy);
        attr.close();

        attr_ds.close();
        attr_str_type.close();

        // call pure virtual function that needs to be specified by
        // the children of this parent class containing detailed information
        // about the actual object that is represented in the file
        this->write_body(file);

        file.close();
    }
        // catch failure caused by the H5File operations
    catch( H5::FileIException error )
    {
        throw SPhaseFileHDF5APIError(error.getDetailMsg());
    }
        // catch failure caused by the DataSet operations
    catch( H5::DataSetIException error )
    {
        throw SPhaseFileHDF5APIError(error.getDetailMsg());
    }
        // catch failure caused by the DataSpace operations
    catch( H5::DataSpaceIException error )
    {
        throw SPhaseFileHDF5APIError(error.getDetailMsg());
    }
}
