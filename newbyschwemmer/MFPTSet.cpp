//
// Created by konstantin on 4/8/20.
//

#include "MFPTSet.h"

MFPTSet::MFPTSet(std::string& modelname, fs::path& configurationfile, fs::path& isochronefile) {
    model_name = modelname;
    ConfigurationFile = configurationfile;
    IsochroneFile = isochronefile;
}

std::shared_ptr<MFPTs> MFPTSet::create(std::string isochrone_group_name) {
    std::shared_ptr<MFPTs> MFPTs_ptr = std::shared_ptr<MFPTs>(new MFPTs(isochrone_group_name));
    MFPTs_ptr_list.push_back(MFPTs_ptr);
    return MFPTs_ptr;
}

void MFPTSet::write_to_file(fs::path filepath) {
    try{
        /*
        * Turn off the auto-printing when failure occurs so that we can
        * handle the errors appropriately
        */
        Exception::dontPrint();

        H5File* file = new H5File(filepath, H5F_ACC_TRUNC);
        /*
         * Write Attributes
         */
        StrType vlst(0, H5T_VARIABLE);
        vlst.setCset(H5T_CSET_UTF8);
        DataSpace* scalarDataSpace = new DataSpace(H5S_SCALAR);
        Attribute* A;
        A = new Attribute(file->createAttribute("class", vlst, *scalarDataSpace));
        A->write(vlst, CLASS_TAG_MFPTS);
        delete A;
        A = new Attribute(file->createAttribute("format", vlst, *scalarDataSpace));
        A->write(vlst, FORMAT_TAG_MFPTS);
        delete A;
        A = new Attribute(file->createAttribute("version", vlst, *scalarDataSpace));
        A->write(vlst, VERSION_TAG_MFPTS);
        delete A;
        A = new Attribute(file->createAttribute("model", vlst, *scalarDataSpace));
        A->write(vlst, this->model_name);
        delete A;
        A = new Attribute(file->createAttribute("ConfigurationFile", vlst, *scalarDataSpace));
        A->write(vlst, this->ConfigurationFile);
        delete A;
        A = new Attribute(file->createAttribute("IsochroneFile", vlst, *scalarDataSpace));
        A->write(vlst, this->IsochroneFile);
        delete A;
        delete scalarDataSpace;

        for (auto MFPT_ptr: MFPTs_ptr_list) {
            Group* isochrone_group = new Group(file->createGroup(MFPT_ptr->get_isochrone_name()));

            DataSet* dSet;
            hsize_t dim_sizes[1];
            DataSpace* dataSpace;

            /*
             * Initial Positions
             */
            dim_sizes[0] = MFPT_ptr->get_sample_size();
            dataSpace = new DataSpace(1, dim_sizes);
            Group* initial_pos_group = new Group(isochrone_group->createGroup("InitialPositions"));
            dSet = new DataSet(initial_pos_group->createDataSet("Rho",
                    PredType::NATIVE_DOUBLE, *dataSpace));
            dSet->write(MFPT_ptr->get_initial_Rhos_buf_ptr(), PredType::NATIVE_DOUBLE);
            delete dSet;
            dSet = new DataSet(initial_pos_group->createDataSet("Phi",
                    PredType::NATIVE_DOUBLE, *dataSpace));
            dSet->write(MFPT_ptr->get_initial_Phis_buf_ptr(), PredType::NATIVE_DOUBLE);
            delete dSet;
            delete initial_pos_group;

            /*
             * MFPTs
             */
            dSet = new DataSet(isochrone_group->createDataSet("MFPT",
                    PredType::NATIVE_DOUBLE, *dataSpace));
            dSet->write(MFPT_ptr->get_MFPTs_buf_ptr(), PredType::NATIVE_DOUBLE);
            delete dSet;
            /*
             * VarFPTs
             */
            dSet = new DataSet(isochrone_group->createDataSet("VarFPT",
                    PredType::NATIVE_DOUBLE, *dataSpace));
            dSet->write(MFPT_ptr->get_VarFPTs_buf_ptr(), PredType::NATIVE_DOUBLE);
            delete dSet;

            /*
             * Mean Period
             */
            Group* mean_period_group = new Group(isochrone_group->createGroup("MeanPeriod"));
            // Tbars
            dSet = new DataSet(mean_period_group->createDataSet("Tbar",
                                                                PredType::NATIVE_DOUBLE, *dataSpace));
            dSet->write(MFPT_ptr->get_Tbars_buf_ptr(), PredType::NATIVE_DOUBLE);
            delete dSet;
            delete dataSpace;

            delete mean_period_group;

            delete isochrone_group;
        }
        delete file;
    }
        // catch failure caused by the H5File operations
    catch( FileIException error ){
        error.printErrorStack();
    }
        // catch failur cause bye the Group operations
    catch( GroupIException error ){
        error.printErrorStack();
    }
        // catch failure caused by the DataSet operations
    catch( DataSetIException error ){
        error.printErrorStack();
    }
        // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error ){
        error.printErrorStack();
    }
        // catch failure caused by the DataSpace operations
    catch( DataTypeIException error ){
        error.printErrorStack();
    }
}
