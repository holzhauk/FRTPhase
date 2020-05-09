//
// Created by konstantin on 4/9/20.
//

#include <string>

#include "H5Cpp.h"
#include "TrajectorySet.h"

TrajectorySet::TrajectorySet(std::map<std::string, double>& p_set, std::string& model_name_tag) {
    parameters = p_set;
    model_name = model_name_tag;
}

TrajectorySet::~TrajectorySet() {
    for (auto Trajectory_ptr : Trajectory_ptr_list){
        delete Trajectory_ptr;
    }
}

void TrajectorySet::push_back(Trajectory* trajectory) {
    Trajectory_ptr_list.push_back(trajectory);
}

void TrajectorySet::write_to_file(fs::path filepath) {
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
        A->write(vlst, CLASS_TAG_TRAJECTORY_SET);
        delete A;
        A = new Attribute(file->createAttribute("format", vlst, *scalarDataSpace));
        A->write(vlst, FORMAT_TAG_TRAJECTORY_SET);
        delete A;
        A = new Attribute(file->createAttribute("version", vlst, *scalarDataSpace));
        A->write(vlst, VERSION_TAG_TRAJECTORY_SET);
        delete A;
        A = new Attribute(file->createAttribute("model", vlst, *scalarDataSpace));
        A->write(vlst, this->model_name);
        delete A;
        delete scalarDataSpace;

        DataSet* dSet;
        DataSpace* dataSpace;
        hsize_t dim_size[1] = {1};
        /*
         * Write Parameterlist
         */
        Group* parameter_group = new Group(file->createGroup("Parameters"));
        dataSpace = new DataSpace(1, dim_size);
        for (auto P: parameters) {
            dSet = new DataSet(parameter_group->createDataSet(P.first, PredType::NATIVE_DOUBLE, *dataSpace));
            dSet->write(&P.second, PredType::NATIVE_DOUBLE);
            delete dSet;
        }
        delete dataSpace;
        delete parameter_group;

        Group* trajectories_group = new Group(file->createGroup("Trajectories"));
        int i = 0; // Counter

        for (auto Trajectory_ptr: Trajectory_ptr_list) {

            Group* trajectory_group = new Group(trajectories_group->createGroup("T" + std::to_string(i)));
            /*
             * State Vectors
             */
            dim_size[0] = Trajectory_ptr->get_Rhos().size();
            dataSpace = new DataSpace(1, dim_size);
            dSet = new DataSet(trajectory_group->createDataSet("Rho",
                                                              PredType::NATIVE_DOUBLE, *dataSpace));
            dSet->write(Trajectory_ptr->get_Rho_buf_ptr(), PredType::NATIVE_DOUBLE);
            delete dSet;
            dSet = new DataSet(trajectory_group->createDataSet("Phi",
                                                              PredType::NATIVE_DOUBLE, *dataSpace));
            dSet->write(Trajectory_ptr->get_Phi_buf_ptr(), PredType::NATIVE_DOUBLE);
            delete dSet;
            dSet = new DataSet(trajectory_group->createDataSet("t",
                                                               PredType::NATIVE_DOUBLE, *dataSpace));
            dSet->write(Trajectory_ptr->get_t_buf_ptr(), PredType::NATIVE_DOUBLE);
            delete dSet;
            delete dataSpace;

            i++;

            delete trajectory_group;
        }
        delete trajectories_group;
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