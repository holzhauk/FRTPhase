//
// Created by konstantin on 1/13/21.
//

#include "MPIIsoSurfaceFile.h"

void MPI_Share(int world_rank, IsoSurfaceFile& isoSurfaceFile) {

    int root = 0;
    bool is_master = (world_rank == root);

    struct IsoSurfaceFileExt_t {
        unsigned int NoSurfaces;
        unsigned int modelNameSize;
        MPI_Datatype get_mpi_type() {
            // registration of mpi datatype for communication
            // of extension of data to be shared across the nodes
            MPI_Datatype MPI_IsoSurfaceFileExt_t;
            int blengths[2] = {1, 1};
            MPI_Datatype types[2] = {MPI_UNSIGNED, MPI_UNSIGNED};
            MPI_Aint offsets[2];
            offsets[0] = offsetof(IsoSurfaceFileExt_t, NoSurfaces);
            offsets[1] = offsetof(IsoSurfaceFileExt_t, modelNameSize);
            MPI_Type_create_struct(2, blengths, offsets, types, &MPI_IsoSurfaceFileExt_t);
            MPI_Type_commit(&MPI_IsoSurfaceFileExt_t);
            // return the composit type registered with mpi
            return MPI_IsoSurfaceFileExt_t;
        };
    };

    struct InterpolatedCurveExt_t {
        unsigned int nameSize;
        unsigned int NoNodes;
        unsigned int NoParameters;
        MPI_Datatype get_mpi_type() {
            // registration of mpi datatype for communication
            // of extension of data to be shared across the nodes
            MPI_Datatype MPI_InterpolatedCurveExt_t;
            int blengths[3] = {1, 1, 1};
            MPI_Datatype types[3] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED};
            MPI_Aint offsets[3];
            offsets[0] = offsetof(InterpolatedCurveExt_t, nameSize);
            offsets[1] = offsetof(InterpolatedCurveExt_t, NoNodes);
            offsets[2] = offsetof(InterpolatedCurveExt_t, NoParameters);
            MPI_Type_create_struct(3, blengths, offsets, types, &MPI_InterpolatedCurveExt_t);
            MPI_Type_commit(&MPI_InterpolatedCurveExt_t);
            // return the composit type registered with mpi
            return MPI_InterpolatedCurveExt_t;
        }
    };

    IsoSurfaceFileExt_t isoSurfaceFileExt;
    InterpolatedCurveExt_t interpolatedCurveExt;

    if (is_master) {
        string modelName = isoSurfaceFile.get_modelName();
        isoSurfaceFileExt.NoSurfaces = isoSurfaceFile.get_NoSurfaces();
        isoSurfaceFileExt.modelNameSize = modelName.size();
        MPI_Datatype MPI_IsoSurfaceFileExt_t = isoSurfaceFileExt.get_mpi_type();
        MPI_Bcast(&isoSurfaceFileExt, 1, MPI_IsoSurfaceFileExt_t, root, MPI_COMM_WORLD);
        MPI_Type_free(&MPI_IsoSurfaceFileExt_t);
        MPI_Bcast(modelName.data(), isoSurfaceFileExt.modelNameSize, MPI_CHAR, root, MPI_COMM_WORLD);

        for (auto curve: isoSurfaceFile){
            string curveName = curve.get_name();
            auto [rhos, phis] = curve.get_nodes();
            auto [keys, vals] = curve.get_parameters();
            interpolatedCurveExt.nameSize = curveName.size();
            interpolatedCurveExt.NoNodes = rhos.size();
            interpolatedCurveExt.NoParameters = keys.size();

            MPI_Datatype MPI_InterpolatedCurveExt_t = interpolatedCurveExt.get_mpi_type();
            MPI_Bcast(&interpolatedCurveExt, 1, MPI_InterpolatedCurveExt_t, root, MPI_COMM_WORLD);
            MPI_Type_free(&MPI_InterpolatedCurveExt_t);
            MPI_Bcast(curveName.data(), interpolatedCurveExt.nameSize, MPI_CHAR, root, MPI_COMM_WORLD);
            MPI_Bcast(rhos.data(), interpolatedCurveExt.NoNodes, MPI_DOUBLE, root, MPI_COMM_WORLD);
            MPI_Bcast(phis.data(), interpolatedCurveExt.NoNodes, MPI_DOUBLE, root, MPI_COMM_WORLD);

            // transmit the parameters
            vector<string>::iterator keyIt;
            vector<double>::iterator valIt;
            for (keyIt = keys.begin(), valIt = vals.begin();
                    keyIt != keys.end(), valIt != vals.end();
                    ++keyIt, ++valIt){
                size_t keySize = keyIt->size();
                MPI_Bcast(&keySize, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);
                MPI_Bcast(keyIt->data(), keySize, MPI_CHAR, root, MPI_COMM_WORLD);
                double val = *valIt;
                MPI_Bcast(&val, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
            }

        }
    } else {
        MPI_Datatype MPI_IsoSurfaceFileExt_t = isoSurfaceFileExt.get_mpi_type();
        MPI_Bcast(&isoSurfaceFileExt, 1, MPI_IsoSurfaceFileExt_t, root, MPI_COMM_WORLD);
        MPI_Type_free(&MPI_IsoSurfaceFileExt_t);
        char* str_buf = new char [isoSurfaceFileExt.modelNameSize + 1];
        MPI_Bcast(str_buf, isoSurfaceFileExt.modelNameSize, MPI_CHAR, root, MPI_COMM_WORLD);
        str_buf[isoSurfaceFileExt.modelNameSize] = '\0';
        string modelName(str_buf);
        isoSurfaceFile = IsoSurfaceFile(modelName);
        delete [] str_buf;

        for (int s = 0; s < isoSurfaceFileExt.NoSurfaces; s++){
            MPI_Datatype MPI_InterpolatedCurveExt_t = interpolatedCurveExt.get_mpi_type();
            MPI_Bcast(&interpolatedCurveExt, 1, MPI_InterpolatedCurveExt_t, root, MPI_COMM_WORLD);
            MPI_Type_free(&MPI_InterpolatedCurveExt_t);
            // receive curve name and create curve
            char* str_buf = new char [interpolatedCurveExt.nameSize + 1];
            MPI_Bcast(str_buf, interpolatedCurveExt.nameSize, MPI_CHAR, root, MPI_COMM_WORLD);
            str_buf[interpolatedCurveExt.nameSize] = '\0';
            string curveName(str_buf);
            delete [] str_buf;
            InterpolatedCurve& curve = isoSurfaceFile.createInterpolatedCurve(curveName);
            // receive notes and assign to the curve
            vector<double> rhos(interpolatedCurveExt.NoNodes);
            vector<double> phis(interpolatedCurveExt.NoNodes);
            MPI_Bcast(rhos.data(), interpolatedCurveExt.NoNodes, MPI_DOUBLE, root, MPI_COMM_WORLD);
            MPI_Bcast(phis.data(), interpolatedCurveExt.NoNodes, MPI_DOUBLE, root, MPI_COMM_WORLD);
            curve.set_nodes(rhos, phis);

            // receive parameter set and add to the curve
            for (int p = 0; p < interpolatedCurveExt.NoParameters; p++){
                unsigned int keySize;
                MPI_Bcast(&keySize, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);
                char* buf = new char [keySize + 1];
                MPI_Bcast(buf, keySize, MPI_CHAR, root, MPI_COMM_WORLD);
                buf[keySize] = '\0';
                string key(buf);
                delete [] buf;
                double val;
                MPI_Bcast(&val, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
                curve.add_parameter(key, val);
            }
        }
    }

}