import h5py
import os
from pathlib import Path
import json
import numpy as np
from .Exceptions import *


class SPhaseFile():
    FORMAT = "SPhaseFile"
    VERSION = "0.0.1"
    CLASSTAG = ""

    def __init__(self, model_name: str, dtype: np.dtype = np.double):
        self.MODELNAME = model_name
        self.dtype = dtype

    def __read_body__(self, file):
        pass

    def __write_body__(self, file):
        pass

    def read(self, file_path):
        if not os.path.isfile(file_path):
            raise SPhaseFileIOError(file_path, "the file does not exist")

        file = h5py.File(file_path, "r+")

        if file.attrs.get("format") != self.FORMAT:
            raise SPhaseFileWrongFormat(file_path, "is not" + self.FORMAT)
        if file.attrs.get("version") != self.VERSION:
            raise SPhaseFileVersionConflict(file_path, \
                                            "file version: " + file.attrs.get("version") + \
                                            " is not version " + self.VERSION)
        if file.attrs.get("class") != self.CLASSTAG:
            raise SPhaseFileClassConflict(file_path, \
                                          "file is of class: " + file.attrs.get("class"))
        if file.attrs.get("model") != self.MODELNAME:
            raise SPhaseFileModelConflict(file_path, \
                                          "data in file is associated to a different model than the one specified")

        self.__read_body__(file)
        file.close()

    def write(self, file_path):
        file = h5py.File(file_path, "w")
        file.attrs.create("format", self.FORMAT)
        file.attrs.create("version", self.VERSION)
        file.attrs.create("class", self.CLASSTAG)
        file.attrs.create("model", self.MODELNAME)
        self.__write_body__(file)
        file.close()


class Curve():

    def __init__(self, key, dtype=np.double):
        self.dtype = dtype
        self.key = key
        self.parameterSet = {}
        self.omegaBar = 0.0
        self.rho = np.array([], dtype=self.dtype)
        self.phi = np.array([], dtype=self.dtype)

    def __eq__(self, other):
        if not isinstance(other, Curve):
            return NotImplemented
        are_equal = (self.dtype == other.dtype)
        are_equal = are_equal and (self.key == other.key)
        are_equal = are_equal and (self.omegaBar == other.omegaBar)
        are_equal = are_equal and (self.parameterSet == other.parameterSet)
        are_equal = are_equal and np.array_equal(self.rho, other.rho)
        are_equal = are_equal and np.array_equal(self.phi, other.phi)
        return are_equal


class IsoSurfaceFile(SPhaseFile):
    CLASSTAG = "IsoSurfaceFile"

    def __init__(self, model_name):
        super(IsoSurfaceFile, self).__init__(model_name)
        self.curveSet = []

    def __write_body__(self, file):
        for curve in self.curveSet:
            IsoSurfaceGroup = file.create_group(curve.key)
            Parameter = IsoSurfaceGroup.create_group("parameters")
            for key in curve.parameterSet:
                Parameter.create_dataset(key, shape=(1,), \
                                         dtype=self.dtype, \
                                         data=curve.parameterSet[key])

            IsoSurfaceGroup.create_dataset("omega_bar", shape=(1,), \
                                           dtype=self.dtype, \
                                           data=curve.omegaBar)

            CurveGroup = IsoSurfaceGroup.create_group("curve")
            CurveGroup.create_dataset("rho", shape=curve.rho.shape, \
                                      dtype=self.dtype, \
                                      data=curve.rho)
            CurveGroup.create_dataset("phi", shape=curve.phi.shape, \
                                      dtype=self.dtype, \
                                      data=curve.phi)

    def __read_body__(self, file):
        for group_key in file:
            IsoSurfaceGroup = file[group_key]
            curve = self.createCurve(group_key)

            scalar_dummy = np.zeros((1,), dtype=self.dtype)
            Parameter = IsoSurfaceGroup["parameters"]
            for key in Parameter:
                Parameter[key].read_direct(scalar_dummy)
                curve.parameterSet[key] = scalar_dummy.copy()

            dSet_omega_bar = IsoSurfaceGroup.get("omega_bar")
            with dSet_omega_bar.astype(self.dtype):
                dSet_omega_bar.read_direct(scalar_dummy)
                curve.omegaBar = scalar_dummy.copy()

            CurveGroup = IsoSurfaceGroup["curve"]
            dSet_rho = CurveGroup["rho"]
            with dSet_rho.astype(self.dtype):
                curve.rho = dSet_rho[:]
            dSet_phi = CurveGroup["phi"]
            with dSet_phi.astype(self.dtype):
                curve.phi = dSet_phi[:]

    def __iter__(self):
        self.IsoSurfaceIterator = iter(self.curveSet)
        return self.IsoSurfaceIterator

    def __next__(self):
        return next(self.IsoSurfaceIterator)

    def __eq__(self, other):
        if not isinstance(other, IsoSurfaceFile):
            return NotImplemented

        are_equal = True
        for (c1, c2) in zip(self.curveSet, other.curveSet):
            are_equal = are_equal and (c1 == c2)
        return are_equal

    def createCurve(self, key):
        newCurve = Curve(key)
        self.curveSet.append(newCurve)
        return newCurve


class FRTData():
    def __init__(self, key, dtype=np.double):
        self.dtype = dtype
        self.key = key
        self.rho0 = np.array([], dtype=self.dtype)
        self.phi0 = np.array([], dtype=self.dtype)
        self.mFRT = np.array([], dtype=self.dtype)
        self.varFRT = np.array([], dtype=self.dtype)

    def __eq__(self, other):
        if not isinstance(other, FRTData):
            return NotImplemented

        are_equal = (self.dtype == other.dtype)
        are_equal = are_equal and (self.key == other.key)
        are_equal = are_equal and np.array_equal(self.rho0, other.rho0)
        are_equal = are_equal and np.array_equal(self.phi0, other.phi0)
        are_equal = are_equal and np.array_equal(self.mFRT, other.mFRT)
        are_equal = are_equal and np.array_equal(self.varFRT, other.varFRT)
        return are_equal


class FRTDataFile(SPhaseFile):
    CLASSTAG = "FRTDataFile"

    def __init__(self, model_name):
        super(FRTDataFile, self).__init__(model_name)
        self.isoSurfaceFilePath = ""
        self.configFilePath = ""
        self.dataSet = []

    def __write_body__(self, file):
        str_dtype = h5py.string_dtype(encoding='utf-8')
        file.create_dataset("isosurface_file", \
                            dtype=str_dtype, \
                            data=self.isoSurfaceFilePath.encode('utf-8'))
        file.create_dataset("configuration_file", \
                            dtype=str_dtype, \
                            data=self.configFilePath.encode('utf-8'))
        for frtData in self.dataSet:
            IsoSurfaceGroup = file.create_group(frtData.key)

            XinitGroup = IsoSurfaceGroup.create_group("initial_position")
            XinitGroup.create_dataset("rho", shape=frtData.rho0.shape, \
                                      dtype=self.dtype, \
                                      data=frtData.rho0)
            XinitGroup.create_dataset("phi", shape=frtData.phi0.shape, \
                                      dtype=self.dtype, \
                                      data=frtData.phi0)

            FRTGroup = IsoSurfaceGroup.create_group("first_return_time")
            FRTGroup.create_dataset("mFRT", shape=frtData.mFRT.shape, \
                                    dtype=self.dtype, \
                                    data=frtData.mFRT)
            FRTGroup.create_dataset("varFRT", shape=frtData.varFRT.shape, \
                                    dtype=self.dtype, \
                                    data=frtData.varFRT)

    def __read_body__(self, file):
        dSet_isoSurfaceFilePath = file["isosurface_file"]
        self.isoSurfaceFilePath = dSet_isoSurfaceFilePath[()]
        dSet_configFilePath = file["configuration_file"]
        self.configFilePath = dSet_configFilePath[()]

        for key in file:
            if (key != "isosurface_file") and (key != "configuration_file"):
                IsoSurfaceGroup = file[key]
                frtData = self.createFRTData(key)

                XinitGroup = IsoSurfaceGroup["initial_position"]
                dSet_rho = XinitGroup["rho"]
                with dSet_rho.astype(self.dtype):
                    frtData.rho0 = dSet_rho[:]
                dSet_phi = XinitGroup["phi"]
                with dSet_phi.astype(self.dtype):
                    frtData.phi0 = dSet_phi[:]

                FRTGroup = IsoSurfaceGroup["first_return_time"]
                dSet_mFRT = FRTGroup["mFRT"]
                with dSet_mFRT.astype(self.dtype):
                    frtData.mFRT = dSet_mFRT[:]
                dSet_varFRT = FRTGroup["varFRT"]
                with dSet_varFRT.astype(self.dtype):
                    frtData.varFRT = dSet_varFRT[:]

    def __iter__(self):
        self.dataIterator = iter(self.dataSet)
        return self.dataIterator

    def __next__(self):
        return next(self.dataIterator)

    def __eq__(self, other):
        if not isinstance(other, FRTDataFile):
            return NotImplemented

        are_equal = (self.isoSurfaceFilePath == other.isoSurfaceFilePath)
        are_equal = are_equal and (self.configFilePath == other.configFilePath)
        for (dSet1, dSet2) in zip(self.dataSet, other.dataSet):
            are_equal = are_equal and (dSet1 == dSet2)
        return are_equal

    def createFRTData(self, key):
        frtData = FRTData(key)
        self.dataSet.append(frtData)
        return frtData


class TbarData():
    def __init__(self, key, dtype=np.double):
        self.dtype = dtype
        self.key = key
        self.Tbar = 0.0

    def __eq__(self, other):
        if not isinstance(other, TbarData):
            return NotImplemented

        are_equal = (self.dtype == other.dtype)
        are_equal = are_equal and (self.key == other.key)
        are_equal = are_equal  and (self.Tbar == other.Tbar)
        return are_equal


class TbarDataFile(SPhaseFile):
    CLASSTAG = "TbarDataFile"

    def __init__(self, model_name):
        super(TbarDataFile, self).__init__(model_name)
        self.isoSurfaceFilePath = ""
        self.configFilePath = ""
        self.dataSet = []

    def __write_body__(self, file):
        str_dtype = h5py.string_dtype(encoding='utf-8')
        file.create_dataset("isosurface_file", \
                            dtype=str_dtype, \
                            data=self.isoSurfaceFilePath.encode('utf-8'))
        file.create_dataset("configuration_file", \
                            dtype=str_dtype, \
                            data=self.configFilePath.encode('utf-8'))
        for tbarData in self.dataSet:
            IsoSurfaceGroup = file.create_group(tbarData.key)

            TbarGroup = IsoSurfaceGroup.create_group("mean_period")
            TbarGroup.create_dataset("Tbar", shape=(1, ), \
                                     dtype=self.dtype, \
                                     data=tbarData.Tbar)

    def __read_body__(self, file):
        dSet_isoSurfaceFilePath = file["isosurface_file"]
        self.isoSurfaceFilePath = dSet_isoSurfaceFilePath[()]
        dSet_configFilePath = file["configuration_file"]
        self.configFilePath = dSet_configFilePath[()]

        for key in file:
            if (key != "isosurface_file") and (key != "configuration_file"):
                IsoSurfaceGroup = file[key]
                tbarData = self.createTbarData(key)

                scalar_dummy = np.zeros((1,), dtype=self.dtype)

                tbarGroup = IsoSurfaceGroup["mean_period"]
                dSet_tbar = tbarGroup.get("Tbar")
                with dSet_tbar.astype(self.dtype):
                    dSet_tbar.read_direct(scalar_dummy)
                    tbarData.Tbar = scalar_dummy.copy()


    def __iter__(self):
        self.dataIterator = iter(self.dataSet)
        return self.dataIterator

    def __next__(self):
        return next(self.dataIterator)

    def __eq__(self, other):
        if not isinstance(other, TbarDataFile):
            return NotImplemented

        are_equal = (self.isoSurfaceFilePath == other.isoSurfaceFilePath)
        are_equal = are_equal and (self.configFilePath == other.configFilePath)
        for (dSet1, dSet2) in zip(self.dataSet, other.dataSet):
            are_equal = are_equal and (dSet1 == dSet2)
        return are_equal

    def createTbarData(self, key):
        tbarData = TbarData(key)
        self.dataSet.append(tbarData)
        return tbarData


class IsoSurfaceCorr():
    def __init__(self, key, dtype=np.double):
        self.dtype = np.double
        self.key = key
        self.N = np.array([], dtype=np.uint32)
        self.offset = np.array([], dtype=np.uint32)
        self.sub_pop_size = np.array([], dtype=np.uint32)
        self.cv = np.array([], dtype=self.dtype)
        self.Err_cv = np.array([], dtype=self.dtype)
        self.rho_k = np.array([], dtype=self.dtype)
        self.Err_rho_k = np.array([], dtype=self.dtype)

    def __eq__(self, other):
        if not isinstance(other, IsoSurfaceCorr):
            return NotImplemented
        are_equal = (self.dtype == other.dtype)
        are_equal = are_equal and (self.key == other.key)
        are_equal = are_equal and np.array_equal(self.N, other.N)
        are_equal = are_equal and np.array_equal(self.offset, other.offset)
        are_equal = are_equal and np.array_equal(self.sub_pop_size, other.sub_pop_size)
        are_equal = are_equal and np.array_equal(self.cv, other.cv)
        are_equal = are_equal and np.array_equal(self.Err_cv, other.Err_cv)
        are_equal = are_equal and np.array_equal(self.rho_k, other.rho_k)
        are_equal = are_equal and np.array_equal(self.Err_rho_k, other.Err_rho_k)
        return are_equal


class SerialCorrFile(SPhaseFile):
    CLASSTAG = "SerialCorrelationFile"

    def __init__(self, model_name):
        super(SerialCorrFile, self).__init__(model_name)
        self.isoSurfaceFilePath = ""
        self.configFilePath = ""
        self.dataSet = []

    def __write_body__(self, file):
        str_dtype = h5py.string_dtype(encoding='utf-8')
        file.create_dataset("isosurface_file", \
                            dtype=str_dtype, \
                            data=self.isoSurfaceFilePath.encode('utf-8'))
        file.create_dataset("configuration_file", \
                            dtype=str_dtype, \
                            data=self.configFilePath.encode('utf-8'))
        for frtCorrData in self.dataSet:
            IsoSurfaceGroup = file.create_group(frtCorrData.key)
            IsoSurfaceGroup.create_dataset("total_number_of_intervals", shape=frtCorrData.N.shape, \
                                           dtype=np.uint32, \
                                           data=frtCorrData.N)
            IsoSurfaceGroup.create_dataset("stationary_offset", shape=frtCorrData.offset.shape, \
                                           dtype=np.uint32, \
                                           data=frtCorrData.offset)
            IsoSurfaceGroup.create_dataset("sub_population_size", shape=frtCorrData.sub_pop_size.shape, \
                                           dtype=np.uint32, \
                                           data=frtCorrData.sub_pop_size)
            IsoSurfaceGroup.create_dataset("coefficient_of_variation", shape=frtCorrData.cv.shape, \
                                           dtype=self.dtype, \
                                           data=frtCorrData.cv)
            IsoSurfaceGroup.create_dataset("Err_coefficient_of_variation", shape=frtCorrData.Err_cv.shape, \
                                           dtype=self.dtype, \
                                           data=frtCorrData.Err_cv)
            IsoSurfaceGroup.create_dataset("rho_k", shape=frtCorrData.rho_k.shape, \
                                           dtype=self.dtype, \
                                           data=frtCorrData.rho_k)
            IsoSurfaceGroup.create_dataset("Err_rho_k", shape=frtCorrData.Err_rho_k.shape, \
                                           dtype=self.dtype, \
                                           data=frtCorrData.Err_rho_k)

    def __read_body__(self, file):
        dSet_isoSurfaceFilePath = file["isosurface_file"]
        self.isoSurfaceFilePath = dSet_isoSurfaceFilePath[()]
        dSet_configFilePath = file["configuration_file"]
        self.configFilePath = dSet_configFilePath[()]

        for key in file:
            if (key != "isosurface_file") and (key != "configuration_file"):
                IsoSurfaceGroup = file[key]
                frtCorrData = self.createIsoSurfaceCorr(key)

                dSet_N = IsoSurfaceGroup.get("total_number_of_intervals")
                with dSet_N.astype(np.uint32):
                    frtCorrData.N = np.array(dSet_N, dtype=np.uint32)

                dSet_offset = IsoSurfaceGroup.get("stationary_offset")
                with dSet_offset.astype(np.uint32):
                    frtCorrData.offset = np.array(dSet_offset, dtype=np.uint32)

                dSet_sub_pop_size = IsoSurfaceGroup.get("sub_population_size")
                with dSet_sub_pop_size.astype(np.uint32):
                    frtCorrData.sub_pop_size = np.array(dSet_sub_pop_size, dtype=np.uint32)

                dSet_cv = IsoSurfaceGroup.get("coefficient_of_variation")
                with dSet_cv.astype(self.dtype):
                    frtCorrData.cv = np.array(dSet_cv, dtype=self.dtype)

                dSet_Err_cv = IsoSurfaceGroup.get("Err_coefficient_of_variation")
                with dSet_Err_cv.astype(self.dtype):
                    frtCorrData.Err_cv = np.array(dSet_Err_cv, dtype=self.dtype)

                dSet_rho_k = IsoSurfaceGroup.get("rho_k")
                with dSet_rho_k.astype(self.dtype):
                    frtCorrData.rho_k = np.array(dSet_rho_k, dtype=self.dtype)

                dSet_Err_rho_k = IsoSurfaceGroup.get("Err_rho_k")
                with dSet_Err_rho_k.astype(self.dtype):
                    frtCorrData.Err_rho_k = np.array(dSet_Err_rho_k, dtype=self.dtype)

    def __iter__(self):
        self.dataIterator = iter(self.dataSet)
        return self.dataIterator

    def __next__(self):
        return next(self.dataIterator)

    def __eq__(self, other):
        are_equal = True
        for (d1, d2) in zip(self.dataSet, other.dataSet):
            are_equal = are_equal and (d1 == d2)
        return are_equal

    def createIsoSurfaceCorr(self, key):
        frtCorrData = IsoSurfaceCorr(key)
        self.dataSet.append(frtCorrData)
        return frtCorrData


class SimConfigFile:
    def __init__(self):
        self.modelName = ""
        self.pSets = []
        self.domainName = ""
        self.domainDims = []
        self.paths = {
            "In": Path(""),
            "Out": Path("")
        }
        self.paths_are_not_absolute = True
        self.simulation = {
            "dt": 0.0,
            "t0": 0.0,
            "T": 0.0,
            "Ensemble Size": 0,
            "Sample Size": 0,
        }

    def read(self, file_path: Path):
        with open(file_path, "r") as file:
            dict_bf = json.load(file)

        modelG = dict_bf["Model"]
        self.modelName = modelG["Name"]
        for pSet in modelG["Parameterizations"]:
            if type(pSet) is dict:
                P = {}
                for key in pSet:
                    if type(pSet[key]) is float:
                        P[key] = pSet[key]
                    else:
                        raise SPhaseFileWrongFormat(file_path, "in Model: " \
                                                    + self.modelName + " Parameterizations: Parameter " \
                                                    + key + " has non-float value")
                self.pSets.append(P)
            else:
                raise SPhaseFileWrongFormat(file_path, "in Model: " \
                                            + self.modelName + "Parameterizations: not in dict format")

        domainG = dict_bf["Domain"]
        self.domainName = domainG["Name"]
        for dims in domainG["Dimensions"]:
            if type(dims) is dict:
                D = {}
                for key in dims:
                    if type(dims[key]) is float:
                        D[key] = dims[key]
                    else:
                        raise SPhaseFileWrongFormat(file_path, "in Domain: " \
                                                    + self.domainName + " Dimensions: Dimension " \
                                                    + key + " has non-float value")
                self.domainDims.append(D)
            else:
                raise SPhaseFileWrongFormat(file_path, "in Domain : " \
                                            + self.domainName + " Dimensions: not in dict format")

        filepath_input = Path(dict_bf["Paths"]["In"]["filepath"])
        if filepath_input == Path(""):
            filepath_input = Path(".")
        self.paths["In"] = filepath_input / \
                           Path(dict_bf["Paths"]["In"]["filename"])

        filepath_output = Path(dict_bf["Paths"]["Out"]["filepath"])
        if filepath_output == Path(""):
            filepath_output = Path(".")
        self.paths["Out"] = filepath_output / \
                            Path(dict_bf["Paths"]["Out"]["filename"])

        self.paths_are_not_absolute = \
            not (self.paths["In"].is_absolute() and self.paths["Out"].is_absolute())

        if not self.paths["In"].is_absolute():
            base = Path(file_path).absolute().parent
            self.paths["In"] =  base / self.paths["In"]
            self.paths["In"] = self.paths["In"].resolve()

        if not self.paths["Out"].is_absolute():
            base = Path(file_path).absolute().parent
            self.paths["Out"] = base / self.paths["Out"]
            self.paths["Out"] = self.paths["Out"].resolve()

        self.simulation = dict_bf["Simulation"]

    def print_contents(self):
        print("SimConfig")
        print("\t ModelName: ", self.modelName)
        print("\t Parameterizations:")
        for pSet in self.pSets:
            for key in pSet:
                print("\t \t Parameter: " + key + " = " + pSet[key])

        print("\t Paths: \n \t \t In: ", self.paths["In"],
              "\n \t \t Out: ", self.paths["Out"])
        print("\t Simulation: \n \t \t", self.simulation)
