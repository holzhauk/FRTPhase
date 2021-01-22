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

    def __init__(self, model_name: str, dtype: np.dtype =np.double):
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
                        "file version: " + file.attrs.get("version") +\
                        " is not version " + self.VERSION)
        if file.attrs.get("class") != self.CLASSTAG:
            raise SPhaseFileClassConflict(file_path,\
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
        self.rho = np.array([], dtype=self.dtype)
        self.phi = np.array([], dtype=self.dtype)

    def __eq__(self, other):
        if not isinstance(other, Curve):
            return NotImplemented
        are_equal = True
        are_equal = are_equal and (self.dtype == other.dtype)
        are_equal = are_equal and (self.key == other.key)
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
                Parameter.create_dataset(key, shape=(1,),\
                                         dtype=self.dtype,\
                                         data=curve.parameterSet[key])

            CurveGroup = IsoSurfaceGroup.create_group("curve")
            CurveGroup.create_dataset("rho", shape=curve.rho.shape,\
                                      dtype=self.dtype,\
                                      data=curve.rho)
            CurveGroup.create_dataset("phi", shape=curve.phi.shape,\
                                      dtype=self.dtype,\
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
        self.dtype=dtype
        self.key = key
        self.rho0 = np.array([], dtype=self.dtype)
        self.phi0 = np.array([], dtype=self.dtype)
        self.mPhiT = np.array([], dtype=self.dtype)
        self.varPhiT = np.array([], dtype=self.dtype)
        self.mT = np.array([], dtype=self.dtype)
        self.varT = np.array([], dtype=self.dtype)
        self.mFRT = np.array([], dtype=self.dtype)
        self.varFRT = np.array([], dtype=self.dtype)

class FRTDataFile(SPhaseFile):
    CLASSTAG = "FRTDataFile"
    def __init__(self, model_name):
        super(FRTDataFile, self).__init__(model_name)
        self.isoSurfaceFilePath = ""
        self.configFilePath = ""
        self.dataSet = []

    def __write_body__(self, file):
        file.create_dataset("isosurface_file", shape=(1,),\
                            dtype=str,\
                            data=self.isoSurfaceFilePath)
        file.create_dataset("configuration_file", shape=(1,),\
                            dtype=str,\
                            data=self.configFilePath)
        for frtData in self.dataSet:
            IsoSurfaceGroup = file.create_group(frtData.key)

            XinitGroup = IsoSurfaceGroup.create_group("initial_position")
            XinitGroup.create_dataset("rho", shape=frtData.rho0.shape,\
                                      dtype=self.dtype,\
                                      data=frtData.rho0)
            XinitGroup.create_dataset("phi", shape=frtData.phi0.shape,\
                                      dtype=self.dtype,\
                                      data=frtData.phi0)

            PhiTGroup = IsoSurfaceGroup.create_group("asymptotic_angle")
            PhiTGroup.create_dataset("mPhiT", shape=frtData.mPhiT.shape,\
                                     dtype=self.dtype,\
                                     data=frtData.mPhiT)
            PhiTGroup.create_dataset("varPhiT", shape=frtData.varPhiT.shape,\
                                     dtype=self.dtype,\
                                     data=frtData.varPhiT)

            TbarGroup = IsoSurfaceGroup.create_group("stationary_period")
            TbarGroup.create_dataset("mT", shape=frtData.mT.shape,\
                                     dtype=self.dtype,\
                                     data=frtData.mT)
            TbarGroup.create_dataset("varT", shape=frtData.varT.shape,\
                                     dtype=self.dtype,\
                                     data=frtData.varT)

            FRTGroup = IsoSurfaceGroup.create_group("first_return_time")
            FRTGroup.create_dataset("mFRT", shape=frtData.mFRT.shape,\
                                    dtype=self.dtype,\
                                    data=frtData.mFRT)
            FRTGroup.create_dataset("varFRT", shape=frtData.varFRT.shape,\
                                    dtype=self.dtype,\
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

                PhiTGroup = IsoSurfaceGroup["asymptotic_angle"]
                dSet_mPhiT = PhiTGroup["mPhiT"]
                with dSet_mPhiT.astype(self.dtype):
                    frtData.mPhiT = dSet_mPhiT[:]
                dSet_varPhiT = PhiTGroup["varPhiT"]
                with dSet_varPhiT.astype(self.dtype):
                    frtData.varPhiT = dSet_varPhiT[:]

                TbarGroup = IsoSurfaceGroup["stationary_period"]
                dSet_mT = TbarGroup["mT"]
                with dSet_mT.astype(self.dtype):
                    frtData.mT = dSet_mT[:]
                dSet_varT = TbarGroup["varT"]
                with dSet_varT.astype(self.dtype):
                    frtData.varT = dSet_varT[:]

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

    def createFRTData(self, key):
        frtData = FRTData(key)
        self.dataSet.append(frtData)
        return frtData

class SimConfigFile:
    def __init__(self):
        self.modelName = ""
        self.paths = {
            "In": Path(""),
            "Out": Path("")
        }
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
        self.modelName = dict_bf["Model Name"]
        self.paths["In"] = Path(dict_bf["Paths"]["In"]["filepath"]) / \
                           Path(dict_bf["Paths"]["In"]["filename"])
        self.paths["Out"] = Path(dict_bf["Paths"]["Out"]["filepath"]) / \
                            Path(dict_bf["Paths"]["Out"]["filename"])
        self.simulation = dict_bf["Simulation"]

    def print_contents(self):
        print("SimConfig")
        print("\t ModelName: ", self.modelName)
        print("\t Paths: \n \t \t In: ", self.paths["In"],
              "\n \t \t Out: ", self.paths["Out"])
        print("\t Simulation: \n \t \t", self.simulation)

