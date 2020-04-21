#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:07:22 2020

@author: konstantin
"""

import h5py
import os
import numpy as np


class SOscFileError(Exception):
    pass

class SOscWrongFileFormat(SOscFileError):
    
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
        
class SOscWrongFile(SOscFileError):
    
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message



class SOscFile():
    
    FORMAT = "SOscFile"
    VERSION = "0.0.1"
    CLASSTAG = ""

    def __init__(self, filepath, modelname, dtype=np.double):
        self.FILEPATH = filepath
        self.MODELNAME = modelname
        self.dtype = dtype
        
        if os.path.isfile(self.FILEPATH):
            self.fh = h5py.File(self.FILEPATH, "r+")
            
            if (self.fh.attrs.get("format") != self.FORMAT) or \
                (self.fh.attrs.get("version") != self.VERSION):
                raise SOscWrongFileFormat
            
            if (self.fh.attrs.get("model") != self.MODELNAME) or \
                (self.fh.attrs.get("class") != self.CLASSTAG):
                raise SOscWrongFile
            
        else:
            self.fh = h5py.File(self.FILEPATH, "w")
            self.fh.attrs.create("format", self.FORMAT)
            self.fh.attrs.create("version", self.VERSION)
            self.fh.attrs.create("class", self.CLASSTAG)
            self.fh.attrs.create("model", self.MODELNAME)
        
    
    def __del__(self):
        self.fh.close()
    
    
    def __iter__(self):
        self.fh_iter = iter(self.fh)
        return self

    
    def __next__(self):
        pass


class IsochroneSet(SOscFile):

    CLASSTAG = "IsochroneSet"

    def __get_Isochrone__(self, group_key):

        parameters = {}

        I = self.fh[group_key]

        dummy = np.zeros((1,), dtype=self.dtype)
        P = I["Parameters"]
        for key in P:
            P[key].read_direct(dummy)
            parameters[key] = np.copy(dummy)

        C = I["Curve"]
        dset_rho = C["Rho"]
        with dset_rho.astype(self.dtype):
            Rho = dset_rho[:]
        dset_phi = C["Phi"]
        with dset_phi.astype(self.dtype):
            Phi = dset_phi[:]

        return (parameters, \
                    Rho,\
                    Phi)


    def __next__(self):
        return self.__get_Isochrone__(next(self.fh_iter))

    
    def add_Isochrone(self, parameters, Rho, Phi):
        NofIsochrones = len(self.fh.keys()) # number of isochrones in the file
        I = self.fh.create_group("Isochrone" \
                            + str(NofIsochrones))
        P = I.create_group("Parameters")
        for key in parameters:
            P.create_dataset(key, shape=(1,),\
                                dtype=self.dtype,\
                                data=parameters[key])
        C = I.create_group("Curve")
        C.create_dataset("Rho", shape=Rho.shape,\
                                dtype=self.dtype,\
                                data=Rho)
        C.create_dataset("Phi", shape=Phi.shape,\
                                dtype=self.dtype,\
                                data=Phi)


    def get_Isochrone(self, i):
        if not self.fh.__contains__("Isochrone" + str(i)):
            return ({}, \
                    np.array([], dtype=self.dtype), \
                        np.array([], dtype=self.dtype))
        else:
            return self.__get_Isochrone__("Isochrone" + str(i))
        

class TrajectorySet(SOscFile):
    
    CLASSTAG = "TrajectorySet"
    
    def __get_Parameters__(self):
        
        parameters = {}
        
        dummy = np.zeros((1,), dtype=self.dtype)
        P_group = self.fh["Parameters"]
        for key in P_group:
            P_group[key].read_direct(dummy)
            parameters[key] = np.copy(dummy)
            
        return parameters
    
    
    def __get_Trajectory__(self, group_key):
        
        T_group = self.fh["Trajectories"]
        Tset = T_group[group_key]
        
        dset_rho = Tset["Rho"]
        with dset_rho.astype(self.dtype):
            Rho = dset_rho[:]
        
        dset_phi = Tset["Phi"]
        with dset_phi.astype(self.dtype):
            Phi = dset_phi[:]
            
        dset_t = Tset["t"]
        with dset_phi.astype(self.dtype):
            t = dset_t[:]
            
        return (Rho, \
                Phi, \
                    t)
            
    
    def __iter__(self):
        self.Tset_iter = iter(self.fh["Trajectories"])
        return self
    
    
    def __next__(self):
        return self.__get_Trajectory__(next(self.Tset_iter))
    
    
    def get_Trajectory(self, i):
        if not self.fh["Trajectories"].__contains__("T" + str(i)):
            return (np.array([], dtype=self.dtype), \
                    np.array([], dtype=self.dtype), 
                    np.array([], dtype=self.dtype))
        else:
            return self.__get_Trajectory__("T" + str(i))


class Isochrone:
    
    def __init__(self, dtype=np.double):
        self.dtype = dtype
        self.key = ""
        self.EnsembleSize = 0
        self.Rho_init = np.array([], dtype=self.dtype)
        self.Phi_init = np.array([], dtype=self.dtype)
        self.MFPT = np.array([], dtype=self.dtype)
        self.VarFPT = np.array([], dtype=self.dtype)


class MFPTSet(SOscFile):
    
    CLASSTAG = "MeanFirstPassageTimesSet"
    
    def get_IsochroneFilePath(self):
        return self.fh.attrs.get("IsochroneFile")
    
    
    def __get_Isochrone__(self, group_key):
        
        I = Isochrone();
        I.key = group_key
        
        I_group = self.fh[group_key]

        dummy = np.zeros((1,), dtype=np.int)
        I_group["EnsembleSize"].read_direct(dummy)
        I.EnsembleSize = np.copy(dummy)

        Init_conds = I_group["InitialPositions"]
        
        dset_rho = Init_conds["Rho"]
        with dset_rho.astype(self.dtype):
            I.Rho_init = dset_rho[:]
            
        dset_phi = Init_conds["Phi"]
        with dset_phi.astype(self.dtype):
            I.Phi_init = dset_phi[:]
            
        dset_mfpts = I_group["MFPT"]
        with dset_mfpts.astype(self.dtype):
            I.MFPT = dset_mfpts[:]
            
        dset_VarFPT = I_group["VarFPT"]
        with dset_VarFPT.astype(self.dtype):
            I.VarFPT = dset_VarFPT[:]
            
        return I
    
    
    def __next__(self):
        return self.__get_Isochrone__(next(self.fh_iter))

    

class ScalarSet(SOscFile):

    CLASSTAG = "ScalarSet"

    def __get_Ensemble__(self, group_key):

        model_parameters = {}
        simulation_parameters = {}
        
        E = self.fh[group_key]
        
        dummy = np.zeros((1, ), dtype=self.dtype)
        ModelParameter = E["Parameters/Model"]
        for key in ModelParameter:
            ModelParameter[key].read_direct(dummy)
            model_parameters[key] = np.copy(dummy)
        
        SimulationParameter = E["Parameters/Simulation"]
        for key in SimulationParameter:
            SimulationParameter[key].read_direct(dummy)
            simulation_parameters[key] = np.copy(dummy)
            
        dset = E["Data"]
        with dset.astype(self.dtype):
            data = dset[:]

        return (model_parameters, \
                    simulation_parameters, \
                    data)


    def __next__(self):
        return self.__get_Ensemble__(next(self.fh_iter))
        
    
    def add_Ensemble(self, model_parameters, \
                     simulation_parameters, \
                     data):
        NofEnsembles = len(self.fh.keys())
        E = self.fh.create_group("Ensemble" + str(NofEnsembles))
        P = E.create_group("Parameters")
        ModelParameters = P.create_group("Model")
        for key in model_parameters:
            ModelParameters.create_dataset(key, shape=(1,), \
                                           dtype=self.dtype, \
                                           data=model_parameters[key])
        SimulationParameters = P.create_group("Simulation")
        for key in simulation_parameters:
            SimulationParameters.create_dataset(key, shape=(1,), \
                                                dtype=self.dtype, \
                                                data=simulation_parameters[key])
        E.create_dataset("Data", shape=data.shape, \
                         dtype=self.dtype,\
                         data=data)
            
    
    def get_Ensemble(self, i):
        if not self.fh.__contains__("Ensemble" + str(i)):
            return ({}, {}, np.array([], dtype=self.dtype))
        else:
            return self.__get_Ensemble__("Ensemble" + str(i))
                            
            
if __name__ == "__main__":
    
    ################################
    # Test of the SOscFilePy module
    ################################
    
    print("RUNNING TEST OF THE SOscFilePy module")
    print("author: Konstantin Holzhausen")

    def test_ScalarSet():
        FILENAME = "test.h5"
        
        ModelP = {
            "D": 0.198,
            "omega": 1.0,
            "gamma": 15.0,
            "c": -15.0
            }
        
        SimulationP = {
            "t0": 0.0,
            "dt": 0.001,
            "T": 100.0,
            "Rho0": 1.6,
            "Phi0": np.pi / 2
            }
        
        dtype = np.double
        
        if os.path.isfile(FILENAME):
            print("FILE " + FILENAME + " already exists -> delete FILE")
            os.remove(FILENAME)
        
        print("Create FILE " + FILENAME)
        F = ScalarSet(FILENAME, "testmodel")
        
        print("Add data to FILE")
        Tbar_test = np.random.randn(100).astype(dtype)
        F.add_Ensemble(ModelP, SimulationP, Tbar_test)
        
        print("Read ONE dataset from FILE")
        print(F.get_Ensemble(0))
        
        print("Read ALL datasets from FILE")
        for Ensemble in F:
            print(Ensemble)

        print("Remove Class object ...")
        del F
            
        print("Delete FILE ...")
        os.remove(FILENAME)

    def test_IsochroneSet():
        FILENAME = "test.h5"
        
        P = {
                "D": 0.198,
                "omega": 1.0,
                "gamma": 15.0,
                "c": -15.0
                }

        dtype = np.double

        if os.path.isfile(FILENAME):
            print("FILE " + FILENAME + " already exists -> delete FILE")
            os.remove(FILENAME)
        
        print("Create FILE " + FILENAME)
        F = IsochroneSet(FILENAME, "testmodel")
        
        print("Add data to FILE")
        Phi = np.linspace(0.0, 2*np.pi, num=100)
        Rho = np.ones(Phi.shape)*1.0
        F.add_Isochrone(P, Rho, Phi)
        
        print("Read ONE dataset from FILE")
        print(F.get_Isochrone(0))
        
        print("Read ALL datasets from FILE")
        for Isochrone in F:
            print(Isochrone)

        print("Remove Class object ...")
        del F
            
        print("Delete FILE ...")
        os.remove(FILENAME)
    
    
    test_ScalarSet()
    test_IsochroneSet()
    print("Test endet successfully")
