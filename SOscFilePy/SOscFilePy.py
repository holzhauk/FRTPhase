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
    

class ScalarSet(SOscFile):

    CLASSTAG = "ScalarSet"
    
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
        
        model_parameters = {}
        simulation_parameters = {}
        
        E = self.fh[next(self.fh_iter)]
        
        dummy = np.zeros((1, ), dtype=self.dtype)
        ModelParameter = E["Parameter/Model"]
        for key in ModelParameter:
            ModelParameter[key].read_direct(dummy)
            model_parameters[key] = np.copy(dummy)
        
        SimulationParameter = E["Parameter/Simulation"]
        for key in SimulationParameter:
            SimulationParameter[key].read_direct(dummy)
            simulation_parameters[key] = np.copy(dummy)
            
        dset = E["Data"]
        with dset.astype(self.dtype):
            data = dset[:]
        
        return (model_parameters, \
                    simulation_parameters, \
                    data)
    
    
    def add_Ensemble(self, model_parameters, \
                     simulation_parameters, \
                     data):
        NofEnsembles = len(self.fh.keys())
        E = self.fh.create_group("Ensemble" + str(NofEnsembles))
        P = E.create_group("Parameter")
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
            model_parameters = {}
            simulation_parameters = {}
            
            E = self.fh["Ensemble" + str(i)]
            dummy = np.zeros((1, ), dtype=self.dtype)
            ModelParameter = E["Parameter/Model"]
            for key in ModelParameter:
                ModelParameter[key].read_direct(dummy)
                model_parameters[key] = np.copy(dummy)
            
            SimulationParameter = E["Parameter/Simulation"]
            for key in SimulationParameter:
                SimulationParameter[key].read_direct(dummy)
                simulation_parameters[key] = np.copy(dummy)
                
            dset = E["Data"]
            with dset.astype(self.dtype):
                data = dset[:]
                
            return (model_parameters, \
                    simulation_parameters, \
                    data)
                
            
if __name__ == "__main__":
    
    ################################
    # Test of the SOscFilePy module
    ################################
    
    print("RUNNING TEST OF THE SOscFilePy module")
    print("author: Konstantin Holzhausen")
    
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
    
    print("Test endet successfully")
