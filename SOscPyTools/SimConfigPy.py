#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:18:57 2020

@author: konstantin
"""

import sys, os

from pathlib import Path
import json

class SimConfig:
   
    def __init__(self):
        self.ModelName = ""
        self.Paths = { 
                        "In": "",
                        "Out": ""
                        }
        self.Simulation = {
                            "dt": 0.0, 
                            "t0": 0.0,
                            "T": 0.0,
                            "Ensemble Size": 0,
                            "Sample Size": 0,
                            }
    
    def load_from_file(self, configFilePath):
        with open(Path(configFilePath), "r") as read_file:
            dict_bf = json.load(read_file)
        self.ModelName = dict_bf["Model Name"]
        self.Paths["In"] = configFilePath.parent / \
            Path(dict_bf["Paths"]["In"]["filepath"]) / \
            Path(dict_bf["Paths"]["In"]["filename"])
        self.Paths["Out"] = configFilePath.parent / \
            Path(dict_bf["Paths"]["Out"]["filepath"]) / \
            Path(dict_bf["Paths"]["Out"]["filename"])
        self.Simulation = dict_bf["Simulation"]
        
    def print_contents(self):
        print("SimConfig")
        print("\t ModelName: ", self.ModelName)
        print("\t Paths: \n \t \t In: ", self.Paths["In"],
              "\n \t \t Out: ", self.Paths["Out"])
        print("\t Simulation: \n \t \t", self.Simulation)


if __name__ == "__main__":
    
    if len(sys.argv) != 2:
        print("Usage: ", sys.argv[0], " <SimConfigFile>.json")
        os._exit(os.EX_IOERR)
        
    simConfig = SimConfig()
    simConfig.load_from_file(sys.argv[1])
    simConfig.print_contents()
    
    