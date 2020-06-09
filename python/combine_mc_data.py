# -*- coding: utf-8 -*-

import os
import sys
import json
import numpy as np
from os import path

####
####
####
def addData (dataSets, data, filter):
  ####
  #### unload data
  ####
  L           = data["L"]
  N           = data["N"]
  beta        = data["beta"]
  t           = data["t"]
  J           = data["J"]
  maxNumWorms = data["maxNumWorms"]

  avgSign = data["avgSign"]
  numData = data["numData"]

  ####
  #### filter
  ####
  if "beta"        in filter and filter["beta"]        != beta:        return
  if "t"           in filter and filter["t"]           != t:           return
  if "J"           in filter and filter["J"]           != J:           return
  if "L"           in filter and filter["L"]           != L:           return
  if "maxNumWorms" in filter and filter["maxNumWorms"] != maxNumWorms: return
  if "betas"       in filter and beta not in filter["betas"]:          return


  ####
  #### must have some data
  ####
  if numData < 1e3:
    print("numData =", numData, "is insufficient")
    return

  # halfNumNNs = 2;
  key = "L={}, N={}, tBeta={}, JBeta={},".format(L, N, t * beta, J * beta)
  if key in dataSets:
    dataSets[key]["count"]   += 1
    dataSets[key]["numData"] += [ numData ]
    dataSets[key]["avgSign"] += [ avgSign ]

    if "C1" in data: dataSets[key]["C1"] += [ data["C1"] ]
    if "C2" in data: dataSets[key]["C2"] += [ data["C2"] ]
    if "C3" in data: dataSets[key]["C3"] += [ data["C3"] ]

    if "avgKinetEnergy"   in data: dataSets[key]["avgKinetEnergy"]   = dataSets[key].get("avgKinetEnergy",   []) + [ data["avgKinetEnergy"] ]
    if "avgExchaEnergy"   in data: dataSets[key]["avgExchaEnergy"]   = dataSets[key].get("avgExchaEnergy",   []) + [ data["avgExchaEnergy"] ]
    if "avgPotenEnergy"   in data: dataSets[key]["avgPotenEnergy"]   = dataSets[key].get("avgPotenEnergy",   []) + [ data["avgPotenEnergy"] ]
    if "avgInterEnergy"   in data: dataSets[key]["avgInterEnergy"]   = dataSets[key].get("avgInterEnergy",   []) + [ data["avgInterEnergy"] ]
    if "avgNnInterEnergy" in data: dataSets[key]["avgNnInterEnergy"] = dataSets[key].get("avgNnInterEnergy", []) + [ data["avgNnInterEnergy"] ]

  else:
    dataSets[key] = {
      "L":        L,
      "numHoles": N,
      "beta":     beta,
      "t":        t,
      "J":        J,
      "count":    1,
      "numData":  [ numData ],
      "avgSign":  [ avgSign ]
    }

    if "C1" in data: dataSets[key]["C1"] = [ data["C1"] ]
    if "C2" in data: dataSets[key]["C2"] = [ data["C2"] ]
    if "C3" in data: dataSets[key]["C3"] = [ data["C3"] ]

    if "avgKinetEnergy"   in data: dataSets[key]["avgKinetEnergy"]   = [ data["avgKinetEnergy"] ]
    if "avgExchaEnergy"   in data: dataSets[key]["avgExchaEnergy"]   = [ data["avgExchaEnergy"] ]
    if "avgPotenEnergy"   in data: dataSets[key]["avgPotenEnergy"]   = [ data["avgPotenEnergy"] ]
    if "avgInterEnergy"   in data: dataSets[key]["avgInterEnergy"]   = [ data["avgInterEnergy"] ]
    if "avgNnInterEnergy" in data: dataSets[key]["avgNnInterEnergy"] = [ data["avgNnInterEnergy"] ]


####
####
####
def combine_mc_data (jobNames, filter):

  if type(jobNames) is not list:
    print("\"jobNames\" not a list")
    return

  dataSets = {}
  for jobName in jobNames:
    jobRoot = os.getcwd() + "/data/" + jobName

    jobs = []
    print("processing:", jobRoot)
    for root, dirs, files in os.walk(jobRoot):

      # no log directory
      if "logs" in root: continue
      # look for statistics.json
      if "statistics.json" in files:jobs += [root]


    for jobFolder in jobs:

      # add "/"
      jobFolder += "/"

      ####
      #### fetch parameters from json file
      ####
      with open(jobFolder + "../parameters.json", "r") as parameterFile:
        parameters = json.load(parameterFile)

      # add parameters from the other parameter file
      with open(jobFolder + "used-parameters.json", "r") as usedParameterFile:
        usedParameters = json.load(usedParameterFile)


      ####
      #### fetch statistics
      ####
      jsonPath = jobFolder + "statistics.json"
      if os.path.isfile(jsonPath): statistics = json.load(open(jsonPath))
      else:                        continue


      ####
      #### version
      ####
      version = statistics["version"] if "version" in statistics else 0


      ####
      #### parameters
      ####
      data = {}
      data["L"]           = int(parameters["lattice"]["size"][0])
      data["N"]           = int(parameters["model"]["Ns"][1])
      data["beta"]        = parameters["model"]["beta"]
      data["t"]           = parameters["model"]["ts"][0]
      data["J"]           = parameters["model"]["Js"][0]
      data["maxNumWorms"] = int(parameters["model"]["maxNumWorms"]) if "maxNumWorms" in parameters["model"] else 2

      data["avgSign"]   = statistics["avgSign"]
      data["numData"]   = statistics["numData"]

      halfNumNNs = 2;
      shape = (data["L"], data["L"], halfNumNNs)
      if version >= 5:
        if "C1" in statistics: data["C1"] = np.resize(statistics["C1"], shape)
        if "C2" in statistics: data["C2"] = np.resize(statistics["C2"], shape)
        if "C3" in statistics: data["C3"] = np.resize(statistics["C3"], shape)

      if version >= 8:
        if "avgKinetEnergy"   in statistics: data["avgKinetEnergy"]   = statistics["avgKinetEnergy"]
        if "avgExchaEnergy"   in statistics: data["avgExchaEnergy"]   = statistics["avgExchaEnergy"]
        if "avgPotenEnergy"   in statistics: data["avgPotenEnergy"]   = statistics["avgPotenEnergy"]
        if "avgInterEnergy"   in statistics: data["avgInterEnergy"]   = statistics["avgInterEnergy"]
        if "avgNnInterEnergy" in statistics: data["avgNnInterEnergy"] = statistics["avgNnInterEnergy"]


      addData(dataSets, data, filter)

  return dataSets
