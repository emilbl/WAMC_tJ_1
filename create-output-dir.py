import os
import sys
import json
from time import localtime, strftime   # to get the local time in a readable format
import importlib.util   # dynamically import modules
import shutil
import re
import traceback   # print traceback of error


####
#### function determining the folder name given input parameters
####
def namePattern (parameters):

  # convert to list or number string
  def toString (strOrList, delimiter = "x"):
    if isinstance(strOrList, list):
      return delimiter.join(map(str, strOrList))
    else:
      return str(strOrList)

  # canonical components and their particle number
  if isinstance(parameters["model"]["is-canonical"], list):
    temp = []
    for i in range(0, len(parameters["model"]["is-canonical"])):
      temp += [str(parameters["model"]["Ns"][i]) if parameters["model"]["is-canonical"][i] else "G"]
    canonical = toString(temp)
  else:
    canonical = str(parameters["model"]["Ns"]) if parameters["model"]["is-canonical"] else "G"

  # the name formula
  return parameters["lattice"]["name"] + "-" + toString(parameters["lattice"]["size"]) + "_" \
       + parameters["model"]["name"] + "-" + canonical + "_" \
       + "beta=" + str(parameters["model"]["beta"]) + "_" \
       + "ts=" + toString(parameters["model"]["ts"]) + "_" \
       + ("Us=" + toString(parameters["model"]["Us"]) + "_" if "Us" in parameters["model"] else "") \
       + ("Js=" + toString(parameters["model"]["Js"]) + "_" if "Js" in parameters["model"] else "") \
       + "mus=" + toString(parameters["model"]["mus"]) + "_" \
       + strftime("%y-%m-%d_%H:%M:%S", localtime())

####
#### job folder containing all parameters sets
####
rootFolder = os.getcwd() + "/" # + "/data/job-parameter-sets/"

####
#### path to the parameter file relative the jobs folder
####
if len(sys.argv) < 2:
  print("create-output-dir.py: no parameter file supplied  ->  EXIT")
  # exit with "error" (important so that the shell script will exit itself)
  sys.exit("error")
parameterFilePath = sys.argv[1]

####
#### check if parameter file exist
####
if not os.path.isfile(rootFolder + parameterFilePath):
  print("create-output-dir.py: the parameter file does not exist  ->  EXIT")
  print(rootFolder)
  print(parameterFilePath)
  print(rootFolder + parameterFilePath)
  # exit with "error" (important so that the shell script will exit itself)
  sys.exit("error")


####
#### Ensure that the folder and required contents of previous job exist
####
if len(sys.argv) == 3:
  prevJobFolder = sys.argv[2]
  prevJobUsedParameters = prevJobFolder + ("" if prevJobFolder.endswith("/") else "/") + "used-parameters.json"
  prevJobSavedConfiguration = prevJobFolder + ("" if prevJobFolder.endswith("/") else "/") + "saved-configuration.json"

  # check if configuration file exist
  if not os.path.isdir(rootFolder + prevJobFolder):
    print("create-output-dir.py: the previous job folder does not exist  ->  EXIT")
    # exit with "error" (important so that the shell script will exit itself)
    sys.exit("error")

    # check content
    if not os.path.isfile(rootFolder + prevJobUsedParameters):
      print("create-output-dir.py: the \"used-parameters.json\" of the previous job does not exist  ->  EXIT")
      # exit with "error" (important so that the shell script will exit itself)
      sys.exit("error")

    # check content
    if not os.path.isfile(rootFolder + prevJobSavedConfiguration):
      print("create-output-dir.py: the \"saved-configuration.json\" of the previous job does not exist  ->  EXIT")
      # exit with "error" (important so that the shell script will exit itself)
      sys.exit("error")

####
#### store the parameter set(s) as a json object
####
if parameterFilePath.endswith('.py'):
  # is a python file which will generate one or more parameter sets

  ####
  #### import python module
  ####
  spec = importlib.util.spec_from_file_location("generate", rootFolder + parameterFilePath)
  module = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(module)

  # generate parameters sets
  parameterSets = module.generate()

else:
  # is a static json file containing the parameter set

  with open(rootFolder + parameterFilePath, "r") as jsonFile:
    # obtain the content
    content = jsonFile.read()

    # remove comments since the json parse cant handle those
    content = re.sub(r'//.*\n', '\n', content)

    # parse string as json
    parameterSets = [json.loads(content)]

####
#### data output folder
####
jobsFolder = os.getcwd() + "/data/jobs/"

####
#### clear previous dev directory if development mode is activated
####
if parameterSets[0]["development"]:
  devMode = True

  if os.path.exists(jobsFolder + "dev"):
    # remove sudo
    os.umask(0)
    # rm
    shutil.rmtree(jobsFolder + "dev", ignore_errors=True)
else:
  devMode = False

####
#### create wrapping folder if it does not already exist
####
date = strftime("%y-%m-%d", localtime())
wrappingFolder = jobsFolder + ("dev" if devMode else date)
if not os.path.exists(wrappingFolder):
  # remove sudo
  os.umask(0)
  # mkdir
  os.makedirs(wrappingFolder)


####
#### map parameter sets to file names
####
_names = {}
for i, parameters in enumerate(parameterSets):
  name = namePattern(parameters)

  if name in _names: _names[name] += [(i, json.dumps(parameters))]
  else:              _names[name]  = [(i, json.dumps(parameters))]


####
#### fix subdirectories to ensure unique paths for each parameter set
####
names = {}
for name, _parametersSets in _names.items():
  ####
  #### prepare a dictionary of unique sets
  ####
  sets = {}
  for _parametersSet in _parametersSets:
    parameters = _parametersSet[1]
    if not parameters in sets: sets[parameters] = None

  numUniqueJobs = len(sets)

  # if numUniqueJobs == 1:
  #   # print(_parametersSets)
  #   names[name] = {}
  #   names[name]["jobIds"] = [_parametersSets[0][0]]
  #   names[name]["parameters"] = _parametersSets[0][1]
  # else:

  count = 0
  for _parametersSet in _parametersSets:
    parametersSetIndex = _parametersSet[0]
    parameters = _parametersSet[1]



    if sets[parameters] == None:
      subdir = count
      key = name + ("/" + str(subdir) if numUniqueJobs > 1 else "")

      names[key] = {
        "jobIds":     [parametersSetIndex],
        "parameters": parameters
      }

      # set subdir
      sets[parameters] = subdir

      # increment counter
      count += 1
    else:
      subdir = sets[parameters]
      key = name + ("/" + str(subdir) if numUniqueJobs > 1 else "")

      names[key]["jobIds"] += [parametersSetIndex]

####
#### create directories etc.
####
print("---------[spawning jobs]---------")
jobDirPaths = []
for jobName, val in names.items():

  jobIds = val["jobIds"]
  parameters = val["parameters"]

  print(jobName, jobIds)

  # create directory
  jobDirPath = wrappingFolder + "/" + jobName
  if not os.path.exists(jobDirPath):
    # remove sudo
    os.umask(0)
    # mkdir
    os.makedirs(jobDirPath)

  # dump parameter file in this directory
  newParameterFilePath = jobDirPath + "/parameters.json"
  with open(newParameterFilePath, 'w') as outfile:

    ####
    #### first we should replace the chemical potential with the optimized one, if provided
    ####
    if len(sys.argv) == 3:
      with open(rootFolder + prevJobUsedParameters) as json_file:
        usedParams = json.load(json_file)

        # print(usedParams["model"]["mu_ia"])
        # print(parameters["model"]["mus"])

        # load
        parameters = json.loads(parameters)

        # overwrite
        parameters["model"]["mus"] = usedParams["model"]["mu_ia"]

        # make string
        parameters = json.dumps(parameters)


    json.dump(json.loads(parameters), outfile, sort_keys = True, indent = 2)

  # populate output list
  for id in jobIds:
    jobDirPaths += [jobDirPath + ";" + str(id)]


  ####
  #### copy also input configuration if provided
  ####
  if len(sys.argv) == 3:

    # copy file
    shutil.copyfile(rootFolder + prevJobSavedConfiguration, jobDirPath + "/initial-configuration.json")

print("---------------------------------")

sys.exit("\n".join(jobDirPaths))