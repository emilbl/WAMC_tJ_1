#include <iostream>

#pragma GCC diagnostic push                   // disables warning temporarily
#pragma GCC diagnostic ignored "-Wpedantic"   //
#include <nlohmann/json.hpp>                  // https://github.com/nlohmann/json
#pragma GCC diagnostic pop                    // enables warning once again

#include "etc/cout.h"
#include "settings/settings.h"
#include "Write-to-file/Write-to-file.h"
#include "Pseudorandom/Pseudorandom.h"
#include "Lattice/Lattice.h"
#include "Model/Model.h"
#include "Hamiltonian/Hamiltonian.h"
#include "Worm/Worm.h"
#include "Date-and-time/Date-and-time.h"
#include "Arithmetic-parser/Arithmetic-parser.h"


using namespace std;
using json = nlohmann::json;

using namespace settings::model;

string jobDirPath;
string parametersFilePath;

////
//// $ bin/worm <absolute path to job directory> <unique job id>
////
int main (
  int argc,
  char ** argv
) {

  ////
  //// make sure a job directory path is supplied
  ////
  if (argc < 2) {
    cout << "main: no parameter file supplied  ->  EXIT" << endl;
    exit(EXIT_FAILURE);
  }
  jobDirPath = argv[1];

  // cout parameter file path
  cout << "job directory path:" << endl
       << jobDirPath << endl
       << "---------------------------" << endl;

  ////
  //// make sure the parameter file exist
  ////
  parametersFilePath = jobDirPath + "/parameters.json";
  ifstream file{parametersFilePath};
  if ( ! file.good()) {
    cout << "main: parameter file \"" << parametersFilePath << "\" does not exist  ->  EXIT" << endl;
    exit(EXIT_FAILURE);
  }

  ////
  //// load parameter file
  ////
  json _parameters;
  file >> _parameters;

  ////
  //// add placeholders before making constant
  ////
  if (_parameters["model"].find("Js") == _parameters["model"].end()) {
    _parameters["model"]["Js"] = vector<double>{};
  }
  if (_parameters["model"].find("Us") == _parameters["model"].end()) {
    _parameters["model"]["Us"] = vector<double>{};
  }

  // make constant
  const auto parameters =_parameters;

  // cout input parameter values
  cout << "input parameters:" << endl
       << parameters.dump(2) << endl
       << "---------------------------" << endl;

  ////
  //// make sure a job id is supplied
  ////
  if (argc < 3) {
    cout << "main: no job id supplied  ->  EXIT" << endl;
    exit(EXIT_FAILURE);
  }

  // add to the job id also the rank
  const string jobId = argv[2];



  // cout job id
  cout << "job id: " << jobId << endl
       << "---------------------------" << endl;


  ////
  //// try load configuration file
  ////
  json _configuration;
  if (argc == 4) {
    // a configuration is also sent
    const string configFilePath = jobDirPath + "/" + argv[3];

    // output file path
    cout << "configuration file:" << endl
         << configFilePath << endl
         << "---------------------------" << endl;

    // make sure the file exists
    ifstream file{configFilePath};
    if ( ! file.good()) {
      cout << "main: the configuration file " << configFilePath << " does not exist  ->  EXIT" << endl;
      exit(EXIT_FAILURE);
    }

    file >> _configuration;
  }
  const auto configuration = _configuration;

  ////
  //// announce settings
  ////
  if (debug_major) {
    cout << settings::cout::enterRed
         << "###########################################" << endl
         << "#                                         #" << endl
         << "#   settings::mode::debug_major == true   #" << endl
         << "#                                         #" << endl
         << "###########################################" << endl
         << settings::cout::resetStyle;
  }
  if (debug) {
    cout << settings::cout::enterRed
         << "#####################################" << endl
         << "#                                   #" << endl
         << "#   settings::mode::debug == true   #" << endl
         << "#                                   #" << endl
         << "#####################################" << endl
         << settings::cout::resetStyle;
  }
  if (shutItDown) {
    cout << settings::cout::enterRed
         << "##########################################" << endl
         << "#                                        #" << endl
         << "#   settings::mode::shutItDown == true   #" << endl
         << "#                                        #" << endl
         << "##########################################" << endl
         << settings::cout::resetStyle;
  }
  if (tMax < 1e10) {
    cout << settings::cout::enterRed
         << "###################################" << endl
         << "#                                 #" << endl
         << "#   settings::worm::tMax < 1e10   #" << endl
         << "#                                 #" << endl
         << "###################################" << endl
         << settings::cout::resetStyle;
  }
  if (parameters["development"]) {
    cout << settings::cout::enterRed
         << "#############################" << endl
         << "#                           #" << endl
         << "#   \"development\" == true   #" << endl
         << "#                           #" << endl
         << "#############################" << endl
         << settings::cout::resetStyle;
  }
  if ( ! parameters["seed"]["unique"]) {
    cout << settings::cout::enterRed
         << "################################" << endl
         << "#                              #" << endl
         << "#   \"seed.unique\" == false     #" << endl
         << "#                              #" << endl
         << "################################" << endl
         << settings::cout::resetStyle;
  }


  ////
  //// Parse run time
  ////
  ArithmeticParser arithmeticParser;
  const unsigned numSecWarmUp  = arithmeticParser.parse<unsigned>(parameters["runOptions"]["numSecWarmUp"].get<string>());    // number of seconds to equilibrate the system using simulated annealing
  const unsigned numSecOptChem = arithmeticParser.parse<unsigned>(parameters["runOptions"]["numSecOptChem"].get<string>());   // number of seconds used to optimize the chemical potentials
  const unsigned numSecPreSim  = arithmeticParser.parse<unsigned>(parameters["runOptions"]["numSecPreSim"].get<string>());    // number of seconds used to estimate what should be the bin frequency
  const unsigned numSecSim     = arithmeticParser.parse<unsigned>(parameters["runOptions"]["numSecSim"].get<string>());       // number of seconds used to gather the data
  const unsigned long numData  = parameters["runOptions"]["numData"].get<unsigned long>();                                   // how many data points should be gathered


  ////
  //// start timer
  ////
  Tic startTime{};

  ////
  //// initial time and date
  ////
  const auto dateAndTime = DateAndTime::getBoth();

  ////
  //// check that the number of components is correct
  ////
  if (parameters["model"]["components"] != numComps) {
    cout << settings::cout::enterRed
         << "parameters[\"model\"][\"components\"] != numComps: "
         << parameters["model"]["components"] << " != " << numComps << "  ->  EXIT"
         << settings::cout::resetStyle << endl;
    exit (EXIT_SUCCESS);
  }

  ////
  //// check that the model is correct
  ////
  if (parameters["model"]["name"] != modelName) {
    cout << settings::cout::enterRed
         << "parameters[\"model\"][\"name\"] != modelName: "
         << parameters["model"]["name"] << " != " << modelName << "  ->  EXIT"
         << settings::cout::resetStyle << endl;
    exit (EXIT_SUCCESS);
  }

  ////
  //// initiate file saver which creates output dir
  //// (also append thread id in case of multi core jobs)
  ////
  WriteToFile writeToFile(jobDirPath + "/" + jobId);

  ////
  //// initiate json object in which parameters and coefficients will be stored
  ////
  json usedParameters;

  // add some static parameters from settings
  usedParameters["tMax"]               = tMax;
  usedParameters["allowMultiCompWorm"] = allowMultiCompWorm;
  usedParameters["allowAntiWorm"]      = allowAntiWorm;
  usedParameters["expDistrEnabled"]    = expDistrEnabled;
  usedParameters["maxNumWorms"]        = parameters["model"]["maxNumWorms"];

  ////
  //// initiate and (reseed pseudorandom number generator)
  ////
  Pseudorandom pseudoRandom{parameters["seed"]["value"],
                            parameters["seed"]["unique"]};
  usedParameters["seed"]["value"] = pseudoRandom.getSeed();

  ////
  //// initiate lattice
  ////
  Lattice lattice{parameters["lattice"]["name"],
                  parameters["lattice"]["size"],
                  parameters["lattice"]["parameters"],
                  parameters["lattice"]["isPeriodic"]};
  usedParameters["lattice"]["sites"] = lattice.getNumSites();
  usedParameters["lattice"]["dimensions"] = lattice.getNumDimensions();


  ////
  //// initiate model by choosing which one and supplying the proper parameters
  //// (do this in a temporary structure inside the constructor of Model and then initialize const arrays from that one)
  ////
  array<bool,     numComps> isBosonic;
  array<bool,     numComps> isCanonical;
  array<unsigned, numComps> Ns;
  array<double,   numComps> maxNsDiff;
  array<bool,     numComps> muOptimization;
  array<double,   numComps> avgNsDiff;
  if (parameters["model"]["equivalent-components"]) {
    isBosonic.fill(     parameters["model"]["is-bosonic"][0]);
    isCanonical.fill(   parameters["model"]["is-canonical"][0]);
    Ns.fill(            parameters["model"]["Ns"][0]);
    maxNsDiff.fill(     parameters["model"]["maxNsDiff"][0]);
    muOptimization.fill(parameters["model"]["muOptimization"][0]);
    avgNsDiff.fill(     parameters["model"]["avgNsDiff"][0]);
  } else {
    copy_n(parameters["model"]["is-bosonic"].begin(),     numComps, isBosonic.begin());
    copy_n(parameters["model"]["is-canonical"].begin(),   numComps, isCanonical.begin());
    copy_n(parameters["model"]["Ns"].begin(),             numComps, Ns.begin());
    copy_n(parameters["model"]["maxNsDiff"].begin(),      numComps, maxNsDiff.begin());
    copy_n(parameters["model"]["muOptimization"].begin(), numComps, muOptimization.begin());
    copy_n(parameters["model"]["avgNsDiff"].begin(),      numComps, avgNsDiff.begin());
  }

  ////
  //// ensure that "maxNsDiff, avgNsDiff <= beta" for the canonical components
  ////
  for (unsigned a = 0; a < numComps; a ++) {
    if (isCanonical[a]) {
      maxNsDiff[a] = min(parameters["model"]["beta"].get<double>(), maxNsDiff[a]);
      avgNsDiff[a] = min(parameters["model"]["beta"].get<double>(), avgNsDiff[a]);
    }
  }
  usedParameters["model"]["maxNsDiff"]      = maxNsDiff;
  usedParameters["model"]["muOptimization"] = muOptimization;
  usedParameters["model"]["avgNsDiff"]      = avgNsDiff;

  Model model{parameters["model"]["name"],
              isBosonic,
              isCanonical,
              Ns,
              parameters["model"]["ts"],
              parameters["model"]["Js"],
              parameters["model"]["Us"],
              parameters["model"]["mus"],
              parameters["model"]["etas"],
              parameters["model"]["beta"],
              lattice.getNumSites(),
              lattice.getNumDimensions(),
              parameters["model"]["uniform"],
              parameters["model"]["equivalent-components"],
              parameters["model"]["maxNumWorms"],
              pseudoRandom};


  ////
  //// initiate Hamiltonian
  ////
  Hamiltonian H{lattice, model};
  usedParameters["model"]["t_ij"]       = H.gett_ija();
  usedParameters["model"]["J_ij"]       = H.getJ_ijab();
  usedParameters["model"]["mu_ia"]      = H.getMu_ia();
  usedParameters["model"]["eta_a"]      = H.getEta_a();
  usedParameters["model"]["U_a"]        = H.getU_a();
  usedParameters["model"]["U_inter_ab"] = H.getU_inter_ab();
  usedParameters["model"]["V_ia"]       = H.getV_ia();
  usedParameters["model"]["U_ab"]       = H.getU_ab();
  usedParameters["model"]["U_nn_ab"]    = H.getU_nn_ab();


  ////
  //// write the used parameters to file as a prettified json object
  ////
  stringstream ss;
  ss << setw(2) << usedParameters << endl;
  writeToFile.aString(ss.str(), "used-parameters.json");


  ////
  //// initiate worm
  ////
  Worm worm{H,
            lattice,
            model,
            pseudoRandom,
            writeToFile,
            parameters["samplingOptions"],
            parameters["computeOptions"],
            parameters["saveOptions"],
            maxNsDiff,
            muOptimization,
            avgNsDiff,
            parameters["development"]};

  ////
  //// warm up or load configuration
  ////
  if (configuration.empty()) {

    // warm up by performing a simulated annealing
    worm.warmUp(numSecWarmUp, 0, numSecOptChem);

    // update used parameters in case of modified chemical potential
    usedParameters["model"]["mu_ia"] = H.getMu_ia();
    usedParameters["model"]["V_ia"]  = H.getV_ia();
    {
      stringstream ss;
      ss << setw(2) << usedParameters << endl;
      writeToFile.aString(ss.str(), "used-parameters.json");
    }
  } else {
    worm.loadConfiguration(configuration);

    const double beta_target = parameters["model"]["beta"].get<double>();
    const double beta_beg    = configuration["beta"].get<double>();

    // warm up by performing a simulated annealing
    if (beta_beg != beta_target) {
      worm.warmUp(numSecWarmUp, beta_beg, numSecOptChem);
    }

  }

  // optimize chemical potential, perform pre-simulation and the actual simulation
  worm.multipleUpdates(numSecPreSim, numSecSim, numData);

  ////
  //// display elapsed time
  ////
  printf("[Program finished in %.2fs]\n", startTime.toc());
  fflush(stdout);

  return 0;
}