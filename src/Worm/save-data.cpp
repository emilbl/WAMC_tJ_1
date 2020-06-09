#include "Worm.h"

using namespace std;
using json = nlohmann::json;


void Worm::saveData (
  const unsigned long long & prevNumData,
  const unsigned long long & newNumData,
  const double &             beta,
  const string &             suffix
) const {
  // PROFILING
  // profiler("Worm::saveData");

  ////
  //// store binary data
  ////
  this->saveBinaryFiles(prevNumData, newNumData, suffix);


  ////
  //// statistics
  ////
  json statistics;
  json statistics_raw;
  statistics["version"] = 17;

  // add simulated annealing statistics
  statistics["warmup"] = this->warmupStatistics;

  // update procedure acceptance ratios
  statistics["updates"] = this->updateStatistics;

  // average segments length for each component
  vector<vector<double> > container;
  this->calculateAverageSegmentLengths(beta, container);
  statistics["segLengHist"] = container;


  ////
  //// some useful quantities
  ////
  const auto numData      = newNumData;
  const auto numLargeData = newNumData / (double) this->largeDataSamplingPeriod;


  ////
  //// some other parameters
  ////
  statistics["numData"]                 = newNumData;
  statistics["numLargeData"]            = numLargeData;
  statistics["dataSamplingPeriod"]      = this->samplingPeriod;
  statistics["largeDataSamplingPeriod"] = largeDataSamplingPeriod;


  ////
  //// information about greens function
  ////
  if (this->sample_GreensFunction) {
    statistics["G"] = json{
      {"N_t", G_hist_N}
    };
  }

  ////
  //// information about the four point correlator
  ////
  if (this->sample_fourPointCorrelator) {

    // which components are being sampled
    vector<unsigned> components;
    for (unsigned a = 0; a < numComps; a++) {
      if (this->fpc_hist[a].size()) components.push_back(a);
    }

    statistics["fpc"] = json{
      {"components", components},
      {"N_t",        fpc_hist_N}
    };
  }

  ////
  //// information about the density-density correlator
  ////
  if (this->sample_densityDensityCorrelator) {

    // which components are being sampled
    vector<unsigned> components;
    for (unsigned a = 0; a < numComps; a++) {
      if (this->ddc_hist[a].size()) components.push_back(a);
    }

    statistics["ddc"] = json{
      {"components", components},
      {"N_o",        ddc_hist_N_o},
      {"N_d",        ddc_hist_N_d}
    };
  }




  ////
  //// compress sampled quantities
  ////
  this->compressData(statistics);


  ////
  //// computed quantities
  ////

  // average sign
  const double avgSign = ((double) this->posSignCount - (double) this->negSignCount) / (this->posSignCount + this->negSignCount);
  statistics["avgSign"] = avgSign;

  if (this->compute_numParticles) {
    // single
    vector<double> avgNumParticlesOne_normalized(this->avgNumParticlesOne.size() / 4);
    for (unsigned i = 0; i < avgNumParticlesOne_normalized.size(); i++) {
      avgNumParticlesOne_normalized[i] = (
           (double) this->avgNumParticlesOne[0 * avgNumParticlesOne_normalized.size() + i]
         - (double) this->avgNumParticlesOne[1 * avgNumParticlesOne_normalized.size() + i]
         - (double) this->avgNumParticlesOne[2 * avgNumParticlesOne_normalized.size() + i]
         + (double) this->avgNumParticlesOne[3 * avgNumParticlesOne_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgNumParticlesOne"] = avgNumParticlesOne_normalized;
    statistics_raw["avgNumParticlesOne"] = this->avgNumParticlesOne;

    // pairwise
    vector<double> avgNumParticlesTwo_normalized(this->avgNumParticlesTwo.size() / 4);
    for (unsigned i = 0; i < avgNumParticlesTwo_normalized.size(); i++) {
      avgNumParticlesTwo_normalized[i] = (
           (double) this->avgNumParticlesTwo[0 * avgNumParticlesTwo_normalized.size() + i]
         - (double) this->avgNumParticlesTwo[1 * avgNumParticlesTwo_normalized.size() + i]
         - (double) this->avgNumParticlesTwo[2 * avgNumParticlesTwo_normalized.size() + i]
         + (double) this->avgNumParticlesTwo[3 * avgNumParticlesTwo_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgNumParticlesTwo"] = avgNumParticlesTwo_normalized;
    statistics_raw["avgNumParticlesTwo"] = this->avgNumParticlesTwo;

    // all
    vector<double> avgNumParticlesAll_normalized(this->avgNumParticlesAll.size() / 4);
    for (unsigned i = 0; i < avgNumParticlesAll_normalized.size(); i++) {
      avgNumParticlesAll_normalized[i] = (
           (double) this->avgNumParticlesAll[0 * avgNumParticlesAll_normalized.size() + i]
         - (double) this->avgNumParticlesAll[1 * avgNumParticlesAll_normalized.size() + i]
         - (double) this->avgNumParticlesAll[2 * avgNumParticlesAll_normalized.size() + i]
         + (double) this->avgNumParticlesAll[3 * avgNumParticlesAll_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgNumParticlesAll"] = avgNumParticlesAll_normalized;
    statistics_raw["avgNumParticlesAll"] = this->avgNumParticlesAll;
  }

  if (this->compute_numWinds) {
    vector<double> avgNumWinds_normalized(this->avgNumWinds.size() / 4);
    for (unsigned i = 0; i < avgNumWinds_normalized.size(); i++) {
      avgNumWinds_normalized[i] = (
           (double) this->avgNumWinds[0 * avgNumWinds_normalized.size() + i]
         - (double) this->avgNumWinds[1 * avgNumWinds_normalized.size() + i]
         - (double) this->avgNumWinds[2 * avgNumWinds_normalized.size() + i]
         + (double) this->avgNumWinds[3 * avgNumWinds_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgNumWinds"] = avgNumWinds_normalized;
    statistics_raw["avgNumWinds"] = this->avgNumWinds;

    // cyclic
    vector<double> avgNumWindsCyclic_normalized(this->avgNumWindsCyclic.size() / 4);
    for (unsigned i = 0; i < avgNumWindsCyclic_normalized.size(); i++) {
      avgNumWindsCyclic_normalized[i] = (
           (double) this->avgNumWindsCyclic[0 * avgNumWindsCyclic_normalized.size() + i]
         - (double) this->avgNumWindsCyclic[1 * avgNumWindsCyclic_normalized.size() + i]
         - (double) this->avgNumWindsCyclic[2 * avgNumWindsCyclic_normalized.size() + i]
         + (double) this->avgNumWindsCyclic[3 * avgNumWindsCyclic_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgNumWindsCyclic"] = avgNumWindsCyclic_normalized;
    statistics_raw["avgNumWindsCyclic"] = this->avgNumWindsCyclic;

    ////
    //// counter- and coflow
    ////
    vector<double> avgPairwiseCounterflow_normalized(this->avgPairwiseCounterflow.size() / 4);
    for (unsigned i = 0; i < avgPairwiseCounterflow_normalized.size(); i++) {
      avgPairwiseCounterflow_normalized[i] = (
           (double) this->avgPairwiseCounterflow[0 * avgPairwiseCounterflow_normalized.size() + i]
         - (double) this->avgPairwiseCounterflow[1 * avgPairwiseCounterflow_normalized.size() + i]
         - (double) this->avgPairwiseCounterflow[2 * avgPairwiseCounterflow_normalized.size() + i]
         + (double) this->avgPairwiseCounterflow[3 * avgPairwiseCounterflow_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgPairwiseCounterflow"] = avgPairwiseCounterflow_normalized;
    statistics_raw["avgPairwiseCounterflow"] = this->avgPairwiseCounterflow;


    vector<double> avgCounterflow_normalized(this->avgCounterflow.size() / 4);
    for (unsigned i = 0; i < avgCounterflow_normalized.size(); i++) {
      avgCounterflow_normalized[i] = (
           (double) this->avgCounterflow[0 * avgCounterflow_normalized.size() + i]
         - (double) this->avgCounterflow[1 * avgCounterflow_normalized.size() + i]
         - (double) this->avgCounterflow[2 * avgCounterflow_normalized.size() + i]
         + (double) this->avgCounterflow[3 * avgCounterflow_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgCounterflow"] = avgCounterflow_normalized;
    statistics_raw["avgCounterflow"] = this->avgCounterflow;


    vector<double> avgCoCounterflow_normalized(this->avgCoCounterflow.size() / 4);
    for (unsigned i = 0; i < avgCoCounterflow_normalized.size(); i++) {
      avgCoCounterflow_normalized[i] = (
           (double) this->avgCoCounterflow[0 * avgCoCounterflow_normalized.size() + i]
         - (double) this->avgCoCounterflow[1 * avgCoCounterflow_normalized.size() + i]
         - (double) this->avgCoCounterflow[2 * avgCoCounterflow_normalized.size() + i]
         + (double) this->avgCoCounterflow[3 * avgCoCounterflow_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgCoCounterflow"] = avgCoCounterflow_normalized;
    statistics_raw["avgCoCounterflow"] = this->avgCoCounterflow;


    vector<double> avgCoflow_normalized(this->avgCoflow.size() / 4);
    for (unsigned i = 0; i < avgCoflow_normalized.size(); i++) {
      avgCoflow_normalized[i] = (
           (double) this->avgCoflow[0 * avgCoflow_normalized.size() + i]
         - (double) this->avgCoflow[1 * avgCoflow_normalized.size() + i]
         - (double) this->avgCoflow[2 * avgCoflow_normalized.size() + i]
         + (double) this->avgCoflow[3 * avgCoflow_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgCoflow"] = avgCoflow_normalized;
    statistics_raw["avgCoflow"] = this->avgCoflow;
  }

  if (this->compute_kinetEnergy) {
    vector<double> avgKinetEnergy_normalized(this->avgKinetEnergy.size() / 4);
    for (unsigned i = 0; i < avgKinetEnergy_normalized.size(); i++) {
      avgKinetEnergy_normalized[i] = (
           (double) this->avgKinetEnergy[0 * avgKinetEnergy_normalized.size() + i]
         - (double) this->avgKinetEnergy[1 * avgKinetEnergy_normalized.size() + i]
         - (double) this->avgKinetEnergy[2 * avgKinetEnergy_normalized.size() + i]
         + (double) this->avgKinetEnergy[3 * avgKinetEnergy_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgKinetEnergy"] = avgKinetEnergy_normalized;
    statistics_raw["avgKinetEnergy"] = this->avgKinetEnergy;
  }

  if (has_J && this->compute_exchaEnergy) {
    vector<double> avgExchaEnergy_normalized(this->avgExchaEnergy.size() / 4);
    for (unsigned i = 0; i < avgExchaEnergy_normalized.size(); i++) {
      avgExchaEnergy_normalized[i] = (
           (double) this->avgExchaEnergy[0 * avgExchaEnergy_normalized.size() + i]
         - (double) this->avgExchaEnergy[1 * avgExchaEnergy_normalized.size() + i]
         - (double) this->avgExchaEnergy[2 * avgExchaEnergy_normalized.size() + i]
         + (double) this->avgExchaEnergy[3 * avgExchaEnergy_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgExchaEnergy"] = avgExchaEnergy_normalized;
    statistics_raw["avgExchaEnergy"] = this->avgExchaEnergy;
  }

  if (this->compute_potenEnergy) {
    vector<double> avgPotenEnergy_normalized(this->avgPotenEnergy.size() / 4);
    for (unsigned i = 0; i < avgPotenEnergy_normalized.size(); i++) {
      avgPotenEnergy_normalized[i] = (
           (double) this->avgPotenEnergy[0 * avgPotenEnergy_normalized.size() + i]
         - (double) this->avgPotenEnergy[1 * avgPotenEnergy_normalized.size() + i]
         - (double) this->avgPotenEnergy[2 * avgPotenEnergy_normalized.size() + i]
         + (double) this->avgPotenEnergy[3 * avgPotenEnergy_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgPotenEnergy"] = avgPotenEnergy_normalized;
    statistics_raw["avgPotenEnergy"] = this->avgPotenEnergy;
  }

  if (has_U && this->compute_interEnergy) {
    vector<double> avgInterEnergy_normalized(this->avgInterEnergy.size() / 4);
    for (unsigned i = 0; i < avgInterEnergy_normalized.size(); i++) {
      avgInterEnergy_normalized[i] = (
           (double) this->avgInterEnergy[0 * avgInterEnergy_normalized.size() + i]
         - (double) this->avgInterEnergy[1 * avgInterEnergy_normalized.size() + i]
         - (double) this->avgInterEnergy[2 * avgInterEnergy_normalized.size() + i]
         + (double) this->avgInterEnergy[3 * avgInterEnergy_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgInterEnergy"] = avgInterEnergy_normalized;
    statistics_raw["avgInterEnergy"] = this->avgInterEnergy;
  }

  if (has_U_nn && this->compute_nnInterEnergy) {
    vector<double> avgNnInterEnergy_normalized(this->avgNnInterEnergy.size() / 4);
    for (unsigned i = 0; i < avgNnInterEnergy_normalized.size(); i++) {
      avgNnInterEnergy_normalized[i] = (
           (double) this->avgNnInterEnergy[0 * avgNnInterEnergy_normalized.size() + i]
         - (double) this->avgNnInterEnergy[1 * avgNnInterEnergy_normalized.size() + i]
         - (double) this->avgNnInterEnergy[2 * avgNnInterEnergy_normalized.size() + i]
         + (double) this->avgNnInterEnergy[3 * avgNnInterEnergy_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgNnInterEnergy"] = avgNnInterEnergy_normalized;
    statistics_raw["avgNnInterEnergy"] = this->avgNnInterEnergy;
  }

  if (this->compute_osParticleProd) {
    vector<double> osParticleProd_normalized(this->osParticleProd.size() / 4);
    for (unsigned i = 0; i < osParticleProd_normalized.size(); i++) {
      osParticleProd_normalized[i] = (
           (double) this->osParticleProd[0 * osParticleProd_normalized.size() + i]
         - (double) this->osParticleProd[1 * osParticleProd_normalized.size() + i]
         - (double) this->osParticleProd[2 * osParticleProd_normalized.size() + i]
         + (double) this->osParticleProd[3 * osParticleProd_normalized.size() + i]
      ) / numData / avgSign;
    }
    statistics["avgOsParticleProd"] = osParticleProd_normalized;
    statistics_raw["avgOsParticleProd"] = this->osParticleProd;
  }


  if (this->compute_nnParticleProd) {
    vector<double> nnParticleProd_normalized(this->nnParticleProd.size() / 4);
    for (unsigned i = 0; i < nnParticleProd_normalized.size(); i++) {
      nnParticleProd_normalized[i] = (
           (double) this->nnParticleProd[0 * nnParticleProd_normalized.size() + i]
         - (double) this->nnParticleProd[1 * nnParticleProd_normalized.size() + i]
         - (double) this->nnParticleProd[2 * nnParticleProd_normalized.size() + i]
         + (double) this->nnParticleProd[3 * nnParticleProd_normalized.size() + i]
      ) / numData / avgSign
      *
      0.5;   // to remove double counting
    }
    statistics["avgNnParticleProd"] = nnParticleProd_normalized;
    statistics_raw["avgNnParticleProd"] = this->nnParticleProd;
  }


  if (this->compute_nnnParticleProd) {
    vector<double> nnnParticleProd_normalized(this->nnnParticleProd.size() / 4);
    for (unsigned i = 0; i < nnnParticleProd_normalized.size(); i++) {
      nnnParticleProd_normalized[i] = (
           (double) this->nnnParticleProd[0 * nnnParticleProd_normalized.size() + i]
         - (double) this->nnnParticleProd[1 * nnnParticleProd_normalized.size() + i]
         - (double) this->nnnParticleProd[2 * nnnParticleProd_normalized.size() + i]
         + (double) this->nnnParticleProd[3 * nnnParticleProd_normalized.size() + i]
      ) / numData / avgSign
      *
      0.5;   // to remove double counting
    }
    statistics["avgNnnParticleProd"] = nnnParticleProd_normalized;
    statistics_raw["avgNnnParticleProd"] = this->nnnParticleProd;
  }




  ////
  //// average nearest neighbor particle number correlations
  //// and average particle numbers, both centered about the hole
  ////
  statistics_raw["numLargeData"] = numLargeData;

  if (modelType == tJ && this->compute_C1 && this->Ns[1] == 1) {
    vector<double> avgC1_normalized(this->avgC1.size() / 4);
    for (unsigned i = 0; i < avgC1_normalized.size(); i++) {
      avgC1_normalized[i] = (
           (double) this->avgC1[0 * avgC1_normalized.size() + i]
         - (double) this->avgC1[1 * avgC1_normalized.size() + i]
         - (double) this->avgC1[2 * avgC1_normalized.size() + i]
         + (double) this->avgC1[3 * avgC1_normalized.size() + i]
      ) / numLargeData / avgSign;
    }
    statistics["C1"] = avgC1_normalized;
    statistics_raw["C1"] = this->avgC1;
  }

  if (modelType == tJ && this->compute_C2 && this->Ns[1] == 1) {
    vector<double> avgC2_normalized(this->avgC2.size() / 4);
    for (unsigned i = 0; i < avgC2_normalized.size(); i++) {
      avgC2_normalized[i] = (
           (double) this->avgC2[0 * avgC2_normalized.size() + i]
         - (double) this->avgC2[1 * avgC2_normalized.size() + i]
         - (double) this->avgC2[2 * avgC2_normalized.size() + i]
         + (double) this->avgC2[3 * avgC2_normalized.size() + i]
      ) / numLargeData / avgSign;
    }
    statistics["C2"] = avgC2_normalized;
    statistics_raw["C2"] = this->avgC2;
  }

  if (modelType == tJ && this->compute_C3 && this->Ns[1] == 1) {
    vector<double> avgC3_normalized(this->avgC3.size() / 4);
    for (unsigned i = 0; i < avgC3_normalized.size(); i++) {
      avgC3_normalized[i] = (
           (double) this->avgC3[0 * avgC3_normalized.size() + i]
         - (double) this->avgC3[1 * avgC3_normalized.size() + i]
         - (double) this->avgC3[2 * avgC3_normalized.size() + i]
         + (double) this->avgC3[3 * avgC3_normalized.size() + i]
      ) / numLargeData / avgSign;
    }
    statistics["C3"] = avgC3_normalized;
    statistics_raw["C3"] = this->avgC3;
  }

  if (modelType == tJ && this->compute_particleCorr) {

    vector<double> avgParticleCorr_normalized(this->avgParticleCorr.size() / 4);
    for (unsigned i = 0; i < avgParticleCorr_normalized.size(); i++) {
      avgParticleCorr_normalized[i] = (
           (double) this->avgParticleCorr[0 * avgParticleCorr_normalized.size() + i]
         - (double) this->avgParticleCorr[1 * avgParticleCorr_normalized.size() + i]
         - (double) this->avgParticleCorr[2 * avgParticleCorr_normalized.size() + i]
         + (double) this->avgParticleCorr[3 * avgParticleCorr_normalized.size() + i]
      ) / numLargeData / avgSign;
    }
    statistics["avgParticleCorr"] = avgParticleCorr_normalized;

    statistics_raw["avgParticleCorr"] = this->avgParticleCorr;
  }

  if (modelType == tJ && this->compute_particleConv) {

    vector<double> avgParticleConv_normalized(this->avgParticleConv.size() / 4);
    for (unsigned i = 0; i < avgParticleConv_normalized.size(); i++) {
      avgParticleConv_normalized[i] = (
           (double) this->avgParticleConv[0 * avgParticleConv_normalized.size() + i]
         - (double) this->avgParticleConv[1 * avgParticleConv_normalized.size() + i]
         - (double) this->avgParticleConv[2 * avgParticleConv_normalized.size() + i]
         + (double) this->avgParticleConv[3 * avgParticleConv_normalized.size() + i]
      ) / numLargeData / avgSign;
    }
    statistics["avgParticleConv"] = avgParticleConv_normalized;

    statistics_raw["avgParticleConv"] = this->avgParticleConv;
  }


  ////
  //// append raw statistics
  ////
  statistics["raw"] = statistics_raw;

  ////
  //// write prettified json object to file
  ////
  stringstream ss;
  ss << setw(2) << statistics << endl;
  this->writeToFile.aString(ss.str(), "statistics.json");
}