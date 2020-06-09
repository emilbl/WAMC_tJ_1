#include "Worm.h"

using namespace std;


void Worm::resetDataContainers () {
  ////
  //// small data
  ////
  fill(this->signHist_small.begin(),    this->signHist_small.end(),    0);
  fill(this->numParticlesHist.begin(),  this->numParticlesHist.end(),  0);
  fill(this->numWindsHist.begin(),      this->numWindsHist.end(),      0);
  fill(this->kinetEnergyHist.begin(),   this->kinetEnergyHist.end(),   0);
  fill(this->potenEnergyHist.begin(),   this->potenEnergyHist.end(),   0);
  fill(this->totEnergyHist.begin(),     this->totEnergyHist.end(),     0);
  fill(this->updateStatistics.begin(),  this->updateStatistics.end(),  0);
  fill(this->exchaEnergyHist.begin(),   this->exchaEnergyHist.end(),   0);
  fill(this->interEnergyHist.begin(),   this->interEnergyHist.end(),   0);
  fill(this->nnInterEnergyHist.begin(), this->nnInterEnergyHist.end(), 0);
  fill(this->osParticleProd.begin(),    this->osParticleProd.end(),    0);
  fill(this->nnParticleProd.begin(),    this->nnParticleProd.end(),    0);
  fill(this->nnnParticleProd.begin(),   this->nnnParticleProd.end(),   0);


  ////
  //// large data
  ////
  fill(this->signHist_large.begin(),             this->signHist_large.end(),             0);
  fill(this->numParticlesAtSiteHist.begin(),     this->numParticlesAtSiteHist.end(),     0);
  fill(this->flowHist.begin(),                   this->flowHist.end(),                   0);
  fill(this->instantParticleNumHist.begin(),     this->instantParticleNumHist.end(),     0);
  fill(this->instantParticleNumConvHist.begin(), this->instantParticleNumConvHist.end(), 0);


  ////
  //// compute containers
  ////
  fill(this->avgNumParticlesOne.begin(), this->avgNumParticlesOne.end(), 0);
  fill(this->avgNumParticlesTwo.begin(), this->avgNumParticlesTwo.end(), 0);
  fill(this->avgNumParticlesAll.begin(), this->avgNumParticlesAll.end(), 0);

  fill(this->avgNumWinds.begin(),       this->avgNumWinds.end(),       0);
  fill(this->avgNumWindsCyclic.begin(), this->avgNumWindsCyclic.end(), 0);

  fill(this->avgPairwiseCounterflow.begin(), this->avgPairwiseCounterflow.end(), 0);
  fill(this->avgCounterflow.begin(),         this->avgCounterflow.end(),         0);
  fill(this->avgCoCounterflow.begin(),       this->avgCoCounterflow.end(),       0);
  fill(this->avgCoflow.begin(),              this->avgCoflow.end(),              0);

  fill(this->avgKinetEnergy.begin(),   this->avgKinetEnergy.end(),   0);
  fill(this->avgExchaEnergy.begin(),   this->avgExchaEnergy.end(),   0);
  fill(this->avgPotenEnergy.begin(),   this->avgPotenEnergy.end(),   0);
  fill(this->avgInterEnergy.begin(),   this->avgInterEnergy.end(),   0);
  fill(this->avgNnInterEnergy.begin(), this->avgNnInterEnergy.end(), 0);

  fill(this->avgC1.begin(), this->avgC1.end(), 0);
  fill(this->avgC2.begin(), this->avgC2.end(), 0);
  fill(this->avgC3.begin(), this->avgC3.end(), 0);
  fill(this->avgParticleCorr.begin(), this->avgParticleCorr.end(), 0);
  fill(this->avgParticleConv.begin(), this->avgParticleConv.end(), 0);




  ////
  //// partition function bin
  ////
  this->Z_bin = 0;
  this->Z_bin_count = 0;


  ////
  //// Green's function histogram
  ////
  for (auto & G : this->G_hist) {
    fill(G.begin(), G.end(), 0);
  }
  for (auto & G : this->G_hist_count) {
    fill(G.begin(), G.end(), 0);
  }

  ////
  //// four point correlator histograms
  ////
  for (auto & h : this->fpc_hist) {
    fill(h.begin(), h.end(), 0);
  }
  for (auto & h : this->fpc_hist_count) {
    fill(h.begin(), h.end(), 0);
  }

  ////
  //// density-density correlator histograms
  ////
  for (auto & h : this->ddc_hist) {
    fill(h.begin(), h.end(), 0);
  }
  for (auto & h : this->ddc_hist_count) {
    fill(h.begin(), h.end(), 0);
  }

  ////
  //// reset counter
  ////
  this->currBinCount0 = 0;
  this->currBinCount1 = 0;
  this->currBinCount2 = 0;


  ////
  //// reset sign counters
  ////
  this->posSignCount = 0;
  this->negSignCount = 0;

  this->wormCount.fill(0);


  ////
  //// reset other containers
  ////
  this->avgParticleDiff.fill(0);
  this->avgParticleDiff_count.fill(0);
}