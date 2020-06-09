#include "Worm.h"

using namespace std;

void Worm::update ( ) {
  // PROFILING
  // profiler("Worm::update");

  // update counter
  this->currUpdateCount++;

  if (verbose) {
    cout << settings::cout::enterBlue
         << "---------------------" << this->currUpdateCount << "---------------------" << endl
         << settings::cout::resetStyle;
  }

  if (debug && debugFrom == this->currUpdateCount) {
    cout << settings::cout::enterBlue
         << "[starting to debug]" << endl
         << settings::cout::resetStyle;
  }


  ////
  //// decide whether or not to try and insert a new worm
  ////
  const auto insertProbability = this->getInsertProposalProbability(this->actComps);
  bool tryInsert = false;
  if (
    insertProbability == 1 ||
    (insertProbability > 0 && insertProbability > this->pseudoRandom.template U<double>(0, 1))
  ) {
    tryInsert = true;
  }

  // cout << "Worm::update: tryInsert=" << tryInsert << endl;


  ////
  //// perform update procedure
  ////
  if (tryInsert) {
    // Z sector
    this->insert(this->beta, insertProbability, false);
  } else {
    // G sector

    // choose worm uniformly
    const bool wormIndex = this->numWorms == 1 ?
                           0 :
                           this->pseudoRandom.template Uint<unsigned>(0, this->numWorms - 1);

    // choose which update procedure
    unsigned updateIndex = this->pseudoRandom.template Uint<unsigned>(0, this->Wtot - 1);

    (this->*this->weightedUpdateProcedures[updateIndex])(this->beta, wormIndex, insertProbability, false);
  }


  ////
  //// try verify worm configuration and quantities
  ////
  if (
    (debug_major && debugFrom <= this->currUpdateCount) ||
    this->currUpdateCount % verifyEvery == 0
  ) {
    this->checkWormConfiguration();
    this->checkWormQuantities();
  }


  ////
  //// increment counters
  ////
  if (this->numWorms) {
    this->currGsCount++;
  } else {
    this->currZsCount++;
  }


  ////
  //// worm counter
  ////
  if (allowMultiCompWorm) { cout << "Worm::update: implement me!" << endl; }
  for (auto a : this->actComps) this->wormCount[a]++;


  ////
  ////
  ////
  if (allowMultiCompWorm) { cout << "Worm::update: implement me!" << endl; }
  if (this->numWorms == 1) {
    if (this->muOptimization[this->actComps[0]]) {
      const long diff = ((long) this->numParticles[this->actComps[0]]) - this->Ns[this->actComps[0]] * tMax;
      this->avgParticleDiff[this->actComps[0]] += diff;
      this->avgParticleDiff_count[this->actComps[0]]++;
    }
  } else if (this->numWorms == 2) {
    if (this->muOptimization[this->actComps[0]]) {
      const long diff = ((long) this->numParticles[this->actComps[0]]) - this->Ns[this->actComps[0]] * tMax;
      this->avgParticleDiff[this->actComps[0]] += diff;
      this->avgParticleDiff_count[this->actComps[0]]++;
    }
    if (this->muOptimization[this->actComps[1]] && this->actComps[0] != this->actComps[1]) {
      const long diff = ((long) this->numParticles[this->actComps[1]]) - this->Ns[this->actComps[1]] * tMax;
      this->avgParticleDiff[this->actComps[1]] += diff;
      this->avgParticleDiff_count[this->actComps[1]]++;
    }
  }


  ////
  //// print some data every S:th second
  ////
  if (this->currUpdateCount % numUpdatesPerTimeCheck == 0) {
    if (this->printTimer.toc() > printInterval) {
      // print some shit
      printf("[#G=%llu, #Z=%llu, #U=%llu, ΔS=%.2fs]\n",  this->currGsCount, this->currZsCount, this->currUpdateCount, this->printTimer.toc());
      // printf("[#G+ = %llu, #G- = %llu -> (#G+ - #G-)/(#G+ + #G-) = %.6f]\n",
      //        this->numGsCount[3],
      //        this->numGsCount[2],
      //        (this->numGsCount[3] - (double) this->numGsCount[2]) / (double) (this->numGsCount[3] + this->numGsCount[2]));


      ////
      //// TEMP: store also the configuration
      ////
      this->saveConfiguration(this->beta);

      // reset timer
      this->printTimer.reset();
    }
  }
}