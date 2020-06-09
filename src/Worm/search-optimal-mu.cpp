#include "Worm.h"

using namespace std;


void Worm::searchOptimalMu (
  const double &        beta,
  const unsigned        comp,
  double &              mu_lower,
  double &              mu_upper,
  const double &        simulationTime,
  const long unsigned & numUpdatesPerTimeCheck,
  const unsigned        numIntervalReductions,
  long unsigned &       numZcounter,
  long unsigned &       numBinnedZcounter,
  const long unsigned & maxNumBinsPerBeta,
  const unsigned        samplingPeriod
) {
  ////
  //// instead, lets do a binary search for the intersection with y=0
  ////
  const auto InduvidualSimulationTime = simulationTime / (numIntervalReductions + 2);

  // get value of initial boundaries
  auto val_lower = this->testMu(beta, comp, mu_lower, numUpdatesPerTimeCheck, InduvidualSimulationTime, numZcounter, numBinnedZcounter, maxNumBinsPerBeta, samplingPeriod),
       val_upper = this->testMu(beta, comp, mu_upper, numUpdatesPerTimeCheck, InduvidualSimulationTime, numZcounter, numBinnedZcounter, maxNumBinsPerBeta, samplingPeriod);


  ////
  ////  val
  ////
  ////   |
  ////   |                                         ___
  ////   |                                    ____|
  ////   |                                 __|
  ////   |                           _____|
  ////   |                       ___|
  ////   |------------------------------------------> mu
  ////   |                 __|
  ////   |          ______|
  ////   |      ___|
  ////   |_____|
  ////   |
  ////
  ////


  for (unsigned i = 0; i < numIntervalReductions; i++) {

    if (false) {
      cout << (i + 1) << " of " << numIntervalReductions
           << ": a = " << comp << ":  mu_lims = [" << mu_lower << ", " << mu_upper << "] => val = [" << val_lower << ", " << val_upper << "]"
           << endl;
    }



    ////
    //// first make sure the intersection lies within the interval
    ////
    if (val_lower > 0 && val_upper > 0) {
      // more Gs with +1 than Gs with -1 particles
      // -> intersection occurs for smaller mu
      const double intervalSize = max(mu_upper - mu_lower, 1e-10);  // in case the interval size is zero
      mu_lower  -= intervalSize;
      val_lower = this->testMu(beta, comp, mu_lower, numUpdatesPerTimeCheck, InduvidualSimulationTime, numZcounter, numBinnedZcounter, maxNumBinsPerBeta, samplingPeriod);
      continue;
    }
    if (val_lower < 0 && val_upper < 0) {
      // more Gs with -1 than Gs with +1 particles
      // -> intersection occurs for larger mu
      const double intervalSize = max(mu_upper - mu_lower, 1e-10);  // in case the interval size is zero
      mu_upper  += intervalSize;
      val_upper = this->testMu(beta, comp, mu_upper, numUpdatesPerTimeCheck, InduvidualSimulationTime, numZcounter, numBinnedZcounter, maxNumBinsPerBeta, samplingPeriod);
      continue;
    }


    ////
    //// interval lies in between
    ////
    const auto mu_middle = 0.5 * (mu_lower + mu_upper);
    const auto val_middle = this->testMu(beta, comp, mu_middle, numUpdatesPerTimeCheck, InduvidualSimulationTime, numZcounter, numBinnedZcounter, maxNumBinsPerBeta, samplingPeriod);

    if (val_middle < 0) {
      // val_lower > 0, val_middle, val_upper < 0
      // [mu_lower, mu_upper] -> [mu_lower, mu_middle]
      mu_lower  = mu_middle;
      val_lower = val_middle;
      continue;
    }
    if (val_middle > 0) {
      // val_lower, val_middle > 0, val_upper < 0
      // [mu_lower, mu_upper] -> [mu_middle, mu_upper]
      mu_upper  = mu_middle;
      val_upper = val_middle;
      continue;
    }
  }

  cout << "a = " << comp << ":  mu_lims = [" << mu_lower << ", " << mu_upper << "] => val = [" << val_lower << ", " << val_upper << "]" << endl;
  ////
  //// set the chemical potential to lie in the middle of the interval
  ////
  this->H.setMu(comp, 0.5 * (mu_lower + mu_upper));
}

void Worm::searchOptimalMu (
  const double &        beta,
  const unsigned        comp,
  double &              mu_lower,
  double &              mu_upper,
  const double &        simulationTime,
  const long unsigned & numUpdatesPerTimeCheck,
  const unsigned        numIntervalReductions
) {
  long unsigned numZcounter       = 0;
  long unsigned numBinnedZcounter = 0;
  long unsigned maxNumBinsPerBeta = 0;
  unsigned samplingPeriod         = 1;

  this->searchOptimalMu(beta,
                        comp,
                        mu_lower,
                        mu_upper,
                        simulationTime,
                        numUpdatesPerTimeCheck,
                        numIntervalReductions,
                        numZcounter,
                        numBinnedZcounter,
                        maxNumBinsPerBeta,
                        samplingPeriod);
}


double Worm::testMu (
  const double &        beta,
  const unsigned        comp,
  const double &        mu,
  const long unsigned & numUpdatesPerTimeCheck,
  const double &        simulationTime,
  long unsigned &       numZcounter,
  long unsigned &       numBinnedZcounter,
  const long unsigned & maxNumBinsPerBeta,
  const unsigned        samplingPeriod
) {
  // set beta
  this->setBeta(beta);

  // set mu
  this->H.setMu(comp, mu);


  // simulate the specified amount of time
  Tic tic{};
  do {
    for (long unsigned i = 0; i < numUpdatesPerTimeCheck; i++) {
      this->update();

      // try to bin every "samplingPeriod"-th Z configurations but no more than "maxNumBinsPerBeta" times
      if (
        ! this->numWorms &&
        ! (numZcounter % samplingPeriod) &&
        numBinnedZcounter < maxNumBinsPerBeta
      ) {
        this->trySampleData(beta);
        numBinnedZcounter++;
      }

      // increase the counter if the current configuration is a Z
      if ( ! this->numWorms) numZcounter++;
    }
  } while (tic.toc() < simulationTime);

  // compute the relative difference in visited Greens function states
  double avgParticleDiff, relDiff;
  if (this->H.model.hasEquivalentComponents()) {
    avgParticleDiff = accumulate(this->avgParticleDiff.begin(), this->avgParticleDiff.end(), 0.0)
                    / accumulate(this->avgParticleDiff_count.begin(), this->avgParticleDiff_count.end(), 0.0)
                    * beta / tMax;
    relDiff = (avgParticleDiff - this->targetAvgNsDiff[0]) / this->maxNsDiff[0];
  } else {
    avgParticleDiff = this->avgParticleDiff[comp] / this->avgParticleDiff_count[comp] * beta / tMax;
    relDiff = (avgParticleDiff - this->targetAvgNsDiff[comp]) / this->maxNsDiff[comp];
  }

  cout << "----------------" << endl;
  cout << "mu              = " << mu << endl;
  cout << "relDiff         = " << relDiff << endl;
  cout << "avgParticleDiff = " << avgParticleDiff << endl;
  cout << "----------------" << endl;

  // reset
  this->avgParticleDiff.fill(0);
  this->avgParticleDiff_count.fill(0);

  return relDiff;
}