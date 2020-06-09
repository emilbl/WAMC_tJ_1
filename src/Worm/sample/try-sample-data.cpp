#include "../Worm.h"

using namespace std;


void Worm::trySampleData (
  const double & beta
) {
  if (this->numWorms == 0) {
    // increment counter
    this->currBinCount0++;

    // partition function quantities
    this->sampleZ();

  } else if (this->numWorms == 1) {
    // increment counter
    this->currBinCount1++;

    // Greens function
    if (this->sample_GreensFunction && this->isTranslInvariant) this->sampleG();

  } else if (this->numWorms == 2) {
    // increment counter
    this->currBinCount2++;

    // four point correlator
    if (this->sample_fourPointCorrelator) this->sampleFPC();

    // density-density correlator
    if (this->sample_densityDensityCorrelator) this->sampleDDC();

  }
}