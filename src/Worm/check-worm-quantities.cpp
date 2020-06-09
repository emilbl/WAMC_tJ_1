#include "Worm.h"

using namespace std;

void Worm::checkWormQuantities () const {
  // PROFILING
  // profiler("Worm::checkWormQuantities");

  if ( ! this->numWorms) {
    // Z-sector

    ////
    //// for the canonical components, make sure the particle number is conserved
    ////
    for (unsigned a = 0; a < numComps; a++) {
      if (this->isCanonical[a]) {
        const long diff = this->numParticles[a] - tMax * this->Ns[a];
        if (diff != 0) {
          cout << settings::cout::enterRed
               << "Worm::checkWormQuantities: particle number of canonical component " << a << " != " << this->Ns[a]
               << settings::cout::resetStyle << endl;
          if (settings::mode::shutItDown) this->shutDown();
        }
      }
    }

    // total ingoing number of particles and tot number of particles
    array<unsigned long long, numComps> totIngoingNumParticles = {},
                                        totNumParticlesPerSiteSummed = {};
    for (unsigned i = 0; i < this->numSites; i++) {
      for (unsigned a = 0; a < numComps; a++) {
        totIngoingNumParticles[a] += this->sites[i][0]->pop[a];
        totNumParticlesPerSiteSummed[a] += this->numParticlesAtSite[i * numComps + a];
      }
    }

    // check so that the particle number is exactly integer valued
    // need to be modified once a spin-orbit coupling is introduced
    for (unsigned a = 0; a < numComps; a++) {
      if (tMax * totIngoingNumParticles[a] != totNumParticlesPerSiteSummed[a]) {
        cout << settings::cout::enterRed << "Worm::checkWormQuantities: particle number in component " << a << " not integer valued in Z-sector (" << tMax * totIngoingNumParticles[a] - totNumParticlesPerSiteSummed[a] << ")" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
    }
  } else {
    // G-sector

    ////
    //// for the canonical components, make sure the particle number is within +- min(tMax, maxNsDiff / beta_target * tMax)
    ////
    for (unsigned a = 0; a < numComps; a++) {
      if (this->isCanonical[a]) {
        const long diff = this->numParticles[a] - tMax * this->Ns[a];
        const long diffMax = min(tMax, (long) (this->maxNsDiff[a] / this->beta_target * tMax));
        if (abs(diff) > diffMax) {
          cout << settings::cout::enterRed
               << "Worm::checkWormQuantities: particle number of canonical component " << a << " (" << this->numParticles[a] / (double) tMax
               << ") differs from " << this->Ns[a] << " by more than " << diffMax / (double) tMax
               << settings::cout::resetStyle << endl;
          if (settings::mode::shutItDown) this->shutDown();
        }
      }
    }

  }

  ////
  //// for the canonical components, make sure the particle number is within ±tMax
  ////
  for (unsigned a = 0; a < numComps; a++) {
    if (this->isCanonical[a]) {
      const long diff = this->numParticles[a] - tMax * this->Ns[a];
      if (abs(diff) >= tMax) {
        cout << settings::cout::enterRed
             << "Worm::checkWormQuantities: particle number of canonical component " << a << " out of bound (diff = " << diff << ")."
             << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
    }
  }


  ////
  //// compute the worm quantities
  ////
  unsigned                             actualNumKinks;
  double                               actualU_nn;
  array<unsigned, numComps>            actualNumJumps;
  array<unsigned, numExternalInters>   actualNumExchanges;
  array<unsigned long long, numComps>  actualNumParticles;
  array<unsigned long long, numInters> actualNumParticlesSquared;
  array<vector<int>, numComps>         actualNumWinds;
  vector<double>                       actualFlow;
  vector<unsigned long long>           actualNumParticlesAtSite;
  this->computeWormQuantities(actualNumKinks,
                              actualU_nn,
                              actualNumJumps,
                              actualNumExchanges,
                              actualNumParticles,
                              actualNumParticlesSquared,
                              actualNumWinds,
                              actualFlow,
                              actualNumParticlesAtSite);



  // check number of jumps
  if (actualNumJumps != this->numJumps) {
    cout << settings::cout::enterRed << "Worm::checkWormQuantities: incorrect value of this->numJumps:" << endl
         << "actual:  " << actualNumJumps << endl
         << "current: " << this->numJumps << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }

  // check number of kinks
  if (actualNumKinks != this->numKinks) {
    cout << settings::cout::enterRed << "Worm::checkWormQuantities: incorrect value of this->numKinks:" << endl
         << "actual:  " << actualNumKinks << endl
         << "current: " << this->numKinks << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }

  // check number of exchanges
  if (has_J) {
    if (actualNumExchanges != this->numExchanges) {
      cout << settings::cout::enterRed << "Worm::checkWormQuantities: incorrect value of this->numExchanges:" << endl
           << "actual:  " << actualNumExchanges << endl
           << "current: " << this->numExchanges << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) this->shutDown();
    }
  }

  // check number of particles per site
  if (actualNumParticlesAtSite != this->numParticlesAtSite) {
    cout << settings::cout::enterRed << "Worm::checkWormQuantities: incorrect value of this->numParticlesAtSite:" << endl
         << "actual:  " << actualNumParticlesAtSite << endl
         << "current: " << this->numParticlesAtSite << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }

  // check number of particles
  if (actualNumParticles != this->numParticles) {
    cout << settings::cout::enterRed << "Worm::checkWormQuantities: incorrect value of this->numParticles:" << endl
         << "actual:  " << actualNumParticles << endl
         << "current: " << this->numParticles << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }

  // check number of particle squared
  if (actualNumParticlesSquared != this->numParticlesSquared) {
    cout << settings::cout::enterRed << "Worm::checkWormQuantities: incorrect value of this->numParticlesSquared:" << endl
         << "actual:  " << actualNumParticlesSquared << endl
         << "current: " << this->numParticlesSquared << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }

  // check nearest neighbor interaction
  if (has_U_nn) {
    if (abs(actualU_nn - this->U_nn) > diffTresh) {
      cout << settings::cout::enterRed
           << "Worm::checkWormQuantities: incorrect value of this->U_nn:" << endl
           << abs(actualU_nn - this->U_nn) << " = abs(actualU_nn - this->U_nn) > diffTresh = " << diffTresh << endl
           << "actual:  " << actualU_nn << endl
           << "current: " << this->U_nn << endl
           << settings::cout::resetStyle;
      if (settings::mode::shutItDown) this->shutDown();
    }
  }

  // check winding number
  if (actualNumWinds != this->numWinds) {
    cout << settings::cout::enterRed << "Worm::checkWormQuantities: incorrect value of this->numWinds:" << endl
         << "actual:  " << actualNumWinds << endl
         << "current: " << this->numWinds << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }

  // check flow
  if (actualFlow != this->flow) {
    cout << settings::cout::enterRed << "Worm::checkWormQuantities: incorrect value of this->flow:" << endl
         << "actual:  " << actualFlow << endl
         << "current: " << this->flow << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }

}