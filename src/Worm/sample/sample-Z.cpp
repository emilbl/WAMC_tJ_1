#include "../Worm.h"

using namespace std;


void Worm::sampleZ () {
  // PROFILING
  // profiler("Worm::sampleZ");

  ////
  //// compute additional sign contribution stemming from fermionic exchange
  ////
  int fermionicSign = 1;
  if (this->sample_GreensFunction) {
    // compute fermionic exchange sign
    if (this->hasFermionicExchange) {
      fermionicSign = this->computeFermionicExchangeSign();
    }

    ////
    //// increment normalization bin
    ////
    this->Z_bin       += this->sign * fermionicSign;
    this->Z_bin_count += 1;
  }


  // bin period time
  if (++this->currTriedBinCount / this->samplingPeriod) {
    // reset counter

    if ( ! this->sample_GreensFunction && this->hasFermionicExchange) {
      // has not already computed fermionic exchange sign
      fermionicSign = this->computeFermionicExchangeSign();
    }

    ////
    //// compute fermionic sign shift
    ////
    const unsigned signShift = fermionicSign == 1 ?
                               (this->sign == 1 ? 0 : 1) :
                               (this->sign == 1 ? 2 : 3);


    ////
    //// store sign
    ////
    if (this->sign * fermionicSign == 1) this->posSignCount++;
    else                                 this->negSignCount++;


    ////
    //// compute average particle number
    ////
    if (this->compute_numParticles) {
      for (unsigned a = 0; a < numComps; a++) {
        // single component
        this->avgNumParticlesOne[signShift * numComps + a] += this->numParticles[a] / (double) tMax;

        // pairwise
        for (unsigned b = a; b < numComps; b++) {
          this->avgNumParticlesTwo[signShift * numInters + this->i_ab(a, b)] += (this->numParticles[a] / (double) tMax) * (this->numParticles[b] / (double) tMax);
        }
      }

      // all
      double all = 1;
      for (unsigned a = 0; a < numComps; a++) all *= this->numParticles[a] / (double) tMax;
      this->avgNumParticlesAll[signShift] += all;
    }


    ////
    //// compute average kinetic energy
    ////
    if (this->compute_kinetEnergy) {
      for (unsigned a = 0; a < numComps; a++) {
        this->avgKinetEnergy[signShift * numComps + a] += - (double) this->numJumps[a] / this->beta_target;
      }
    }

    ////
    //// compute average exchange energy
    ////
    if (has_J && this->compute_exchaEnergy) {
      for (unsigned ab = 0; ab < numExternalInters; ab++) {
        this->avgExchaEnergy[signShift * numExternalInters + ab] += - (double) this->numExchanges[ab] / this->beta_target;
      }
    }

    ////
    //// compute average potential energy
    ////
    if (this->compute_potenEnergy) {
      if (this->H.isUniform()) {
        this->H.calculatePotenEnergies_hom(this->numParticles,
                                           this->avgPotenEnergy.begin() + signShift * numComps,
                                           true);
      } else {
        this->H.calculatePotenEnergies_inhom(this->numParticlesAtSite,
                                             this->avgPotenEnergy.begin() + signShift * numComps,
                                             true);
      }
    }

    ////
    //// compute average on-site interaction energy
    ////
    if (has_U && this->compute_interEnergy) {
      this->H.calculateInterEnergies(this->numParticlesSquared,
                                     this->numParticles,
                                     this->avgInterEnergy.begin() + signShift * numInters,
                                     true);
    }

    ////
    //// compute nearest neighbor interaction energy
    ////
    if (has_U_nn && this->compute_nnInterEnergy) {
      this->avgNnInterEnergy[signShift] += this->U_nn;
    }


    ////
    //// compute windings
    ////
    if (compute_numWinds) {
      for (unsigned a = 0; a < numComps; a++) {
        for (unsigned b = a; b < numComps; b++) {
          for (unsigned d = 0; d < this->numDims; d++) {
            this->avgNumWinds[signShift * numInters + this->i_ab(a, b)] += this->numWinds[a][d] * this->numWinds[b][d];
          }
        }
      }

      // cyclic
      for (unsigned a = 0; a < numComps; a++) {
        double cyclicProd = 1;
        for (unsigned shift = 0; shift < numComps; shift++) {
          double innerProd = 0;
          for (unsigned d = 0; d < this->numDims; d++) innerProd += this->numWinds[a][d] * this->numWinds[(a + shift) % numComps][d];
          cyclicProd *= innerProd;
        }
        this->avgNumWindsCyclic[signShift * numComps + a] += cyclicProd;
      }


      ////
      //// counterflow order parameters
      ////
      if (numComps == 3 && this->numDims == 2) {

        ////
        //// need all inner products since if for one component W = 0 it is
        //// impossible to deduce the inner product between the other two
        ////

        // winding scalar products
        const int _w01 =  inner_product(this->numWinds[0].begin(), this->numWinds[0].end(), this->numWinds[1].begin(), 0); //this->numWinds[0][0] * this->numWinds[1][0] +  this->numWinds[0][1] * this->numWinds[1][1];
        const int _w02 =  inner_product(this->numWinds[0].begin(), this->numWinds[0].end(), this->numWinds[2].begin(), 0); //this->numWinds[0][0] * this->numWinds[2][0] +  this->numWinds[0][1] * this->numWinds[2][1];
        const int _w12 =  inner_product(this->numWinds[1].begin(), this->numWinds[1].end(), this->numWinds[2].begin(), 0); //this->numWinds[1][0] * this->numWinds[2][0] +  this->numWinds[1][1] * this->numWinds[2][1];

        // normalize
        const double w0 = sqrt(inner_product(this->numWinds[0].begin(), this->numWinds[0].end(), this->numWinds[0].begin(), 0));
        const double w1 = sqrt(inner_product(this->numWinds[1].begin(), this->numWinds[1].end(), this->numWinds[1].begin(), 0));
        const double w2 = sqrt(inner_product(this->numWinds[2].begin(), this->numWinds[2].end(), this->numWinds[2].begin(), 0));
        const double w01 = _w01 == 0 ? 0 : _w01 / w0 / w1;
        const double w02 = _w02 == 0 ? 0 : _w02 / w0 / w2;
        const double w12 = _w12 == 0 ? 0 : _w12 / w1 / w2;

        ////
        //// store also normalization correlators
        ////
        //// < (  [w - ...]^2   [w - ...]^2   [w - ...]^2  )^3 >
        //// ---------------------------------------------------
        ////   ( <[w - ...]^2> <[w - ...]^2> <[w - ...]^2> )^3
        ////
        #pragma message("store also normalization correlators")

        this->avgCounterflow[signShift] += pow(pow(w01 + 0.5, 2) + pow(w02 + 0.5, 2) + pow(w12 + 0.5, 2), 3);

        this->avgCoflow[signShift] += pow(pow(w01 - 1, 2) + pow(w02 - 1, 2) + pow(w12 - 1, 2), 3);

        this->avgPairwiseCounterflow[signShift] += pow(  (pow(w01 + 1, 2) + pow(w02, 2)     + pow(w12, 2)    )
                                                       * (pow(w01, 2)     + pow(w02 + 1, 2) + pow(w12, 2)    )
                                                       * (pow(w01, 2)     + pow(w02, 2)     + pow(w12 + 1, 2)), 3);

        this->avgCoCounterflow[signShift] += pow(  (pow(w01 - 1, 2) + pow(w02 + 1, 2) + pow(w12 + 1, 2))
                                                 * (pow(w01 + 1, 2) + pow(w02 - 1, 2) + pow(w12 + 1, 2))
                                                 * (pow(w01 + 1, 2) + pow(w02 + 1, 2) + pow(w12 - 1, 2)), 3);

      }
    }


    ////
    //// compute average NNs particle correlation and average particle number centered about the hole
    ////
    if (modelType == tJ && this->compute_spatialCorrelators) {
      this->sampleSpatialCorrelators();
    }


    ////
    //// on-site particle product
    ////
    if (this->compute_osParticleProd) {
      for (const auto & sites : this->sites) {
        const auto & seg = sites[0];
        unsigned i = 0;
        for (unsigned a = 0; a < numComps; a++) {
          for (unsigned b = a; b < numComps; b++) {
            this->osParticleProd[signShift * numInters + i++] += seg->pop[a] * seg->pop[b];
          }
        }
      }
    }


    ////
    //// nearest-neighbor particle product
    ////
    if (this->compute_nnParticleProd) {
      for (const auto & sites : this->sites) {
        const auto & seg = sites[0];
        for (const auto i : this->lattice.NNs[seg->siteIndex]) {
          const auto & _seg = this->sites[i][0];
          unsigned j = 0;
          for (unsigned a = 0; a < numComps; a++) {
            for (unsigned b = a; b < numComps; b++) {
              this->nnParticleProd[signShift * numInters + j++] += seg->pop[a] * _seg->pop[b];
            }
          }
        }
      }
    }


    ////
    //// next nearest-neighbor particle product
    ////
    if (this->compute_nnnParticleProd) {
      for (const auto & sites : this->sites) {
        const auto & seg = sites[0];
        for (const auto i : this->lattice.N2Ns[seg->siteIndex]) {
          const auto & _seg = this->sites[i][0];
          unsigned j = 0;
          for (unsigned a = 0; a < numComps; a++) {
            for (unsigned b = a; b < numComps; b++) {
              this->nnnParticleProd[signShift * numInters + j++] += seg->pop[a] * _seg->pop[b];
            }
          }
        }
      }
    }



    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////// reconstruct distribution /////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    ////
    //// store small data
    ////

    ////
    //// sign
    ////
    if (this->sample_sign) {
      this->signHist_small[this->currBinCount0] = this->sign * fermionicSign;
    }


    ////
    //// number of particles per component
    //// (must occur prior to potential energy sampling due to the dependence on "N")
    ////
    if (this->sample_numParticles) {
      for (unsigned a = 0; a < numComps; a++) {
        this->numParticlesHist[this->currBinCount0 * numComps + a] = this->numParticles[a] / (double) tMax;
      }
    }


    ////
    //// kinetic energy
    ////
    if (this->sample_kinetEnergy) {
      for (unsigned a = 0; a < numComps; a++) {
        this->kinetEnergyHist[this->currBinCount0 * numComps + a] = - (double) this->numJumps[a] / this->beta_target;
      }
    }


    ////
    //// exchange energy
    ////
    if (has_J && this->sample_exchaEnergy) {
      for (unsigned ab = 0; ab < numExternalInters; ab++) {
        this->exchaEnergyHist[this->currBinCount0 * numExternalInters + ab] = - (double) this->numExchanges[ab] / this->beta_target;
      }
    }


    ////
    //// potential energy
    ////
    if (this->sample_potenEnergy) {
      if (this->H.isUniform()) {
        this->H.calculatePotenEnergies_hom(this->numParticles, this->potenEnergyHist.begin() + this->currBinCount0 * numComps);
      } else {
        this->H.calculatePotenEnergies_inhom(this->numParticlesAtSite, this->potenEnergyHist.begin() + this->currBinCount0 * numComps);
      }
    }


    ////
    //// on-site interaction energy
    ////
    if (has_U && this->sample_interEnergy) {
      this->H.calculateInterEnergies(this->numParticlesSquared,
                                     this->numParticles,
                                     this->interEnergyHist.begin() + this->currBinCount0 * numInters);
    }


    ////
    //// nearest neighbor interaction energy
    ////
    if (has_U_nn && this->sample_nnInterEnergy) {
      this->nnInterEnergyHist[this->currBinCount0] = this->U_nn;
    }


    ////
    //// total energy
    ////
    if (this->sample_totalEnergy) {
      // kinetic contribution
      for (unsigned a = 0; a < numComps; a++) {
        this->totEnergyHist[this->currBinCount0] += this->kinetEnergyHist[this->currBinCount0 * numComps + a];
      }

      // potential contribution
      for (unsigned a = 0; a < numComps; a++) {
        this->totEnergyHist[this->currBinCount0] += this->potenEnergyHist[this->currBinCount0 * numComps + a];
      }

      // exchange contribution
      if (has_J) {
        for (unsigned ab = 0; ab < numExternalInters; ab++) {
          this->totEnergyHist[this->currBinCount0] += this->exchaEnergyHist[this->currBinCount0 * numExternalInters + ab];
        }
      }

      // on-site interaction contribution
      if (has_U) {
        for (unsigned a = 0; a < numComps; a++) {
          for (unsigned b = a; b < numComps; b++) {
            this->totEnergyHist[this->currBinCount0] += this->interEnergyHist[this->currBinCount0 * numInters + Worm::i_ab(a, b)];
          }
        }
      }

      // nearest neighbor interaction contribution
      if (has_U_nn) {
        this->totEnergyHist[this->currBinCount0] += this->U_nn;
      }
    }


    ////
    //// winding numbers
    ////
    if (this->sample_numWinds) {
      for (unsigned a = 0; a < numComps; a++) {
        for (unsigned d = 0; d < this->numDims; d++) {
          this->numWindsHist[(this->currBinCount0 * numComps + a ) * this->numDims + d] = this->numWinds[a][d];
        }
      }
    }


    ////
    //// store large data
    ////
    if ( ! (this->currBinCount0 % this->largeDataSamplingPeriod)) {
      const unsigned long largeDataStoreCount = this->currBinCount0 / this->largeDataSamplingPeriod;

      if (this->sample_sign) {
        this->signHist_large[largeDataStoreCount] = this->sign * fermionicSign;
      }

      ////
      //// particle number per site
      ////
      if (this->sample_numParticlesAtSite) {
        transform(this->numParticlesAtSite.begin(),
                  this->numParticlesAtSite.end(),
                  this->numParticlesAtSiteHist.begin() + largeDataStoreCount * this->numSites * numComps,
                  [&] (double n) { return n / tMax; });
      }


      ////
      //// particle flow
      ////
      if (this->sample_flow) {
        copy(this->flow.begin(),
             this->flow.end(),
             this->flowHist.begin() + largeDataStoreCount * this->numSites * numComps * this->numDims);
      }


      ////
      //// store the instant particle number
      ////
      if (this->sample_instantParticleNum) {
        const auto beg = this->instantParticleNumHist.begin() + largeDataStoreCount * this->numSites * numComps;
        for (unsigned i = 0; i < this->numSites; i++) {
          copy(this->sites[i][0]->pop.begin(),
               this->sites[i][0]->pop.end(),
               beg + i * numComps);
        }
      }
    }


    ////
    //// check whether or not we have filled up the small containers and they need to be reallocated
    ////
    if (this->sample_hist) {
      if (this->currBinCount0 == this->signHist_small.size()) {
        // temporary store current size
        const auto prevSize = this->signHist_small.size();

        ////
        //// reallocate the vectors
        ////
        this->signHist_small.resize(  2 * this->signHist_small.size());
        this->numParticlesHist.resize(2 * this->numParticlesHist.size());
        this->kinetEnergyHist.resize( 2 * this->kinetEnergyHist.size());
        this->potenEnergyHist.resize( 2 * this->potenEnergyHist.size());
        this->totEnergyHist.resize(   2 * this->totEnergyHist.size());
        this->numWindsHist.resize(    2 * this->numWindsHist.size());
        this->exchaEnergyHist.resize(  2 * this->exchaEnergyHist.size());
        this->interEnergyHist.resize(  2 * this->interEnergyHist.size());
        this->nnInterEnergyHist.resize(2 * this->nnInterEnergyHist.size());

        cout << "Worm::trySample: increasing size of the small history vectors: " << prevSize << " -> " << this->signHist_small.capacity() << endl;
      }

      ////
      //// check whether or not we have filled up the large containers and they need to be reallocated
      ////
      if ( (this->currBinCount0 / this->largeDataSamplingPeriod)  == this->signHist_large.size()) {
        // temporary store current size
        const auto prevSize = this->signHist_large.size();

        // reallocate the vectors
        this->signHist_large.resize(            2 * this->signHist_large.size());
        this->numParticlesAtSiteHist.resize(    2 * this->numParticlesAtSiteHist.size());
        this->flowHist.resize(                  2 * this->flowHist.size());
        this->instantParticleNumHist.resize(    2 * this->instantParticleNumHist.size());
        this->instantParticleNumConvHist.resize(2 * this->instantParticleNumConvHist.size());

        cout << "Worm::trySample: increasing size of the large history vectors: " << prevSize << " -> " << this->signHist_large.capacity() << endl;
      }
    }
  }

}