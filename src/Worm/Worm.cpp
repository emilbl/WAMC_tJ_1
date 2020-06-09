#include "Worm.h"

using namespace std;



static bool _hasFermionicExchange (
  const array<bool, numComps>     isBosonic,
  const array<bool, numComps>     isCanonical,
  const array<unsigned, numComps> Ns
) {
  for (unsigned a = 0; a < numComps; a++) {
    if ( ! isBosonic[a] && ( ! isCanonical[a] || Ns[a] > 1 )) {
      return true;
    }
  }

  return false;
}



Worm::Worm (
  Hamiltonian &                   H,
  const Lattice &                 lattice,
  const Model &                   model,
  Pseudorandom &                  pseudoRandom,
  const WriteToFile &             writeToFile,
  const nlohmann::json &          sampleOptions,
  const nlohmann::json &          computeOptions,
  const nlohmann::json &          saveOptions,
  const array<double, numComps> & maxNsDiff,
  const array<bool, numComps> &   muOptimization,
  const array<double, numComps> & avgNsDiff,
  const bool                      devMode
) :
  devMode{devMode},
  H{H},
  lattice{lattice},
  pseudoRandom{pseudoRandom},
  writeToFile{writeToFile},
  beta_target{model.getBeta()},
  isBosonic{model.getIsBosonic()},
  isCanonical{model.getIsCanonical()},
  Ns{model.getNs()},
  maxNsDiff{maxNsDiff},
  muOptimization{muOptimization},
  targetAvgNsDiff{avgNsDiff},
  hasFermionicExchange{_hasFermionicExchange(isBosonic, isCanonical, Ns)},
  maxNumWorms{model.getMaxNumWorms()},

  //
  largeDataSamplingPeriod{      sampleOptions["largeDataSamplingPeriod"].get<unsigned>()},

  //
  sample_sign{                    sampleOptions["quantities"]["sign"].get<bool>()},
  sample_numParticles{            sampleOptions["quantities"]["numParticles"].get<bool>()},
  sample_numWinds{                sampleOptions["quantities"]["numWinds"].get<bool>()},
  sample_kinetEnergy{             sampleOptions["quantities"]["kinetEnergy"].get<bool>()},
  sample_exchaEnergy{             sampleOptions["quantities"]["exchaEnergy"].get<bool>()},
  sample_potenEnergy{             sampleOptions["quantities"]["potenEnergy"].get<bool>()},
  sample_interEnergy{             sampleOptions["quantities"]["interEnergy"].get<bool>()},
  sample_totEnergy{               sampleOptions["quantities"]["totEnergy"].get<bool>()},
  sample_nnInterEnergy{           sampleOptions["quantities"]["nnInterEnergy"].get<bool>()},
  sample_numParticlesAtSite{      sampleOptions["quantities"]["numParticlesAtSite"].get<bool>()},
  sample_flow{                    sampleOptions["quantities"]["flow"].get<bool>()},
  sample_instantParticleNum{      sampleOptions["quantities"]["instantParticleNum"].get<bool>()},
  sample_instantParticleNumConv{  sampleOptions["quantities"]["instantParticleNumConv"].get<bool>()},

  //
  sample_GreensFunction{sampleOptions["quantities"]["GreensFunction"].get<bool>()},

  //
  sample_fourPointCorrelator{sampleOptions["quantities"]["fourPointCorrelator"].get<bool>()},
  sample_densityDensityCorrelator{sampleOptions["quantities"]["densityDensityCorrelator"].get<bool>()},

  //
  compute_numParticles{   computeOptions["quantities"]["numParticles"]},
  compute_numWinds    {   computeOptions["quantities"]["numWinds"]},
  compute_kinetEnergy{    computeOptions["quantities"]["kinetEnergy"]},
  compute_exchaEnergy{    computeOptions["quantities"]["exchaEnergy"]},
  compute_potenEnergy{    computeOptions["quantities"]["potenEnergy"]},
  compute_interEnergy{    computeOptions["quantities"]["interEnergy"]},
  compute_nnInterEnergy{  computeOptions["quantities"]["nnInterEnergy"]},
  compute_osParticleProd{ computeOptions["quantities"]["osParticleProd"].get<bool>()},
  compute_nnParticleProd{ computeOptions["quantities"]["nnParticleProd"].get<bool>()},
  compute_nnnParticleProd{computeOptions["quantities"]["nnnParticleProd"].get<bool>()},
  compute_C1{             computeOptions["quantities"]["C1"]},
  compute_C2{             computeOptions["quantities"]["C2"]},
  compute_C3{             computeOptions["quantities"]["C3"]},
  compute_particleCorr{   computeOptions["quantities"]["particleCorr"]},
  compute_particleConv{   computeOptions["quantities"]["particleConv"]},
  compute_average{        computeOptions["average"]},

  //
  saveWarmUpData{       saveOptions["saveWarmUpData"].get<bool>()},
  savePresimulationData{saveOptions["savePresimulationData"].get<bool>()}
{
  ////
  //// set initial beta
  ////
  this->setBeta(this->beta_target);

  ////
  //// set up data containers
  ////
  fill(this->numWinds.begin(), this->numWinds.end(), vector<int>(this->numDims));

  ////
  //// set up storage containers
  ////

  // small data history vectors
  if (this->sample_sign)          this->signHist_small.resize(reservedHistSize);
  if (this->sample_numParticles)  this->numParticlesHist.resize( reservedHistSize * numComps);
  if (this->sample_numWinds)      this->numWindsHist.resize(     reservedHistSize * numComps * this->numDims);
  if (this->sample_kinetEnergy)   this->kinetEnergyHist.resize(  reservedHistSize * numComps);
  if (this->sample_exchaEnergy)   this->exchaEnergyHist.resize(  reservedHistSize * numExternalInters);
  if (this->sample_potenEnergy)   this->potenEnergyHist.resize(  reservedHistSize * numComps);
  if (this->sample_interEnergy)   this->interEnergyHist.resize(  reservedHistSize * numInters);
  if (this->sample_totEnergy)     this->totEnergyHist.resize(    reservedHistSize);
  if (this->sample_nnInterEnergy) this->nnInterEnergyHist.resize(reservedHistSize);

  // large data history vectors
  if (this->sample_sign)                   this->signHist_large.resize(            reservedHistSize / this->largeDataSamplingPeriod);
  if (this->sample_numParticlesAtSite)     this->numParticlesAtSiteHist.resize(    reservedHistSize / this->largeDataSamplingPeriod * numSites * numComps);
  if (this->sample_flow)                   this->flowHist.resize(                  reservedHistSize / this->largeDataSamplingPeriod * numSites * numComps * this->numDims);
  if (this->sample_instantParticleNum)     this->instantParticleNumHist.resize(    reservedHistSize / this->largeDataSamplingPeriod * numSites * numComps);
  if (this->sample_instantParticleNumConv) this->instantParticleNumConvHist.resize(reservedHistSize / this->largeDataSamplingPeriod * numSites * numInters);


  // find out possible site differences depending on periodicities
  this->numTransl = 1;
  for (unsigned d = 0; d < this->numDims; d++) {
    if (this->isPeriodic[d]) this->numTransl *= this->latticeSize[d];
  }


  ////
  //// setup C1, C2 and C3 containers
  ////
  if (modelType == tJ) {
    if (this->Ns[1] == 1) {
      // save signs separate
      if (this->compute_C1) this->avgC1.resize(4 * this->numSites * this->halfNumNNs, 0);
      if (this->compute_C2) this->avgC2.resize(4 * this->numSites * this->halfNumNNs, 0);
      if (this->compute_C3) this->avgC3.resize(4 * this->numSites * this->halfNumNNs, 0);
    } else if (this->Ns[1] == 2) {
      if (this->numDims == 2) {
        const auto numConfigs = this->Cs[0] * this->latticeSize[1] + this->Cs[1];
        if (this->compute_C1) this->avgC1.resize(numConfigs * this->numSites * this->halfNumNNs, 0);
        if (this->compute_C2) this->avgC2.resize(numConfigs * this->numSites * this->halfNumNNs, 0);
        if (this->compute_C3) this->avgC3.resize(numConfigs * this->numSites * this->halfNumNNs, 0);
      } else {
        cout << "Worm::Worm:: implement me for numDims == 1" << endl;
        exit(EXIT_SUCCESS);
      }
    }
  }


  ////
  //// setup Green's function histogram
  ////
  if (this->sample_GreensFunction) {
    if (this->isTranslInvariant) {
      // the lattice must be periodic in any direction since the Green's
      // function is assumed to be a function of the differences in positions

      for (unsigned a = 0; a < numComps; a++) {
        this->G_hist[a].resize(      this->numTransl * G_hist_N, 0);
        this->G_hist_count[a].resize(this->numTransl * G_hist_N, 0);
      }
    } else {
      cout << settings::cout::enterRed
           << "Worm::Worm WARNING: cannot sample Green's functions unless the lattice is translation invariant" << endl
           << settings::cout::resetStyle;
    }
  }


  ////
  //// setup four point correlator histogram
  ////
  if (this->sample_fourPointCorrelator) {
    for (unsigned a = 0; a < numComps; a++) {
      // must live in the canonical ensemble and have zero particles
      if (this->isCanonical[a] && this->Ns[a] == 0) {
        this->fpc_hist[a].resize(      fpc_hist_N * fpc_hist_N * fpc_hist_N, 0);
        this->fpc_hist_count[a].resize(fpc_hist_N * fpc_hist_N * fpc_hist_N, 0);
      }
    }
  }


  ////
  //// setup density-density correlator histogram
  ////
  if (this->sample_densityDensityCorrelator) {
    for (unsigned a = 0; a < numComps; a++) {
      // must live in the canonical ensemble and have zero particles
      if (this->isCanonical[a] && this->Ns[a] == 0) {
        this->ddc_hist[a].resize(      ddc_hist_N_o * ddc_hist_N_d * this->numTransl, 0);
        this->ddc_hist_count[a].resize(ddc_hist_N_o * ddc_hist_N_d * this->numTransl, 0);
      }
    }
  }



  if (this->compute_nnParticleProd && ! this->isTranslInvariant) {
    cout << "Worm::sampleZ: ERROR: implement \"nnParticleProd\" for translational invariant systems" << endl;
    exit(EXIT_SUCCESS);
  }
  if (this->compute_nnParticleProd && ! (this->lattice.name != "1d-chain" || this->lattice.name != "2d-rectangular")) {
    cout << "Worm::sampleZ: ERROR: implement \"nnParticleProd\" for translational invariant systems" << endl;
    exit(EXIT_SUCCESS);
  }


  ////
  //// prepare a vector with G-sector update procedure functions
  ////
  #pragma message("prepare a vector with G-sector update procedure functions")



  // weights for choosing the corresponding update procedure (to be read from file)
  this->Wremove        = 1;
  this->Wjump          = 1;
  this->WantiJump      = 1;
  this->Wreconnect     = 1;
  this->WantiReconnect = 1;
  this->WtimeShift     = 1;
  this->Wovertake      = modelType == tJ ? 0 : 1;
  this->Wfollow        = has_J ? 1 : 0;

  // the total weight of the update procedures in the G sector
  this->Wtot = this->Wremove
             + this->Wjump
             + this->WantiJump
             + this->Wreconnect
             + this->WantiReconnect
             + this->WtimeShift
             + this->Wovertake
             + this->Wfollow;

  // populate
  this->weightedUpdateProcedures.reserve(this->Wtot);
  this->weightedUpdateProcedures.insert(weightedUpdateProcedures.begin(), this->Wremove,        &Worm::remove);
  this->weightedUpdateProcedures.insert(weightedUpdateProcedures.begin(), this->Wjump,          &Worm::jump);
  this->weightedUpdateProcedures.insert(weightedUpdateProcedures.begin(), this->WantiJump,      &Worm::antiJump);
  this->weightedUpdateProcedures.insert(weightedUpdateProcedures.begin(), this->Wreconnect,     &Worm::reconnect);
  this->weightedUpdateProcedures.insert(weightedUpdateProcedures.begin(), this->WantiReconnect, &Worm::antiReconnect);
  this->weightedUpdateProcedures.insert(weightedUpdateProcedures.begin(), this->WtimeShift,     &Worm::timeShift);
  this->weightedUpdateProcedures.insert(weightedUpdateProcedures.begin(), this->Wovertake,      &Worm::overtake);
  this->weightedUpdateProcedures.insert(weightedUpdateProcedures.begin(), this->Wfollow,        &Worm::follow);


  ////
  ////
  ////



  ////
  //// create reference to t=0 segments
  ////
  this->sites.resize(this->numSites);
  for (unsigned i = 0; i < this->numSites; i++) {
    // allocate sufficient memory
    this->sites[i].reserve(initialSegmentContainerSize);

    // initial empty site
    array<unsigned, numComps> pop = {};
    auto segment = make_shared<Segment>(i,
                                        pop,
                                        nullptr,
                                        nullptr);

    this->sites[i].push_back(segment);
  }

  ////
  //// in case of canonical component, fill segments according to statistics
  ////
  for (unsigned a = 0; a < numComps; a++) {
    if (this->isCanonical[a]) {

      // make sure Pauli exclusion principle is fulfilled
      if ( ! this->isBosonic[a] && this->Ns[a] > this->numSites) {
        cout << settings::cout::enterRed
             << "Worm::Worm ERROR: unable to fulfill Pauli exclusion principle  ->  EXIT" << endl
             << settings::cout::resetStyle;
        this->shutDown();
      }

      unsigned N = 0;
      while (N < this->Ns[a]) {
        // pick a site on random
        const unsigned i = this->pseudoRandom.template Uint<unsigned>(0, this->numSites - 1);

        // in case of fermions, if the site is occupied already, continue looking elsewhere
        if ( ! this->isBosonic[a] && this->sites[i][0]->pop[a]) {
          continue;
        }

        // update the particle number
        this->numParticlesAtSite[i * numComps + a] += tMax;

        // update the particle number
        this->numParticles[a] += tMax;

        // update the particle number square
        for (unsigned b = 0; b < numComps; b++) {
          const unsigned i_ab = Worm::i_ab(a, b),
                         n_a = this->sites[i][0]->pop[a],
                         n_b = this->sites[i][0]->pop[b];

          if (a == b) {
            this->numParticlesSquared[i_ab] += (  pow(n_a + 1, 2)
                                                - pow(n_a, 2) ) * tMax;
          } else {
            this->numParticlesSquared[i_ab] += (  (n_a + 1) * n_b
                                                - n_a * n_b ) * tMax;
          }
        }

        // add the particle to the current configuration
        this->sites[i][0]->pop[a]++;

        // increase counter
        N++;
      }
    }
  }

  ////
  //// reset data containers
  ////
  this->resetDataContainers();

  ////
  //// always check worm at initialization
  ////
  this->checkWormConfiguration();
  this->checkWormQuantities();
}


const Lattice & Worm::getLattice () const {
  return this->lattice;
}


const vector<shared_ptr<Segment> > & Worm::getSegments (unsigned i) const {
  return this->sites[i];
}


double Worm::getBeta () const {
  return this->beta;
}

double Worm::getBeta_target () const {
  return this->beta_target;
}

void Worm::setBeta (
  const double & beta
) {
  this->beta = beta;

  // recompute max particle difference for the canonical components
  for (unsigned a = 0; a < numComps; a++) {
    if (this->isCanonical[a]) {
      // need to be "beta_target" rather than "beta" since "_maxNsDiff" should be independent of current "beta"
      this->_maxNsDiff[a] = min(tMax, (long) (this->maxNsDiff[a] / this->beta_target * tMax) ) - 1;
    }
  }
}

void Worm::setMaxNsDiff (
  const array<double, numComps> & maxNsDiff
) {
  // recompute max particle difference for the canonical components
  for (unsigned a = 0; a < numComps; a++) {
    if (this->isCanonical[a]) {
      cout << "Worm::setMaxNsDiff: uncertain about if it does what it is supposed to do..." << endl;
      this->maxNsDiff[a]  = min(this->beta, maxNsDiff[a]);
      this->_maxNsDiff[a] = min(tMax, (long) (this->maxNsDiff[a] / this->beta * tMax)) - 1;
    }
  }
}




// long Worm::gettMax () const {
//   return tMax;
// }


unsigned long long Worm::getCurrUpdateCount () const {
  return this->currUpdateCount;
}


array<int, numComps> Worm::getWormPop (
  const unsigned wormIndex
) const {
  if (allowMultiCompWorm) {
    cout << "Worm::getWormPop: implement me" << endl;
    return this->headSegments[wormIndex] ? this->wormPop : array<int, numComps>{};
  } else {
    array<int, numComps> wp = {};
    wp[this->actComps[wormIndex]] = this->actPops[wormIndex];
    return wp;
  }
}

unsigned Worm::getNumWorms () const {
  return this->numWorms;
}


void Worm::timeModulo (long & time) const {
  time = (time % tMax + tMax) % tMax;
}


/* OLD ? */  bool Worm::isHeadSegment (shared_ptr<Segment> segment) const {
/* OLD ? */    return find(this->headSegments.begin(), this->headSegments.end(), segment) != this->headSegments.end();
/* OLD ? */  }
/* OLD ? */
/* OLD ? */
/* OLD ? */  bool Worm::isTailSegment (shared_ptr<Segment> segment) const {
/* OLD ? */    return find(this->tailSegments.begin(), this->tailSegments.end(), segment) != this->tailSegments.end();
/* OLD ? */  }
/* OLD ? */
/* OLD ? */
/* OLD ? */  bool Worm::isActiveComponent (unsigned comp) const {
/* OLD ? */    if (allowMultiCompWorm) {
/* OLD ? */      return !! this->wormPop[comp];
/* OLD ? */    } else {
/* OLD ? */      for (const auto actComp : this->actComps) {
/* OLD ? */        if (actComp == comp) return true;
/* OLD ? */      }
/* OLD ? */    }
/* OLD ? */
/* OLD ? */    return false;
/* OLD ? */  }

shared_ptr<Segment> Worm::getHeadSegment (
  const unsigned w
) const {
  return this->headSegments[w];
}
shared_ptr<Segment> Worm::getTailSegment (
  const unsigned w
) const {
  return this->tailSegments[w];
}




template<typename To, typename Ti1, typename Ti2>
std::array<To, numComps> Worm::diff (
  const std::array<Ti1, numComps> & pop1,
  const std::array<Ti2, numComps> & pop2
) {
  array<To, numComps> diff = {};

  for (unsigned i = 0; i < numComps; i++) {
    diff[i] = (int) pop1[i] - (int) pop2[i];
  }

  return diff;
}
template array<unsigned, numComps> Worm::diff<unsigned, unsigned, int> (const array<unsigned, numComps> &, const array<int, numComps> &);
template array<int, numComps> Worm::diff<int, unsigned, unsigned> (const array<unsigned, numComps> &, const array<unsigned, numComps> &);
template array<int, numComps> Worm::diff<int, unsigned, int> (const array<unsigned, numComps> &, const array<int, numComps> &);
template array<int, numComps> Worm::diff<int, int, int> (const array<int, numComps> &, const array<int, numComps> &);

template<typename To, typename Ti1, typename Ti2>
std::array<To, numComps> Worm::sum (
  const std::array<Ti1, numComps> & pop1,
  const std::array<Ti2, numComps> & pop2
) {
  array<To, numComps> sum = {};

  for (unsigned i = 0; i < numComps; i++) {
    sum[i] = (int) pop1[i] + (int) pop2[i];
  }

  return sum;
}
template array<unsigned, numComps> Worm::sum<unsigned, unsigned, unsigned> (const array<unsigned, numComps> &, const array<unsigned, numComps> &);
template array<unsigned, numComps> Worm::sum<unsigned, unsigned, int> (const array<unsigned, numComps> &, const array<int, numComps> &);
template array<int, numComps> Worm::sum<int, int, int> (const array<int, numComps> &, const array<int, numComps> &);


void Worm::add (
  array<unsigned, numComps> & pop1,
  const array<int, numComps> & pop2
) {
  for (unsigned i = 0; i < numComps; i++) {
    pop1[i] += (unsigned) pop2[i];
  }
}

void Worm::add (
  vector<int> & v1,
  const vector<int> & v2
) {
  for (unsigned d = 0; d < v1.size(); d++) {
    v1[d] += v2[d];
  }
}

void Worm::subtract (
  array<unsigned, numComps> & pop1,
  const array<int, numComps> & pop2
) {
  for (unsigned i = 0; i < numComps; i++) {
    pop1[i] -= (unsigned) pop2[i];
  }
}

void Worm::subtract (
  vector<int> & v1,
  const vector<int> & v2
) {
  for (unsigned d = 0; d < v1.size(); d++) {
    v1[d] -= v2[d];
  }
}


template<typename T1, typename T2>
void Worm::addOrSubtract (
  array<T1, numComps> & arr1,
  const array<T2, numComps> & arr2,
  const int sign
) {
  for (unsigned a = 0; a < numComps; a++) {
    arr1[a] += sign * arr2[a];
  }
}
template void Worm::addOrSubtract<unsigned, int> (
  array<unsigned, numComps> &,
  const array<int, numComps> &,
  const int
);


template<typename T1, typename T2>
void Worm::addOrSubtract (
  vector<T1> & arr1,
  const vector<T2> & arr2,
  const int sign
) {
  for (unsigned i = 0; i < arr1.size(); i++) {
    arr1[i] += sign * arr2[i];
  }
}
template void Worm::addOrSubtract<int, int> (
  vector<int> &,
  const vector<int> &,
  const int
);


bool Worm::mismatchingPopDiff (
  const array<unsigned, numComps> & pop2,
  const array<unsigned, numComps> & pop1,
  const unsigned                    actComp,
  const int                         actPop
) {
  for (unsigned a = 0; a < numComps; a++) {
    const int dPop = (int) pop2[a] - (int) pop1[a];

    if (
      (a == actComp && dPop != actPop) ||
      (a != actComp && dPop != 0)
    ) return true;
  }

  // the difference matched
  return false;
}
bool Worm::mismatchingPopDiff (
  const array<unsigned, numComps> & pop2,
  const array<unsigned, numComps> & pop1,
  const array<int, numComps> &      wormPop
) {
  for (unsigned a = 0; a < numComps; a++) {
    const int dPop = (int) pop2[a] - (int) pop1[a];

    if (dPop != wormPop[a]) return true;
  }

  // the difference matched
  return false;
}


// Fermi what!? Rather Pauli exclusion principle?
bool Worm::wouldBreakFermiPrinc (
  array<unsigned, numComps> & segmentPop,
  array<int, numComps> & wormPop
) {
  for (unsigned a = 0; a < numComps; a++) {
    if ( (! this->isBosonic[a]) && ( (int) (segmentPop[a] * wormPop[a]) != 0 )) {
      return true;
    }
  }

  return false;
}


unsigned Worm::i_ab (
  const unsigned a,
  const unsigned b
) {
  // numInter

  // event though we are using the lower triangular part of the component matrix
  // the stored values of the number of particles squared should only be in the upper half
  const unsigned A = min(a, b),
                 B = max(a, b);

  return 0.5 * A * (2 * numComps - 1 - A) + B;
}

unsigned Worm::i_ab_external (
  const unsigned a,
  const unsigned b
) {
  // numExternalInter

  // like i_ab but without the diagonal
  const unsigned A = min(a, b),
                 B = max(a, b);

  return 0.5 * A * (2 * numComps - 3 - A) + B - 1;
}


double Worm::int2time (
  const double & beta
) const {
  return beta / (double) tMax;
}




void Worm::shutDown () const {

  // list of update function names
  vector<string> updates = {
    "insert",
    "remove",
    "jump",
    "anti-jump",
    "reconnect",
    "anti-reconnect",
    "time-shift",
    "overtake",
    "follow",
    "swap"
  };

  // worm population
  vector<array<int, numComps> > worms = {};
  for (unsigned w = 0; w < this->numWorms; w++) {
    if (allowMultiCompWorm) {
      worms.push_back(this->wormPop);
    } else {
      array<int, numComps> pop = {};
      pop[this->actComps[w]] = this->actPops[w];
      worms.push_back(pop);
    }
  }

  // output some data
  cout << settings::cout::enterRed
       << "Worm::shutDown: shutting down after " << this->currUpdateCount << " updates and " << this->currBinCount0 << " / " << this->currBinCount1 << " / " << this->currBinCount2
       << " binnings. The most recent/ongoing update being \"" << updates[this->currentUpdate] << "\""
       << " with worms " << worms << "."
       << settings::cout::resetStyle << endl;

  // shut down
  exit(EXIT_FAILURE);
}
