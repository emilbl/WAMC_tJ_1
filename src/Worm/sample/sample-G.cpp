#include "../Worm.h"

using namespace std;


void Worm::sampleG () {
  // PROFILING
  // profiler("Worm::sampleG");

  if (allowMultiCompWorm) {
    cout << "Worm::trySample: Implement me!" << endl;
    exit(EXIT_SUCCESS);
  }

  // compute fermionic exchange sign
  int fermionicSign = 1;
  if (this->hasFermionicExchange) {
    fermionicSign = this->computeFermionicExchangeSign();
  }

  const auto actComp = this->actComps[0];

  ////
  //// determine temporal index
  ////
  unsigned temporal_i;
  int fermionicTimeOrderingSign = 1;
  if (this->isCanonical[actComp]) {
    // [-Δ, Δ] -> [0, 2Δ]
    const long shifted = this->numParticles[actComp] + this->_maxNsDiff[actComp] - this->Ns[actComp] * tMax;
    // [0, 2Δ] -> [0, N - 1]
    // (2Δ + 1 -> N)
    temporal_i = (shifted * G_hist_N) / (2 * this->_maxNsDiff[actComp] + 1);

    if (temporal_i >= G_hist_N) {
      cout << "Worm::trySample: canonical ERROR temporal_i >= G_hist_N" << endl;
      exit(EXIT_SUCCESS);
    }

    ////
    //// time ordering
    ////
    if ( ! this->isBosonic[actComp]) {
      // G = <T[a a^†]> = ± <T[a^† a]>
      // negative times should have a sign flip
      if (this->numParticles[actComp] < this->Ns[actComp] * tMax) fermionicTimeOrderingSign = -1;
    }

  } else {
    ////
    //// Using the formula
    ////
    ////   G(t_rem - t_ins < 0) = ± G(t_rem + β - t_ins > 0)
    ////
    //// we may only consider the time-ordered part
    ////

    // figure out if the worm carries a particle or a hole
    const bool particleWorm = this->headSegments[0]->pop[actComp] > this->headSegments[0]->end->outg->pop[actComp];

    // dt = [0, beta]
    const auto dt = particleWorm > 0 ?
                    positiveModulo<long, long>(this->headSegments[0]->end->t - this->tailSegments[0]->beg->t, tMax) :
                    positiveModulo<long, long>(this->tailSegments[0]->beg->t - this->headSegments[0]->end->t, tMax);

    temporal_i = dt * G_hist_N / tMax;

    if (temporal_i >= G_hist_N) {
      cout << "Worm::trySample: grand canonical ERROR temporal_i >= G_hist_N" << endl;
      exit(EXIT_SUCCESS);
    }

    ////
    //// time ordering
    ////
    if ( ! this->isBosonic[actComp]) {
      // G = <T[a a^†]> = ± <T[a^† a]>
      cout << "Worm::trySample: implement fermionic time ordering sign in the case of grand canonical ensemble" << endl;
      exit(EXIT_SUCCESS);
    }
  }

  ////
  //// determine spatial index
  ////
  const unsigned spatial_i = this->spatialIndexDifference(this->headSegments[0]->siteIndex,
                                                          this->tailSegments[0]->siteIndex);

  if (temporal_i < 0 || temporal_i >= G_hist_N) {
    cout << "Worm::trySample: temporal_i < 0 || temporal_i >= G_hist_N" << endl;
    exit(EXIT_SUCCESS);
  } else if (spatial_i < 0 || spatial_i >= this->numSites){
    cout << "Worm::trySample: spatial_i < 0 || spatial_i >= this->spatial_i" << endl;
    exit(EXIT_SUCCESS);
  }

  ////
  //// bin
  ////
  this->G_hist[actComp][spatial_i       * G_hist_N + temporal_i] += this->sign * fermionicSign * fermionicTimeOrderingSign;
  this->G_hist_count[actComp][spatial_i * G_hist_N + temporal_i] += 1;
}