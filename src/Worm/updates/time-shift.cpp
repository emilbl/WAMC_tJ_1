#include "../Worm.h"

using namespace std;

void Worm::timeShift (
  const double & beta,
  const unsigned wormIndex,
  const double & P_insert,
  const bool     test
) {
  // PROFILING
  // profiler("Worm::timeShift");

  // VERBOSE: let the terminal know that a jump has begun
  if (verbose && test) {
    cout << settings::cout::enterYellow << "testing to shift time back (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  } else if (verbose) {
    cout << settings::cout::enterYellow << "about to shift time (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  }

  // store the current update procedure
  constexpr unsigned updateNum = 6;
  this->currentUpdate = updateNum;

  // update statistics
  if ( ! test) this->updateStatistics[3 * updateNum]++;

  ////
  //// DEBUG: worm weight check
  ////
  if (debug_major && debugFrom <= this->currUpdateCount) this->preCheckWormWeight(beta);

  ////
  //// store some variables
  ////
  const auto headSegment = this->headSegments[wormIndex];
  const auto tailSegment = this->tailSegments[wormIndex];
  const unsigned i = headSegment->siteIndex;
  const auto actComp = allowMultiCompWorm ? 0 : this->actComps[wormIndex];
  const auto actPop = allowMultiCompWorm ? 0 : this->actPops[wormIndex];


  ////
  //// interval bounds for the new head time
  ////
  long tHeadMin_gcan = headSegment->beg ?
                       headSegment->beg->t + 1 :                   // not crossing
                       this->sites[i].back()->beg->t + 1 - tMax,   // crossing t=0
       tHeadMax_gcan = headSegment->end->outg->end ?
                       headSegment->end->outg->end->t - 1 :    // not crossing
                       this->sites[i][0]->end->t - 1 + tMax;   // crossing t=beta

  // take into account canonical limitations
  if (allowMultiCompWorm) {
    cout << "Worm::insert: Implement me" << endl;
    // for (unsigned a = 0; a < numComps; a++) {
    //   if (this->isCanonical[a]) {
    //     // whether or not we are beyond the fixed particle number
    //     const auto beyond = (long) this->numParticles[a] - (long) tMax * this->Ns[a];

    //     if (
    //       (beyond > 0 && actPop > 0) ||
    //       (beyond < 0 && actPop < 0)
    //     ) {
    //       const long tHeadMax_can = tailSegment->beg->t < headSegment->end->t ?
    //                                 tailSegment->beg->t + tMax - 1 :
    //                                 tailSegment->beg->t - 1;

    //       tHeadMax_gcan = min(tHeadMax_gcan, tHeadMax_can);
    //     } else if (
    //       (beyond > 0 && actPop < 0) ||
    //       (beyond < 0 && actPop > 0)
    //     ) {
    //       const long tHeadMin_can = tailSegment->beg->t < headSegment->end->t ?
    //                                 tailSegment->beg->t + 1 :
    //                                 tailSegment->beg->t - tMax + 1;

    //       tHeadMin_gcan = max(tHeadMin_gcan, tHeadMin_can);
    //     }
    //   }
    // }
  } else if (this->isCanonical[actComp]) {
    // whether or not we are beyond the fixed particle number
    const auto beyond = (long) this->numParticles[actComp] - (long) tMax * this->Ns[actComp];

    // upper limit
    const long tHeadMax_can = headSegment->end->t + this->_maxNsDiff[actComp] - actPop * beyond;
    tHeadMax_gcan = min(tHeadMax_gcan, tHeadMax_can);

    // lower limit
    const long tHeadMin_can = headSegment->end->t - this->_maxNsDiff[actComp] - actPop * beyond;
    tHeadMin_gcan = max(tHeadMin_gcan, tHeadMin_can);
  }

  // convert to constants
  const long tHeadMin = tHeadMin_gcan,
             tHeadMax = tHeadMax_gcan;

  // in order to proceed the interval must have a nonzero length
  if (tHeadMin >= tHeadMax) {
    if (verbose)
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::reconnect: tHeadMin >= tHeadMax  ->  RETURN"
           << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // TEST: make sure that the bounds agree
  if (debug && test) {
    if (tHeadMin != this->tMin_prev) {
      cout << settings::cout::enterRed << "Worm::timeShift: tHeadMin=" << tHeadMin << " vs tMin_prev=" << this->tMin_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tHeadMax != tMax_prev){
      cout << settings::cout::enterRed << "Worm::timeShift: tHeadMax=" << tHeadMax << " vs tMax_prev=" << this->tMax_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// the weight ratio of the proposed local update: W(t') / W(t)
  ////
  double localWeightRatio_base = 1,
         localWeightRatio_exponent = 0;
  if (allowMultiCompWorm) {
    this->H.potenDiff(i,
                      this->wormPop,
                      1,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
    if (has_U) {
      localWeightRatio_exponent += this->H.interDiff(headSegment->end->outg->pop,
                                                     this->wormPop,
                                                     1);
    }
  } else {
    this->H.potenDiff(i,
                      actComp,
                      actPop,
                      1,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
    if (has_U) {
      localWeightRatio_exponent += this->H.interDiff(headSegment->end->outg->pop,
                                                     actComp,
                                                     actPop,
                                                     1);
    }
  }

  ////
  //// parameters for the tHead distribution
  ////
  const long intLength = tHeadMax - tHeadMin;
  const double intervalLength = intLength * this->int2time(beta),
               lambda = localWeightRatio_exponent;


  if (isnan(lambda)) {
    cout << localWeightRatio_exponent << endl;
    cout << settings::cout::enterRed << "Worm::timeShift: lambda=" << lambda << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  // TEST: make sure that the distribution parameters agree
  if (debug && test) {
    if (intLength != this->intLength_prev) {
      cout << settings::cout::enterRed << "Worm::timeShift: intLength=" << intLength << " vs intLength_prev=" << this->intLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (lambda != this->lambda_prev) {
      cout << settings::cout::enterRed << "Worm::timeShift: lambda=" << lambda << " vs lambda_prev=" << this->lambda_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// propose a new head time
  ////
  long tHead;
  if (debug && test) {
    tHead = this->t_prev;
  } else {
    ////
    //// what type of distribution should be used
    ////
    if (expDistrEnabled) {
      const double dt = this->pseudoRandom.template Exp<double>(intervalLength, lambda);
      tHead = tHeadMin + round(dt / this->int2time(beta));
    } else {
      const double dt = this->pseudoRandom.template U<double>(0, intervalLength);
      tHead = tHeadMin + round(dt / this->int2time(beta));
    }
  }

  // in order to proceed the new head time must not occur on the time boundary
  if (tHead == 0 || tHead == tMax) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::timeShift: tHead == 0 || tHead == tMax  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // DEBUG: make sure that tHead is inside the bounds
  if (debug && debugFrom <= this->currUpdateCount) {
    if (tHead < tHeadMin) {
      cout << settings::cout::enterRed << "Worm::timeShift: tHead=" << tHead << " < " << tHeadMin << "=tHeadMin" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tHead > tHeadMax) {
      cout << settings::cout::enterRed << "Worm::timeShift: tHead=" << tHead << " > " << tHeadMax << "=tHeadMax" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// the segment length being modified by this update
  ////
  const long modifiedLength = tHead - headSegment->end->t;

  ////
  //// the periodicity is taken care of inside "update worm"
  ////


  ////
  //// other weight contributions
  ////
  double exp_U_nn  = 0;

  // nearest neighbor interaction
  if (has_U_nn) {
    if (allowMultiCompWorm) {
      cout << "Worm::timeShift: TODO: implement me!" << endl;
      this->shutDown();
    } else {
      if (modifiedLength > 0) {
        exp_U_nn += computeNNinteractionDiff(i,
                                             headSegment->pop,
                                             headSegment->end->outg->pop,
                                             headSegment->end->t,
                                             tHead);
      } else {
        exp_U_nn += computeNNinteractionDiff(i,
                                             headSegment->end->outg->pop,
                                             headSegment->pop,
                                             tHead,
                                             headSegment->end->t);
      }
    }
  }


  ////
  //// The acceptance ratio
  //// OBS: cancellation of exponential dependence
  ////
  ////                  ╭      W(t')  ╮
  //// A(t -> t') = min │ 1, -------- │   since g(t) = g(t') (symmetric update)
  ////                  ╰      W(t)   ╯
  ////
  const double R = abs(localWeightRatio_base)
                 * exp( - (expDistrEnabled ? 0 : localWeightRatio_exponent * this->int2time(beta) * modifiedLength)
                        - exp_U_nn * this->int2time(beta) ),
               A = min(1., R);

  // DEBUG: the acceptance ration must always be larger than zero
  if (debug && debugFrom <= this->currUpdateCount && A < 0) {
    cout << settings::cout::enterRed << "Worm::timeShift: invalid value of the acceptance ration A=" << A << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  ////
  //// TEST: some final tests
  ////
  if (debug && test) {
    // modified lengths must agree
    if (modifiedLength != -this->modifiedLength_prev) {
      cout << settings::cout::enterRed << "Worm::timeShift: modifiedLength=" << modifiedLength << " vs -modifiedLength_prev=" << -this->modifiedLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }

    // the acceptance ratios must agree
    if (abs(R - 1 / this->R_prev) > diffTresh) {
      cout << settings::cout::enterRed << "Worm::timeShift: R=" << R << " vs 1/R_prev=" << 1 / this->R_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }

    // everything is OK -> proceed with new update procedure
    if (verbose) cout << settings::cout::enterGreen << "...test successful" << settings::cout::resetStyle << endl;
    return;
  }

  // update statistics
  if ( ! test) this->updateStatistics[3 * updateNum + 1]++;

  ////
  //// accept or reject
  ////
  if (A == 1 || A > this->pseudoRandom.template U<double>(0, 1)) {
    // accept proposed update
    if (verbose) cout << settings::cout::enterGreen << "...proposed update accepted (A=" << A << ")..." << settings::cout::resetStyle << endl;

  } else {
    // reject proposed update
    if (verbose) cout << settings::cout::enterOrange << "...proposed update rejected (A=" << A << ")..." << settings::cout::resetStyle << endl;

    return;
  }

  // update statistics
  if ( ! test) this->updateStatistics[3 * updateNum + 2]++;

  ////
  //// update quantities
  ////
  this->sign *= sgn(localWeightRatio_base);
  if (allowMultiCompWorm) {
    for (unsigned a = 0; a < numComps; a++) {
      // update the particle number
      this->numParticlesAtSite[i * numComps + a] += modifiedLength * this->wormPop[a];

      this->numParticles[a] += modifiedLength * this->wormPop[a];

      // update the particle number square
      for (unsigned b = a; b < numComps; b++) {
        this->numParticlesSquared[Worm::i_ab(a, b)] += modifiedLength
                                                     * (  (int) (headSegment->pop[a] * headSegment->pop[b])
                                                        - (int) (headSegment->end->outg->pop[a] * headSegment->end->outg->pop[b]));
      }
    }
  } else {
    // update the particle number
    this->numParticlesAtSite[i * numComps + actComp] += modifiedLength * actPop;

    this->numParticles[actComp] += modifiedLength * actPop;

    // update the particle number square
    for (unsigned b = 0; b < numComps; b++) {
      const unsigned i_ab = Worm::i_ab(actComp, b);

      this->numParticlesSquared[i_ab] += modifiedLength
                                       * (  (int) (  headSegment->pop[actComp]
                                                   * headSegment->pop[b] )
                                          - (int) (  headSegment->end->outg->pop[actComp]
                                                   * headSegment->end->outg->pop[b] ) );
    }
  }

  // update nearest neighbor interaction
  if (has_U_nn) this->U_nn += exp_U_nn / tMax;



  ////
  //// DEBUG: store previous head time before being overwritten
  ////
  if (debug && debugFrom <= this->currUpdateCount) {
    // the previous tHead time taken into account that t=0 might be shifted ± tMax
    this->t_prev = headSegment->end->t - ((tHead > tMax) - (tHead < 0)) * tMax;
  }

  ////
  //// update worm
  ////
  if (tHead < 0) {
    // -1 beta winding

    // split segment previous to the head segment
    const auto newHeadSegment = this->splitSegment(i,
                                                   this->sites[i].size() - 1,
                                                   tHead + tMax);

    // remove former head segment
    this->removeSegment(i, 0);

    // fix population imbalance
    if (allowMultiCompWorm) {
      this->subtract(newHeadSegment->end->outg->pop, this->wormPop);
    } else {
      newHeadSegment->end->outg->pop[actComp] -= actPop;
    }

    // update head
    this->headSegmentIndices[wormIndex] = this->sites[i].size() - 2;
    this->headSegments[wormIndex] = newHeadSegment;
  } else if (tHead > tMax) {
    // +1 beta winding

    // split segment following to the head segment
    const auto newHeadSegment = this->splitSegment(i, 0, tHead - tMax);

    // remove the last segment
    this->removeSegment(i, this->sites[i].size() - 1);

    // fix population imbalance
    if (allowMultiCompWorm) {
      this->add(newHeadSegment->pop, this->wormPop);
    } else {
      newHeadSegment->pop[actComp] += actPop;
    }

    // update head
    this->headSegmentIndices[wormIndex] = 0;
    this->headSegments[wormIndex] = newHeadSegment;
  } else {
    // simply update the head time
    this->headSegments[wormIndex]->end->t = tHead;
  }

  ////
  //// DEBUG: worm weight check
  ////
  if (debug_major && debugFrom <= this->currUpdateCount) {
    this->checkWormWeight(beta,
                          localWeightRatio_base,
                          - localWeightRatio_exponent * this->int2time(beta) * modifiedLength
                          - exp_U_nn * this->int2time(beta) );
  }

  // VERBOSE: let the terminal know this update is complete
  if (verbose) cout << settings::cout::enterGreen << "...time shifted" << settings::cout::resetStyle << endl;

  ////
  //// DEBUG: test an anti-time-shift
  ////
  if (debug && debugFrom <= this->currUpdateCount) {

    // interval bounds for tHead taken into account that t=0 might be shifted ± tMax
    this->tMin_prev = tHeadMin - ((tHead > tMax) - (tHead < 0)) * tMax;
    this->tMax_prev = tHeadMax - ((tHead > tMax) - (tHead < 0)) * tMax;

    // the previous tHead time
    // already stored

    // parameters for the exponential distribution
    this->intLength_prev = intLength;
    this->lambda_prev = lambda;

    // closely related to acceptance ratio
    this->modifiedLength_prev = modifiedLength;
    this->R_prev = R;

    // try going back with the anti-update
    this->timeShift(beta, wormIndex, P_insert, true);
  }
}