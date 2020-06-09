#include "../Worm.h"

using namespace std;

void Worm::remove (
  const double & beta,
  const unsigned wormIndex,
  const double & P_insert,
  const bool     test
) {
  // PROFILING
  // profiler("Worm::remove");

  // VERBOSE: let the terminal know that a removal has begun
  if (verbose && test) {
    cout << settings::cout::enterYellow << "testing to remove back (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  } else if (verbose) {
    cout << settings::cout::enterYellow << "about to remove (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  }

  // store the current update procedure
  constexpr unsigned updateNum = 1;
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
  // const auto tailSegment = this->tailSegments[wormIndex];
  const unsigned i = headSegment->siteIndex;
  const auto actComp = allowMultiCompWorm ? 0 : this->actComps[wormIndex];
  const auto actPop = allowMultiCompWorm ? 0 : this->actPops[wormIndex];
  const auto segAfterHead = headSegment->end->outg;



  ////
  ////
  ////
  unsigned tailIndexAhead  = 0,
           tailIndexBehind = 0;
  shared_ptr<Segment> tailAhead  = nullptr,
                      tailBehind = nullptr;
  if (allowMultiCompWorm) {
    cout << "Worm::remove: implement me" << endl;
  } else {
    for (unsigned wi = 0; wi < this->numWorms; wi++) {
      if (
           this->tailSegments[wi]->siteIndex == i
        && this->actComps[wi] == actComp
      ) {
        // just in front of head
        if (
             (segAfterHead->end ? segAfterHead->end->outg : this->sites[i][1]) == this->tailSegments[wi]
          && headSegment->pop[actComp] == this->tailSegments[wi]->pop[actComp]
        ) {
          tailIndexAhead = wi;
          tailAhead = this->tailSegments[wi];
        }

        // just behind head
        if (
             (headSegment->beg ? headSegment : this->sites[i].back()) == this->tailSegments[wi]
          && this->tailSegments[wi]->beg->inco->pop[actComp] == headSegment->end->outg->pop[actComp]
        ) {
          tailIndexBehind = wi;
          tailBehind = this->tailSegments[wi];
        }
      }
    }
  }




  // ////
  // //// in order to proceed both head and tail must be located on the same lattice site
  // ////
  // if (i != tailSegment->siteIndex) {
  //   if (verbose || (debug && test))
  //     cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
  //          << "Worm::remove: the ends are not on the same site  ->  RETURN" << settings::cout::resetStyle << endl;
  //   if (debug && test && shutItDown) this->shutDown(); else return;
  // }

  // ////
  // //// determine whether or not the worm may be removed by going forward/backward
  // ////
  // // if the worm head may reach the tail moving backward and/or forward
  // const bool forward = segAfterHead->end ?                        // segment after head crossing t=beta?
  //                      segAfterHead->end->outg == tailSegment :   //    n: next segment = tail?
  //                      this->sites[i][1] == tailSegment,          //    y:      -- || --
  //            backward = headSegment->beg ?                      // head segment crossing t=0?
  //                       headSegment == tailSegment :            //    n: head = tail?
  //                       this->sites[i].back() == tailSegment;   //    y:   -- || --

  ////
  //// determine whether or not the worm may be removed by going forward/backward
  ////
  const bool forward  = !! tailAhead,
             backward = !! tailBehind;

  // in order to proceed it must be possible to remove the worm in at least one direction
  if ( ! forward && ! backward) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::remove: ! forward && ! backward  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // // DEBUG: make sure that both forward and backward should be possible
  // //        if and only if the site has three segments (was a flat site)
  // if (debug && debugFrom <= this->currUpdateCount) {
  //   if ((this->sites[i].size() == 3) ^ (forward && backward)) {
  //     cout << settings::cout::enterRed << "Worm::remove: (this->sites[i].size() == 3) ^ (forward && backward)" << settings::cout::resetStyle << endl;
  //     if (shutItDown) this->shutDown();
  //   }
  // }

  // TEST: make sure the worm may be removed in the same way it was created
  if (debug && test) {
    if (this->subtract_prev && ! forward) {
      cout << settings::cout::enterRed << "Worm::remove: this->subtract_prev && ! forward" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if ( ! this->subtract_prev && ! backward) {
      cout << settings::cout::enterRed << "Worm::remove: ! this->subtract_prev && ! backward" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// propose a direction
  ////
  bool goingForward;
  if (debug && test) {
    // insertion by subtraction -> going forward
    // insertion by addition -> going backwards
    goingForward = this->subtract_prev;
  } else {
    goingForward = backward && forward ?
                   !! this->pseudoRandom.template Uint<unsigned>(0, 1) :
                   forward;
  }
  const auto tailIndex = goingForward ? tailIndexAhead : tailIndexBehind;
  const auto tailSegment = goingForward ? tailAhead : tailBehind;

  // inverse probability of proposing the particular direction
  const double W_forwardOrBackward = (forward && backward) ? 2 : 1;

  // VERBOSE: let the terminal know what type of removal is proposed
  if (verbose) {
    array<int, numComps> wp = {};
    if (allowMultiCompWorm) {
      wp = this->wormPop;
    } else {
      wp[actComp] = actPop;
    }

    cout << settings::cout::enterYellow << "...by going "
         << (goingForward ? "forward " : "backwards ") << wp
         << "..." << settings::cout::resetStyle << endl;
  }


  ////
  //// if particle number should be conserved, is such an update allowed
  ////
  if (allowMultiCompWorm) {
    cout << "Worm::remove: implement me" << endl;
    // for (unsigned a = 0; a < numComps; a++) {
    //   if (this->isCanonical[a]) {

    //     const bool beyond = this->numParticles[a] > (unsigned long) tMax * this->Ns[a];

    //     if (
    //       (beyond && this->wormPop[a] > 0 && goingForward) ||
    //       (beyond && this->wormPop[a] < 0 && (! goingForward)) ||
    //       ( (! beyond) && this->wormPop[a] > 0 && (! goingForward) ) ||
    //       ( (! beyond) && this->wormPop[a] < 0 && goingForward )
    //     ) {
    //       if (verbose) {
    //         cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
    //              << "Worm::remove: non-conservation of particle number  ->  RETURN"
    //              << settings::cout::resetStyle << endl;
    //       }
    //       if (debug && test && shutItDown) this->shutDown(); else return;
    //     }
    //   }
    // }
  } else if (this->isCanonical[actComp]) {
    long modifiedLength = (goingForward ? 1 : -1) * (tailSegment->beg->t - headSegment->end->t);

    // take into account tau periodicity
    if (modifiedLength < 0) modifiedLength += tMax;

    // modified number of particles
    const auto dN = modifiedLength * (goingForward ? 1 : -1) * actPop;

    // the overshoot of the particle number from the specified one
    const auto overshoot = abs(((long) this->numParticles[actComp]) + dN - tMax * this->Ns[actComp]);

    // particle number must be inside +- tMax
    if (overshoot > this->_maxNsDiff[actComp]) {
      if (verbose) {
        cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
             << "Worm::remove: particle number out of bound  ->  RETURN"
             << settings::cout::resetStyle << endl;
      }
      if (debug && test && shutItDown) this->shutDown(); else return;
    }
  }


  ////
  //// the population of the modified segment parts after the update
  ////
  const auto segPopAfterUpdate = goingForward ? headSegment->pop : segAfterHead->pop;

  ////
  //// the weight ratio of the proposed local update: W(G) / W(Z)
  ////
  double localWeightRatio_base = 1,
         localWeightRatio_exponent = 0;

  if (allowMultiCompWorm) {
    this->H.potenDiff(i,
                      this->wormPop,
                      goingForward ? -1 : 1,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
    if (has_U) {
      localWeightRatio_exponent += this->H.interDiff(segPopAfterUpdate,
                                                     this->wormPop,
                                                     goingForward ? -1 : 1);
    }
    this->H.discoDiff(segAfterHead->pop,
                      this->wormPop,
                      1,
                      beta,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
  } else {
    this->H.potenDiff(i,
                      actComp,
                      actPop,
                      goingForward ? -1 : 1,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
    if (has_U) {
      localWeightRatio_exponent += this->H.interDiff(segPopAfterUpdate,
                                                     actComp,
                                                     actPop,
                                                     goingForward ? -1 : 1);
    }
    this->H.discoDiff(segAfterHead->pop[actComp],
                      actComp,
                      actPop,
                      1,
                      beta,
                      localWeightRatio_base,
                      localWeightRatio_exponent);


  }

  ////
  //// insert proposal distribution
  ////
  double _W_insert = 1;
  auto proposedActComps = this->actComps;
  proposedActComps.erase(proposedActComps.begin() + wormIndex);
  if (allowMultiCompWorm) {
    cout << "Worm::remove: implement me!" << endl;
  } else {
    _W_insert = 1 / this->getInsertProposalProbability(proposedActComps);
  }
  const double W_insert = _W_insert;

  ////
  //// the inverse probability of having proposed particular site
  ////
  const double W_site = this->numSites;

  ////
  //// propose a worm population and insertion method given the segment
  ////
  double W_worm = 0;
  this->proposeWormPopulationAndInsertionType(proposedActComps,
                                              goingForward ? headSegment->pop : headSegment->end->outg->pop,
                                              actComp,
                                              W_worm);

  // DEBUG: return in case of no possible worms
  if (debug && debugFrom <= this->currUpdateCount && W_worm == 0) {
      cout << settings::cout::enterRed
           << "Worm::remove: W_worm == 0  ->  EXIT"
           << settings::cout::resetStyle << endl;
    this->shutDown();
  }

  ////
  //// the inverse probability of having proposed "t1"
  ////
  const double W_t1 = this->int2time(beta) * (tMax - 2);

  ////
  //// the times and bounds
  ////
  long t1, t2Min_gcan, t2Max_gcan, t2;
  if (goingForward) {
    t1 = headSegment->end->t;   // t1 correspond to the head
    t2Min_gcan = headSegment->end->t + 1;
    t2Max_gcan = segAfterHead->end ?
                 (segAfterHead->end->outg->end ?
                  segAfterHead->end->outg->end->t - 1 :
                  this->sites[i][0]->end->t - 1 + tMax) :
                 this->sites[i][1]->end->t - 1 + tMax;
    t2 = segAfterHead->end ?           // t2 is located on the same side of the t=beta boundary?
         tailSegment->beg->t :         //   y
         tailSegment->beg->t + tMax;   //   n
  } else {
    t1 = tailSegment->beg->t;   // t1 correspond to the tail
    t2Min_gcan = tailSegment->beg->t + 1;
    t2Max_gcan = tailSegment->end ?
                 (tailSegment->end->outg->end ?
                  tailSegment->end->outg->end->t - 1 :
                  this->sites[i][0]->end->t - 1 + tMax) :
                 this->sites[i][1]->end->t - 1 + tMax;
    t2 = headSegment->beg ?            // t2 is located on the same side of the t=beta boundary?
         headSegment->end->t :         //   y
         headSegment->end->t + tMax;   //   n
  }

  ////
  //// the segment length being modified by this update
  ////
  const long modifiedLength = t2 - t1;

  // take into account canonical limitations
  if (allowMultiCompWorm) {
    // for (unsigned a = 0; a < numComps; a++) {
    //   if (this->isCanonical[a]) {
    //     // whether or not we are beyond the fixed particle number
    //     const auto beyond = (long) this->numParticles[a] - (long) tMax * this->Ns[a];
    //     cout << "Worm:remove: implement me!" << endl;

    //     if (
    //       (beyond > 0 && actPop > 0) ||
    //       (beyond < 0 && actPop < 0)
    //     ) {
    //       const long tHeadMax_can = tailSegment->beg->t < headSegment->end->t ?
    //                                 tailSegment->beg->t + tMax - 1 :
    //                                 tailSegment->beg->t - 1;

    //       tHeadMax_gcan = min(tHeadMax_gcan, tHeadMax_can);
    //     }
    //   }
    // }
  } else if (this->isCanonical[actComp]) {
    // whether or not we are beyond the fixed particle number
    const int dPop = (goingForward ? 1 : -1) * actPop;
    const auto beyond = (long) this->numParticles[actComp] + modifiedLength * dPop
                      - (long) tMax * this->Ns[actComp];

    const long t2Max_can = t1 + this->_maxNsDiff[actComp] - (-dPop) * beyond;   // OBS: -pop since insert is the anti update
    t2Max_gcan = min(t2Max_gcan, t2Max_can);
  }

  // convert to constants
  const long t2Min = t2Min_gcan,
             t2Max = t2Max_gcan;


  // in order to proceed the interval must have a nonzero length
  if (t2Min >= t2Max) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::insert: t2Min >= t2Max  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // TEST: make sure the times agree
  if (debug && test) {
    if (t1 != this->t1_prev) {
      cout << settings::cout::enterRed << "Worm::remove: t1=" << t1 << " vs t1_prev=" << this->t1_prev << settings::cout::resetStyle << endl;
    }
    if (t2Min != this->tMin_prev) {
      cout << settings::cout::enterRed << "Worm::remove: t2Min=" << t2Min << " vs tMin_prev=" << this->tMin_prev << settings::cout::resetStyle << endl;
    }
    if (t2Max != tMax_prev) {
      cout << settings::cout::enterRed << "Worm::remove: t2Max=" << t2Max << " vs tMax_prev=" << tMax_prev << settings::cout::resetStyle << endl;
    }
    if (t2 != this->t2_prev) {
      cout << settings::cout::enterRed << "Worm::remove: t2=" << t2 << " vs t2_prev=" << this->t2_prev << settings::cout::resetStyle << endl;
    }
  }


  ////
  //// parameters for the t2 distribution
  ////
  const long intLength = t2Max - t2Min;
  const double lambda = localWeightRatio_exponent,
               intervalLength = intLength * this->int2time(beta);

  // TEST: make sure that the distribution parameters agree
  if (debug && test) {
    if (intLength != this->intLength_prev) {
      cout << settings::cout::enterRed << "Worm::remove: intLength=" << intLength << " vs intLength_prev=" << this->intLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (lambda != this->lambda_prev) {
      cout << settings::cout::enterRed << "Worm::remove: lambda=" << lambda << " vs lambda_prev=" << this->lambda_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  // the inverse probability (not exponential part) of having proposed "t2"
  double W_t2_base;
  if (expDistrEnabled) {
    W_t2_base = abs(lambda) > pow(10., -10.) ?
               (1 - exp(-lambda * intervalLength)) / lambda :
               intervalLength;
  } else {
    W_t2_base = intervalLength;
  }

  ////
  //// the weight ratio of picking the worm
  ////
  const double W_select_worm = this->numWorms;

  ////
  //// the weight ratio of picking this update procedure as compared to the reverse update procedure
  ////
  const double W_remove = 1 / (1 - P_insert)
                        * (double) this->Wtot / (double) this->Wremove;

  ////
  //// the proposal distribution ratio of the proposed update: g(Z -> G) / g(G -> Z)
  //// OBS: cancellation of exponential dependence
  ////
  const double P = W_remove
                 * W_select_worm
                 * W_forwardOrBackward
                 / W_insert
                 / W_site
                 / W_worm
                 / W_t1
                 / W_t2_base;


  ////
  //// other weight contributions
  ////
  double exp_U_nn  = 0;

  // nearest neighbor interaction
  if (has_U_nn) {
    if (allowMultiCompWorm) {
      cout << "Worm::remove: TODO: implement me!" << endl;
      this->shutDown();
    } else {
      exp_U_nn += computeNNinteractionDiff(i,
                                           goingForward ? segAfterHead->pop : headSegment->pop,
                                           goingForward ? headSegment->pop : segAfterHead->pop,
                                           t1,
                                           t2);
    }
  }



  ////
  //// The acceptance ratio of going from state G to proposed state Z
  //// OBS: cancellation of exponential dependence
  ////
  ////                 ╭      W(Z)       g(Z -> G)  ╮
  //// A(G -> Z) = min │ 1, --------  ------------- │
  ////                 ╰      W(G)       g(G -> Z)  ╯
  ////
  const double R = P / abs(localWeightRatio_base)
                 * exp(   (expDistrEnabled ? 0 : localWeightRatio_exponent * this->int2time(beta) * modifiedLength)
                        + exp_U_nn * this->int2time(beta) ),
               A = min(1., R);

  // cout << "remove, beta=" << beta << ": " << R << " (" << 1/W_t2_base << ")" << endl;

  // DEBUG: the acceptance ration must always be larger than zero
  if (debug && debugFrom <= this->currUpdateCount && A < 0) {
    cout << settings::cout::enterRed << "Worm::remove: invalid value of the acceptance ration A=" << A << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  ////
  //// TEST: some final tests
  ////
  if (debug && test) {
    // segment lengths must agree
    if (this->modifiedLength_prev != modifiedLength) {
      cout << settings::cout::enterRed << "Worm::remove: modifiedLength=" << modifiedLength << " vs modifiedLength_prev=" << this->modifiedLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }

    // the acceptance ratios must agree
    if (abs(R - 1 / this->R_prev) > diffTresh) {
      cout << settings::cout::enterRed << "Worm::remove: R=" << R << " vs 1/R_prev=" << 1 / this->R_prev
           << ", abs(R - 1 / R_prev)=" << abs(R - 1 / this->R_prev) << settings::cout::resetStyle << endl;
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
      // particle number
      this->numParticlesAtSite[i * numComps + a] += modifiedLength
                                                  * (goingForward ? 1 : -1)
                                                  * this->wormPop[a];

      this->numParticles[a] += modifiedLength
                             * (goingForward ? 1 : -1)
                             * this->wormPop[a];

      // particle number squared
      for (unsigned b = a; b < numComps; b++) {
        this->numParticlesSquared[Worm::i_ab(a, b)] += modifiedLength
                                                     * (goingForward ? 1 : -1)
                                                     * (  (int) (headSegment->pop[a] * headSegment->pop[b])
                                                        - (int) (segAfterHead->pop[a] * segAfterHead->pop[b]) );
      }
    }
  } else {
    // update the particle number
    this->numParticlesAtSite[i * numComps + actComp] += modifiedLength
                                                      * (goingForward ? 1 : -1)
                                                      * actPop;

    this->numParticles[actComp] += modifiedLength
                                 * (goingForward ? 1 : -1)
                                 * actPop;

    // update the particle number square
    for (unsigned b = 0; b < numComps; b++) {
      const unsigned i_ab = Worm::i_ab(actComp, b);

      this->numParticlesSquared[i_ab] += modifiedLength
                                       * (goingForward ? 1 : -1)
                                       * (  (int) (headSegment->pop[actComp] * headSegment->pop[b])
                                          - (int) (segAfterHead->pop[actComp] * segAfterHead->pop[b]) );
    }
  }

  // update nearest neighbor interaction
  if (has_U_nn) this->U_nn -= exp_U_nn / tMax;


  ////
  //// DEBUG: store the worm type before it is removed
  ////
  if (debug && debugFrom <= this->currUpdateCount) {
    // previous worm type
    if (allowMultiCompWorm) {
      this->wormPop_prev = this->wormPop;
    } else {
      this->actComp_prev = actComp;
      this->actPop_prev = actPop;
    }
  }

  ////
  //// update worm
  ////
  if (goingForward) {
    // eating the tail from behind
    if (segAfterHead->end) {
      // not crossing t=beta
      this->removeSegment(i, this->headSegmentIndices[wormIndex] + 1);
      this->removeSegment(i, this->headSegmentIndices[wormIndex]);
    } else {
      // crossing t=beta
      this->removeSegment(i, this->headSegmentIndices[wormIndex] + 1);
      this->removeSegment(i,  0);
    }
  } else {
    // reversing into the tail
    if (headSegment->beg) {
      // not crossing t=0
      this->removeSegment(i, this->headSegmentIndices[wormIndex] - 1);
      this->removeSegment(i, this->headSegmentIndices[wormIndex]);   // the head segment index is reduced by
                                                                   // one from the previous segment removal
    } else {
      // crossing t=0
      this->removeSegment(i, 0);
      this->removeSegment(i, this->tailSegmentIndices[tailIndex]);
    }
  }

  // remove worm quantities
  this->headSegments.erase(this->headSegments.begin() + wormIndex);
  this->tailSegments.erase(this->tailSegments.begin() + tailIndex);
  this->headSegmentIndices.erase(this->headSegmentIndices.begin() + wormIndex);
  this->tailSegmentIndices.erase(this->tailSegmentIndices.begin() + tailIndex);
  if (allowMultiCompWorm) {
    cout << "Worm::remove: implement me" << endl;
  } else {
    this->actComps.erase(this->actComps.begin() + wormIndex);
    this->actPops.erase(this->actPops.begin() + wormIndex);
  }

  // this must occur after the worm has been removed
  this->numWorms--;

  ////
  //// DEBUG: worm weight check
  ////
  if (debug_major && debugFrom <= this->currUpdateCount) {
    this->checkWormWeight(beta,
                          1 / abs(localWeightRatio_base),
                            localWeightRatio_exponent * this->int2time(beta) * modifiedLength
                          + exp_U_nn * this->int2time(beta) );
  }

  // VERBOSE: let the terminal know this update is complete
  if (verbose) cout << settings::cout::enterGreen << "...removed" << settings::cout::resetStyle << endl;

  ////
  //// DEBUG: test an insert
  ////
  if (debug && debugFrom <= this->currUpdateCount) {
    // at what site the worm once existed
    this->i_prev = i;

    // previous worm type
    // (already stored)

    // interval bounds for t2
    this->tMin_prev = t2Min;
    tMax_prev = t2Max;

    // calculate t1 and t2
    this->t1_prev = t1;
    this->t2_prev = t2;

    // parameters for the exponential distribution
    this->intLength_prev = intLength;
    this->lambda_prev = lambda;

    // how was the worm removed
    this->goingForward_prev = goingForward;

    // closely related to acceptance ratio
    this->modifiedLength_prev = modifiedLength;
    this->R_prev = R;

    // try insert the worm back
    this->insert(beta, this->getInsertProposalProbability(this->actComps), true);
  }
}