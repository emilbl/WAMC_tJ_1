#include "../Worm.h"

using namespace std;



void Worm::proposeWormPopulationAndInsertionType (
  const vector<unsigned> &          activeComps,
  const array<unsigned, numComps> & segmentPop,
  // -------
  array<int, numComps> &            propWormPop,
  unsigned &                        propActComp,
  int &                             propActPop,
  double &                          W_worm,
  // -------
  bool &                            subtract
) const {
  ////
  //// we should not allow two worms of the same type since it would be complicated
  //// to allow for splitting and merging of worms. And by not allowing this would
  //// maybe harm statistics since they would block one another???
  ////

  if (allowMultiCompWorm) {
    cout << "Worm::proposeWormPopulationAndInsertionType: overload in this case" << endl;
  }

  ////
  //// tJ
  ////
  if (modelType == tJ) {
    if (this->maxNumWorms > 1) {
      // take into account current worms in the system
      bool _spinAllowed = true,
           _holeAllowed = true;
      if (activeComps.size() == 1) {
        if (activeComps[0] == 0)  _spinAllowed = false;
        else if (this->Ns[1] > 0) _holeAllowed = false;
      }

      // make constants
      const bool spinAllowed = _spinAllowed,
                 holeAllowed = _holeAllowed;

      // take into account current segment population
      const bool hasSpin = segmentPop[0] == 1;
      const bool hasHole = segmentPop[1] == 1;

      ////
      //// select a worm
      ////
      if (spinAllowed && holeAllowed) {
        // no worms in the system
        if (hasSpin || hasHole) {
          // spin/hole occupied segment, may only insert "spin/hole + subtract" or "anti-spin/anti-hole + add"
          const unsigned dice = allowAntiWorm ?
                                this->pseudoRandom.template Uint<unsigned>(0, 1) :
                                0;

          if (dice == 0) {
            // subtract worm
            propActComp = hasSpin ? 0 : 1;
            propActPop  = 1;
            subtract    = true;
          } else {
            // add anti-worm
            propActComp = hasSpin ? 0 : 1;
            propActPop  = -1;
            subtract    = false;
          }

          W_worm = allowAntiWorm ? 2 : 1;
        } else {
          // empty segment, may insert "spin + add", "anti-spin + subtract", "hole + add" or "anti-hole + subtract"
          const unsigned dice = allowAntiWorm ?
                                this->pseudoRandom.template Uint<unsigned>(0, 3) :
                                this->pseudoRandom.template Uint<unsigned>(0, 1);

          if (dice == 0) {
            propActComp = 0;
            propActPop  = 1;
            subtract    = false;
          } else if (dice == 1) {
            propActComp = 1;
            propActPop  = 1;
            subtract    = false;
          } else if (dice == 2) {
            propActComp = 0;
            propActPop  = -1;
            subtract    = true;
          } else if (dice == 3) {
            propActComp = 1;
            propActPop  = -1;
            subtract    = true;
          }

          W_worm = allowAntiWorm ? 4 : 2;
        }

        return;
      } else if (
        (hasSpin && ! spinAllowed) ||
        (hasHole && ! holeAllowed)
      ) {
        // impossible scenarios
        // must remove certain types of components whose worms are forbidden
        W_worm = 0;
        return;
      } else {
        // a single worm already in the system

        // which component is eligible a worm
        const unsigned allowedComp = spinAllowed ? 0 : 1;

        // if we should add or remove a particle
        const bool removeParticle = hasSpin || hasHole;

        const unsigned dice = allowAntiWorm ?
                              this->pseudoRandom.template Uint<unsigned>(0, 1) :
                              0;

        if (dice == 0) {
          propActComp = allowedComp;
          propActPop  = 1;
          subtract    = removeParticle;
        } else {
          propActComp = allowedComp;
          propActPop  = -1;
          subtract    = ! removeParticle;
        }

        W_worm = allowAntiWorm ? 2 : 1;
        return;
      }


      // propActComp = this->pseudoRandom.template Uint<unsigned>(0, 1);
      // propActPop  = allowAntiWorm ? (this->pseudoRandom.template Uint<unsigned>(0, 1) == 0 ? 1 : -1) : 1;
      // subtract    = !! this->pseudoRandom.template Uint<unsigned>(0, 1);

      // W_worm = allowAntiWorm ? 8 : 4;
      // return;
    } else {

      if (allowAntiWorm) {
        cout << "Worm::proposeWormPopulationAndInsertionType: implement allowAntiWorm" << endl;
        exit(EXIT_SUCCESS);
      }

      if (segmentPop[0] || segmentPop[1]) {
        // has population
        propActComp = segmentPop[0] ? 0 : 1;
        propActPop  = 1;
        subtract    = true;

        W_worm = 1;
        return;
      } else {
        // is empty
        const double W_spin = pow(this->H.getEta(0), 2);
        const double W_hole = pow(this->H.getEta(1), 2);

        // choose component
        propActComp = this->pseudoRandom.template U<double>(0, 1) > W_hole / (W_spin + W_hole) ? 0 : 1;
        propActPop  = 1;
        subtract    = false;

        W_worm = (W_spin + W_hole) / (propActComp == 0 ? W_spin : W_hole);
        return;
      }

    }
  }



  ////
  //// unknown model / bug
  ////
  cout << settings::cout::enterRed
       << "Worm::proposeWormPopulationAndInsertionType: ERROR: unknown model \"" << modelName << "\"  ->  EXIT"<< endl
       << "(this should never happen for implemented models)"
       << settings::cout::resetStyle;
  exit(EXIT_SUCCESS);
}

void Worm::proposeWormPopulationAndInsertionType (
  const vector<unsigned> &          activeComps,
  const array<unsigned, numComps> & segmentPop,
  const unsigned                    actComp,
  double &                          W_worm
) const {
  if (allowMultiCompWorm) {
    cout << "Worm::proposeWormPopulationAndInsertionType: overload in this case" << endl;
  }

  ////
  //// tj
  ////
  if (modelType == tJ) {
    if (this->maxNumWorms > 1) {
      // take into account current worms in the system
      bool _spinAllowed = true,
           _holeAllowed = true;
      if (activeComps.size() == 1) {
        if (activeComps[0] == 0)  _spinAllowed = false;
        else if (this->Ns[1] > 0) _holeAllowed = false;
      }

      // make constants
      const bool spinAllowed = _spinAllowed,
                 holeAllowed = _holeAllowed;

      // take into account current segment population
      const bool hasSpin = segmentPop[0] == 1;
      const bool hasHole = segmentPop[1] == 1;

      ////
      //// select a worm
      ////
      if (spinAllowed && holeAllowed) {
        // no worms in the system
        if (hasSpin || hasHole) {
          // spin/hole occupied segment, may only insert "spin/hole + subtract" or "anti-spin/anti-hole + add"
          W_worm = allowAntiWorm ? 2 : 1;
        } else {
          // empty segment, may insert "spin + add", "anti-spin + subtract", "hole + add" or "anti-hole + subtract"
          W_worm = allowAntiWorm ? 4 : 2;
        }
        return;
      }  else if (
        (hasSpin && ! spinAllowed) ||
        (hasHole && ! holeAllowed)
      ) {
        // impossible scenarios
        // must remove certain types of components whose worms are forbidden
        W_worm = 0;
        return;
      } else {
        // a single worm already in the system
        W_worm = allowAntiWorm ? 2 : 1;
        return;
      }

      // W_worm = allowAntiWorm ? 8 : 4;
      // return;

    } else {

      if (allowAntiWorm) {
        cout << "Worm::proposeWormPopulationAndInsertionType: implement allowAntiWorm" << endl;
        exit(EXIT_SUCCESS);
      }

      if (segmentPop[0] || segmentPop[1]) {
        // has population
        W_worm = 1;
        return;
      } else {
        // is empty
        const double W_spin = pow(this->H.getEta(0), 2);
        const double W_hole = pow(this->H.getEta(1), 2);

        W_worm = (W_spin + W_hole) / (activeComps[0] == 0 ? W_spin : W_hole);
        return;
      }

    }
  }



  ////
  //// unknown model / bug
  ////
  cout << settings::cout::enterRed
       << "Worm::proposeWormPopulationAndInsertionType: ERROR: unknown model \"" << modelName << "\"  ->  EXIT"<< endl
       << "(this should never happen for implemented models)"
       << settings::cout::resetStyle;
  exit(EXIT_SUCCESS);
}



void Worm::insert (
  const double & beta,
  const double & P_insert,
  const bool test
) {
  // PROFILING
  // profiler("Worm::insert");

  // VERBOSE: let the terminal know that an insertion has begun
  if (verbose && test) {
    cout << settings::cout::enterYellow << "testing to insert back..." << settings::cout::resetStyle << endl;
  } else if (verbose) {
    cout << settings::cout::enterYellow << "about to insert..." << settings::cout::resetStyle << endl;
  }

  // store the current update procedure
  constexpr unsigned updateNum = 0;
  this->currentUpdate = updateNum;

  // update statistics
  if ( ! test) this->updateStatistics[3 * updateNum]++;

  ////
  //// DEBUG: worm weight check
  ////
  if (debug_major && debugFrom <= this->currUpdateCount) this->preCheckWormWeight(beta);

  ////
  //// propose a site
  ////
  const unsigned i = debug && test ?
                     this->i_prev :
                     this->pseudoRandom.template Uint<unsigned>(0, this->numSites - 1);

  // the inverse probability of having proposing the particular site
  const double W_site = this->numSites;


  ////
  //// propose the t1 uniformly from 0 to beta
  ////
  long t1;
  if (debug && test) {
    t1 = this->t1_prev;
  } else {
    const double dt = this->pseudoRandom.template U<double>(0, 1);
    t1 = 1 + round(dt * (tMax - 2));
  }

  // the inverse probability of having proposed t1
  const double W_t1 = this->int2time(beta) * (tMax - 2);

  // look up the segment containing this time
  shared_ptr<Segment> contSeg;
  unsigned int contSegIndex;
  this->findSegment(i, t1, contSeg, contSegIndex);

  // make sure that t1 is not located at the end of segment
  if (contSeg->end && t1 == contSeg->end->t) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::insert: t1 must not be equal to any of the segment end times  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }


  ////
  //// propose a worm population and insertion method given the segment
  ////
  array<int, numComps> propWormPop = {};
  unsigned propActComp = 0;
  int propActPop = 0;
  double W_worm = 0;
  bool subtract = false;

  if (debug && test) {
    if (allowMultiCompWorm) {
      propWormPop = this->wormPop_prev;
    } else {
      propActComp = this->actComp_prev;
      propActPop = this->actPop_prev;
    }

    subtract = this->goingForward_prev;

    // get the probabilities
    this->proposeWormPopulationAndInsertionType(this->actComps, contSeg->pop, propActComp, W_worm);
  } else {
    ////
    //// model dependent population selection
    ////
    this->proposeWormPopulationAndInsertionType(this->actComps,
                                                contSeg->pop,
                                                propWormPop,
                                                propActComp,
                                                propActPop,
                                                W_worm,
                                                subtract);
  }

  // return in case of no possible worms
  if (W_worm == 0) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::insert: no possible worms  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // VERBOSE: let the terminal know what type of insertion is proposed
  if (verbose) {
    array<int, numComps> wp = {};
    if (allowMultiCompWorm) {
      wp = propWormPop;
    } else {
      wp[propActComp] = propActPop;
    }

    cout << settings::cout::enterYellow
         << "...by " << (subtract ? "subtracting " : "adding ") << wp << "..."
         << settings::cout::resetStyle << endl;
  }

  // make sure there will not be a negative population on the new segment
  if (
    allowMultiCompWorm ?
    Worm::wouldHaveInvalidPop(contSeg->pop, propWormPop, subtract ? -1 : 1) :
    Worm::wouldHaveInvalidPop(contSeg->pop, propActComp, propActPop, subtract ? -1 : 1)
  ) {
    // cout << "insert: would have invalid pop!" << endl;
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::insert: segment would have invalid population  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // DEBUG: make sure that 0 < t1 < beta
  if (debug && debugFrom <= this->currUpdateCount) {
    if (t1 <= 0) {
      cout << settings::cout::enterRed << "Worm::insert: t1=" << t1 << " <= 0" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (t1 >= tMax) {
      cout << settings::cout::enterRed << "Worm::insert: t1=" << t1 << " >= " << tMax << "=tMax" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// the t2 interval bounds
  ////
  long t2Min_gcan = t1 + 1,
       t2Max_gcan = this->sites[i].size() == 1 ?              // completely flat site?
                    t2Min_gcan + tMax - 2 :                   //   y ("-2" in order for t2 != t1)
                    (contSeg->end ?                           //   n: segment crossing t=beta?
                     contSeg->end->t - 1 :                    //     n
                     this->sites[i][0]->end->t - 1 + tMax);   //     y

  // take into account canonical limitations
  if (allowMultiCompWorm) {
    cout << "Worm::insert: Implement me" << endl;
    // for (unsigned a = 0; a < numComps; a++) {
    //   if (this->isCanonical[a]) {
    //     // whether or not we are beyond the fixed particle number
    //     const auto beyond = (long) this->numParticles[a] - (long) tMax * this->Ns[a];
    //     cout << "Worm:insert: implement me!" << endl;

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
  } else if (this->isCanonical[propActComp]) {
    // whether or not we are beyond the fixed particle number
    const auto beyond = (long) this->numParticles[propActComp] - (long) tMax * this->Ns[propActComp];
    const int dPop = (subtract ? -1 : 1) * propActPop;

    const long t2Max_can = t1 + this->_maxNsDiff[propActComp] - dPop * beyond;
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

  // TEST: make sure that the bounds agree
  if (debug && test) {
    if (t2Min != this->tMin_prev) {
      cout << settings::cout::enterRed << "Worm::insert: t2Min=" << t2Min << " vs tMin_prev=" << this->tMin_prev << settings::cout::resetStyle << endl;
    }
    if (t2Max != tMax_prev) {
      cout << settings::cout::enterRed << "Worm::insert: t2Max=" << t2Max << " vs tMax_prev=" << tMax_prev << settings::cout::resetStyle << endl;
    }
  }



  ////
  //// the weight ratio of the proposed local update: W(G) / W(Z)
  ////
  double localWeightRatio_base = 1,
         localWeightRatio_exponent = 0;
  if (allowMultiCompWorm) {
    this->H.potenDiff(i,
                      propWormPop,
                      subtract ? -1 : 1,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
    if (has_U) {
      localWeightRatio_exponent += this->H.interDiff(contSeg->pop,
                                                     propWormPop,
                                                     subtract ? -1 : 1);
    }
    this->H.discoDiff(contSeg->pop,
                      propWormPop,
                      subtract ? -1 : 1,
                      beta,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
  } else {
    this->H.potenDiff(i,
                      propActComp,
                      propActPop,
                      subtract ? -1 : 1,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
    if (has_U) {
      localWeightRatio_exponent += this->H.interDiff(contSeg->pop,
                                                     propActComp,
                                                     propActPop,
                                                     subtract ? -1 : 1);
    }
    this->H.discoDiff(contSeg->pop[propActComp],
                      propActComp,
                      propActPop,
                      subtract ? -1 : 1,
                      beta,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
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
      cout << settings::cout::enterRed << "Worm::insert: intLength=" << intLength << " vs intLength_prev=" << this->intLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (lambda != this->lambda_prev) {
      cout << settings::cout::enterRed << "Worm::insert: lambda=" << lambda << " vs lambda_prev=" << this->lambda_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// propose the t2 using an exponential distribution
  ////
  long t2;
  if (debug && test) {
    t2 = this->t2_prev;
  } else {
    ////
    //// what type of distribution should be used
    ////
    if (expDistrEnabled) {
      const double dt = this->pseudoRandom.template Exp<double>(intervalLength, lambda);
      t2 = t2Min + round(dt / this->int2time(beta));
    } else {
      const double dt = this->pseudoRandom.template U<double>(0, intervalLength);
      t2 = t2Min + round(dt / this->int2time(beta));
    }
  }

  // the inverse probability (not exponential part) of having proposed t2
  double W_t2_base;
  if (expDistrEnabled) {
    W_t2_base = abs(lambda) > pow(10., -10.) ?
                (1 - exp(-lambda * intervalLength)) / lambda :
                intervalLength;
  } else {
    W_t2_base = intervalLength;
  }

  // DEBUG: make sure that t2 is inside the bounds
  if (debug && debugFrom <= this->currUpdateCount) {
    if (t2 < t2Min) {
      cout << settings::cout::enterRed << "Worm::insert: t2=" << t2 << " < " << t2Min << "=t2Min" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (t2 > t2Max) {
      cout << settings::cout::enterRed << "Worm::insert: t2=" << t2 << " > " << t2Max << "=t2Max" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// the segment length being modified by this update
  ////
  const long modifiedLength = t2 - t1;

  ////
  //// take into account the periodicity in time
  ////
  const long t2Raw = t2;
  this->timeModulo(t2);

  // will not allow for segment to end/begin at boundary
  if (t2 == 0) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::insert: t2 == 0  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  ////
  //// inverse probability of having proposed a direction in which to remove the worm
  ////
  double _W_forwardOrBackward = 1;
  if (this->sites[i].size() == 1) {
    // empty site
    _W_forwardOrBackward = 2;
  } else {
    // loop through all tails to figure
    for (unsigned wi = 0; wi < this->numWorms; wi++) {

      // filter out possible worms
      if (
           this->tailSegments[wi]->siteIndex == i
        && this->actComps[wi] == propActComp
      ) {

        if (subtract) {
          // head in front of tail
          if (
               this->tailSegments[wi]->beg == (contSeg->beg ? contSeg->beg : this->sites[i].back()->beg)
            && this->tailSegments[wi]->beg->inco->pop[propActComp] == contSeg->pop[propActComp] + (subtract ? -1 : 1) * propActPop
          ) {
            _W_forwardOrBackward = 2;
            break;
          }
        } else {
          // head behind tail
          if (
               this->tailSegments[wi]->beg == (contSeg->end ? contSeg->end : this->sites[i][0]->end)
            && this->tailSegments[wi]->pop[propActComp] == contSeg->pop[propActComp] + (subtract ? -1 : 1) * propActPop
            ){
            _W_forwardOrBackward = 2;
            break;
          }
        }

      }
    }
  }
  const double W_forwardOrBackward = _W_forwardOrBackward;


  ////
  //// the weight ratio of picking this update procedure as compared to the reverse update procedure
  ////
  auto propActComps = this->actComps;
  propActComps.push_back(propActComp);
  if (allowMultiCompWorm) {
    cout << "Worm::insert: implement me" << endl;
  }
  const double W_remove = 1 / (1 - this->getInsertProposalProbability(propActComps))
                       * (double) this->Wtot / (double) this->Wremove;

  ////
  //// the weight ratio of picking the worm
  ////
  const double W_select_worm = this->numWorms + 1;

  // ////
  // //// insert proposal weight
  // ////
  // double _Pinsert = 1;
  // if (allowMultiCompWorm) {
  //   cout << "Worm::insert: implement me!" << endl;
  // } else {
  //   _Pinsert = 1 / getInsertProposalProbability(this->actComps);
  // }
  // const double Pinsert = _Pinsert;

  ////
  //// the proposal distribution ratio of the proposed update: g(G -> Z) / g(Z -> G)
  //// OBS: cancellation of exponential dependence
  ////
  const double P = (1 / P_insert) // actual probability and not weight
                 * W_site
                 * W_worm
                 * W_t1
                 * W_t2_base
                 / W_remove
                 / W_select_worm
                 / W_forwardOrBackward;

  ////
  //// other weight contributions
  ////
  double exp_U_nn  = 0;

  // nearest neighbor interaction
  if (has_U_nn) {
    if (allowMultiCompWorm) {
      cout << "Worm::insert: TODO: implement me!" << endl;
      this->shutDown();
    } else {
      exp_U_nn += computeNNinteractionDiff(i,
                                           contSeg->pop,
                                           propActComp,
                                           propActPop * (subtract ? -1 : 1),
                                           t1,
                                           t2);
    }
  }


  ////
  //// The acceptance ratio of going from state Z to proposed state G
  //// OBS: cancellation of exponential dependence
  ////
  ////                 ╭      W(G)       g(G -> Z)  ╮
  //// A(Z -> G) = min │ 1, --------  ------------- │
  ////                 ╰      W(Z)       g(Z -> G)  ╯
  ////
  const double R = abs(localWeightRatio_base) * P
                 * exp( - (expDistrEnabled ? 0 : localWeightRatio_exponent * this->int2time(beta) * modifiedLength)
                        - exp_U_nn * this->int2time(beta) ),
               A = min(1., R);

  // cout << "insert, beta=" << beta << ": " << R << endl;

  // DEBUG: the acceptance ration must always be larger than zero
  if (debug && debugFrom <= this->currUpdateCount && A < 0) {
    cout << settings::cout::enterRed << "Worm::insert: invalid value of the acceptance ration A=" << A << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  ////
  //// TEST: some final tests
  ////
  if (debug && test) {
    // modified lengths must agree
    if (this->modifiedLength_prev != modifiedLength) {
      cout << settings::cout::enterRed << "Worm::insert: modifiedLength=" << modifiedLength << " vs modifiedLength_prev=" << this->modifiedLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }

    // the acceptance ratios must agree
    if (abs(R - 1 / this->R_prev) > diffTresh) {
      cout << settings::cout::enterRed << "Worm::insert: R=" << R << " vs 1/R_prev=" << 1 / this->R_prev
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
      // update the particle number per site
      this->numParticlesAtSite[i * numComps + a] += modifiedLength
                                                  * (subtract ? -1 : 1)
                                                  * propWormPop[a];

      // update the particle number
      this->numParticles[a] += modifiedLength
                             * (subtract ? -1 : 1)
                             * propWormPop[a];

      // update the particle number square
      for (unsigned b = a; b < numComps; b++) {
        this->numParticlesSquared[Worm::i_ab(a, b)] += modifiedLength
                                                     * (  (int) ((contSeg->pop[a] + (subtract ? -propWormPop[a] : propWormPop[a]))
                                                                 * (contSeg->pop[b] + (subtract ? -propWormPop[b] : propWormPop[b])))
                                                        - (int) (contSeg->pop[a] * contSeg->pop[b]) );
      }
    }
  } else {
    // update the particle number
    this->numParticlesAtSite[i * numComps + propActComp] += modifiedLength
                                                          * (subtract ? -1 : 1)
                                                          * propActPop;

    // update the particle number
    this->numParticles[propActComp] += modifiedLength
                                     * (subtract ? -1 : 1)
                                     * propActPop;

    // update the particle number square
    for (unsigned b = 0; b < numComps; b++) {
      const unsigned i_ab = Worm::i_ab(propActComp, b);

      if (b == propActComp) {
        this->numParticlesSquared[i_ab] += modifiedLength
                                         * propActPop
                                         * (  propActPop
                                            + (subtract ? -1 : 1)
                                            * 2 * (int) contSeg->pop[propActComp] );
      } else {
        this->numParticlesSquared[i_ab] += modifiedLength
                                         * (subtract ? -1 : 1)
                                         * (int) contSeg->pop[b] * propActPop;
      }
    }
  }

  // update nearest neighbor interaction
  if (has_U_nn) this->U_nn += exp_U_nn / tMax;



  ////
  //// update worm
  ////
  if (allowMultiCompWorm) {
    this->wormPop = propWormPop;
  } else {
    this->actComps.push_back(propActComp);
    this->actPops.push_back(propActPop);
  }


  if (t2 > t1) {
    // the new worm is not crossing the time boundary

    // need to figure out what segment is the actual container
    if (contSeg->beg && contSeg->beg->t > t2) {
      // the worm is inserted on the other side of the t=beta boundary compared to the original containing segment
      contSegIndex = 0;
    } else if (contSeg->end && contSeg->end->t < t1) {
      // the worm is inserted on the other side of the t=0 boundary compared to the original containing segment
      contSegIndex = this->sites[i].size() - 1;
    }

    // split containing segments
    this->splitSegment(i, contSegIndex, t1);
    this->splitSegment(i, contSegIndex + 1, t2);

    // modify population according to worm
    if (allowMultiCompWorm) {
      Worm::addOrSubtract(this->sites[i][contSegIndex + 1]->pop,
                          propWormPop,
                          subtract ? -1 : 1);
    } else {
      this->sites[i][contSegIndex + 1]->pop[propActComp] += (subtract ? -1 : 1) * propActPop;
    }

    // assign segment indices
    if (subtract) {
      this->headSegmentIndices.push_back(contSegIndex);
      this->tailSegmentIndices.push_back(contSegIndex + 2);
    } else {
      this->headSegmentIndices.push_back(contSegIndex + 1);
      this->tailSegmentIndices.push_back(contSegIndex + 1);
    }
  } else {
    // the new worm is crossing the time boundary

    // The actual containers
    unsigned contSegIndex1 = 0,
             contSegIndex2 = this->sites[i].size() - 1;

    // split containing segments
    this->splitSegment(i, contSegIndex1, t2);
    this->splitSegment(i, contSegIndex2 + 1, t1);

    // modify population according to worm
    if (allowMultiCompWorm) {
      Worm::addOrSubtract(this->sites[i][contSegIndex1]->pop,
                          propWormPop,
                          subtract ? -1 : 1);
      Worm::addOrSubtract(this->sites[i][contSegIndex2 + 2]->pop,
                          propWormPop,
                          subtract ? -1 : 1);
    } else {
      this->sites[i][contSegIndex1]->pop[propActComp] += (subtract ? -1 : 1) * propActPop;
      this->sites[i][contSegIndex2 + 2]->pop[propActComp] += (subtract ? -1 : 1) * propActPop;
    }

    // assign segment indices
    if (subtract) {
      this->headSegmentIndices.push_back(contSegIndex2 + 1);
      this->tailSegmentIndices.push_back(contSegIndex1 + 1);
    } else {
      this->headSegmentIndices.push_back(contSegIndex1);
      this->tailSegmentIndices.push_back(contSegIndex2 + 2);
    }
  }

  // make head and tail
  this->headSegments.push_back(this->sites[i][this->headSegmentIndices.back()]);
  this->tailSegments.push_back(this->sites[i][this->tailSegmentIndices.back()]);

  // this must occur after the new worm has been inserted
  this->numWorms++;


  ////
  //// DEBUG: worm weight check
  ////
  if (debug_major && debugFrom <= this->currUpdateCount) {
    this->checkWormWeight(beta,
                          abs(localWeightRatio_base),
                          - localWeightRatio_exponent * this->int2time(beta) * modifiedLength
                          - exp_U_nn * this->int2time(beta) );
  }

  // VERBOSE: let the terminal know this update is complete
  if (verbose) cout << settings::cout::enterGreen << "...inserted" << settings::cout::resetStyle << endl;

  ////
  //// DEBUG: test a remove
  ////
  if (debug && debugFrom <= this->currUpdateCount) {
    // how is the worm inserted
    this->subtract_prev = subtract;

    // times
    this->t1_prev = t1;
    this->t2_prev = t2Raw;

    // interval bounds
    this->tMin_prev = t2Min;
    tMax_prev = t2Max;

    // parameters for the exponential distribution
    this->intLength_prev = intLength;
    this->lambda_prev = lambda;

    // closely related to acceptance ratio
    this->modifiedLength_prev = modifiedLength;
    this->R_prev = R;

    // try removing the worm
    this->remove(beta, this->numWorms - 1, this->getInsertProposalProbability(this->actComps), true);
  }
}