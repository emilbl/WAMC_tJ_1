#include "../Worm.h"

using namespace std;

void Worm::overtake (
  const double & beta,
  const unsigned wormIndex,
  const double & P_insert,
  const bool     test
) {
  // PROFILING
  // profiler("Worm::overtake");

  // VERBOSE: let the terminal know that a jump has begun
  if (verbose && test) {
    cout << settings::cout::enterYellow << "testing to overtake back (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  } else if (verbose) {
    cout << settings::cout::enterYellow << "about to overtake (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  }

  // store the current update procedure
  constexpr unsigned updateNum = 7;
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
  const auto segAfterHead = headSegment->end->outg;
  const auto actComp = allowMultiCompWorm ? 0 : this->actComps[wormIndex];
  const auto actPop = allowMultiCompWorm ? 0 : this->actPops[wormIndex];

  ////
  //// propose going forward or backwards
  ////
  bool goingForwards;
  if (debug && test) {
    goingForwards = ! this->goingForward_prev;
  } else {
    goingForwards = !! this->pseudoRandom.template Uint<unsigned>(0, 1);
  }

  // VERBOSE: let the terminal know what type of direction is proposed
  if (verbose) {
    array<int, numComps> wp = {};
    if (allowMultiCompWorm) {
      wp = this->wormPop;
    } else {
      wp[actComp] = actPop;
    }

    cout << settings::cout::enterYellow << "...by going "
         << (goingForwards ? "forwards " : "backwards ") << wp
         << "..." << settings::cout::resetStyle << endl;
  }

  ////
  //// containing segment of the new head
  ////
  const auto contSeg = goingForwards ?                       // going forwards in time?
                       (segAfterHead->end ?                  //   y: segment after head is crossing t=beta boundary?
                        segAfterHead->end->outg :            //     n
                        this->sites[i][1]) :                 //     y
                       (headSegment->beg ?                   //   n: head segment is crossing t=0 boundary?
                        headSegment->beg->inco :             //     n
                        this->sites[i].back()->beg->inco);   //     y

  // in order to proceed the population after the proposed update must not be negative
  if (
    allowMultiCompWorm ?
    Worm::wouldHaveInvalidPop(contSeg->pop, this->wormPop, goingForwards ? 1 : -1) :
    Worm::wouldHaveInvalidPop(contSeg->pop, actComp, actPop, goingForwards ? 1 : -1)
  ) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::overtake: segment would have invalid population  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  ////
  //// interval bounds for the new head time
  ////
  long _tHeadMin = goingForwards ?                               // going forwards in time?
                   (segAfterHead->end ?                          //   y: crossed t=beta boundary locating the container segment?
                    contSeg->beg->t + 1 :                        //     n
                    contSeg->beg->t + 1 + tMax) :                //     y
                   (contSeg->beg ?                               //   n: container segment is crossing t=0 boundary?
                    (headSegment->beg ?                          //     n: crossed t=0 boundary locating the container segment?
                     contSeg->beg->t + 1 :                       //       n
                     contSeg->beg->t + 1 - tMax) :               //       y
                    this->sites[i].back()->beg->t + 1 - tMax),   //     y
       _tHeadMax = goingForwards ?                            // going forwards in time?
                   (contSeg->end ?                            //   y: container segment crossing t=beta boundary?
                    (segAfterHead->end ?                      //     n: crossed t=beta boundary locating the container segment?
                     contSeg->end->t - 1 :                    //       n
                     contSeg->end->t - 1  + tMax) :           //       y
                    this->sites[i][0]->end->t - 1 + tMax) :   //     y
                   (headSegment->beg ?                        //  n: crossed t=0 boundary locating the container segment?
                    contSeg->end->t - 1 :                     //    n
                    contSeg->end->t - 1 - tMax);              //    y

  // take into account canonical limitations
  if (allowMultiCompWorm) {
    cout << "Worm::timeShift: implement me!" << endl;
  } else if (this->isCanonical[actComp]) {
    // whether or not we are beyond the fixed particle number
    const auto beyond = (long) this->numParticles[actComp] - (long) tMax * this->Ns[actComp];

    if (goingForwards) {
      // upper limit
      const long tHeadMax_can = headSegment->end->t + this->_maxNsDiff[actComp] - actPop * beyond;
      _tHeadMax = min(_tHeadMax, tHeadMax_can);
    } else {
      // lower limit
      const long tHeadMin_can = headSegment->end->t - this->_maxNsDiff[actComp] - actPop * beyond;
      _tHeadMin = max(_tHeadMin, tHeadMin_can);
    }
  }

  ////
  //// if there is a NN interaction we cannot allow for the head to travel more than tMax
  //// otherwise the overlapping part of the interaction between site i and j wont be computed correctly
  ////
  if (has_U_nn) {
    if (goingForwards) _tHeadMax = min(headSegment->end->t + tMax, _tHeadMax);
    else               _tHeadMin = max(headSegment->end->t - tMax, _tHeadMin);
  }

  // convert to constants
  const long tHeadMin = _tHeadMin,
             tHeadMax = _tHeadMax;

  // in order to proceed the interval must have a nonzero length
  if (tHeadMin >= tHeadMax) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::overtake: tHeadMin >= tHeadMax  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  ////
  //// the part of the weight ratio of the proposed local
  //// update which will depend on the new head time
  ////
  double localWeightRatioNew_base = 1,
         localWeightRatioNew_exponent = 0;
  if (allowMultiCompWorm) {
    this->H.potenDiff(i,
                      this->wormPop,
                      goingForwards ? 1 : -1,
                      localWeightRatioNew_base,
                      localWeightRatioNew_exponent);
    if (has_U) {
      localWeightRatioNew_exponent += this->H.interDiff(contSeg->pop,
                                                        this->wormPop,
                                                        goingForwards ? 1 : -1);
    }
  } else {
    this->H.potenDiff(i,
                      actComp,
                      actPop,
                      goingForwards ? 1 : -1,
                      localWeightRatioNew_base,
                      localWeightRatioNew_exponent);
    if (has_U) {
      localWeightRatioNew_exponent += this->H.interDiff(contSeg->pop,
                                                        actComp,
                                                        actPop,
                                                        goingForwards ? 1 : -1);
    }
  }

  ////
  //// parameters for the tHead distribution
  ////
  const long intLength = tHeadMax - tHeadMin;
  const double intervalLength = intLength * this->int2time(beta),
               lambda = localWeightRatioNew_exponent;

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
    const double dt = expDistrEnabled ?
                      this->pseudoRandom.template Exp<double>(intervalLength, lambda) :
                      this->pseudoRandom.template U<double>(0, intervalLength);

    tHead = goingForwards ?
            tHeadMin + round(dt / this->int2time(beta)) :
            tHeadMax - round(dt / this->int2time(beta));
  }

  // the inverse probability (not exponential part) of having proposed this new head time
  double PtHead_base;
  if (expDistrEnabled) {
    PtHead_base = abs(lambda) > pow(10., -10.) ?
                  (1 - exp(-lambda * intervalLength)) / lambda :
                  intervalLength;
  } else {
    PtHead_base = intervalLength;
  }

  // in order to proceed the new head time must not be located on the time boundary
  if (tHead == 0 || tHead == tMax) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::overtake: tHead == 0 || tHead == tMax  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // DEBUG: make sure that tHead is inside the bounds
  if (debug && debugFrom <= this->currUpdateCount) {
    if (tHead < tHeadMin) {
      cout << settings::cout::enterRed << "Worm::overtake: tHead=" << tHead << " < " << tHeadMin << "=tHeadMin" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tHead > tHeadMax) {
      cout << settings::cout::enterRed << "Worm::overtake: tHead=" << tHead << " > " << tHeadMax << "=tHeadMax" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// take into account the periodicity in time
  ////
  const long tHeadRaw = tHead;
  this->timeModulo(tHead);

  ////
  //// interval bounds for the old (current) head time distribution
  ////
  long _tHeadMin_old = goingForwards && this->sites[i].size() == 3 ?   // special case?
                       tHead + 1 :                                     //   y
                       (headSegment->beg ?                             //   n: is the head segment crossing the t=0 boundary?
                        headSegment->beg->t + 1 :                      //     n
                        this->sites[i].back()->beg->t + 1 - tMax),     //     y
       _tHeadMax_old = ( ! goingForwards) && this->sites[i].size() == 3 ?   // special case?
                       tHead - 1 :                                          //   y
                       (segAfterHead->end ?                                 //   n: is the segment after the head corssing the t=beta boundary?
                        segAfterHead->end->t - 1 :                          //     n
                        this->sites[i][0]->end->t - 1 + tMax);              //     y

  // take into account canonical limitations
  if (allowMultiCompWorm) {
    cout << "Worm::timeShift: implement me!" << endl;
  } else if (this->isCanonical[actComp]) {
    // whether or not we are beyond the fixed particle number
    const auto beyond = (long) this->numParticles[actComp] - (long) tMax * this->Ns[actComp];

    if ( ! goingForwards) {
      // upper limit
      const long tHeadMax_old_can = headSegment->end->t + this->_maxNsDiff[actComp] - actPop * beyond;
      _tHeadMax_old = min(_tHeadMax_old, tHeadMax_old_can);
    } else {
      // lower limit
      const long tHeadMin_old_can = headSegment->end->t - this->_maxNsDiff[actComp] - actPop * beyond;
      _tHeadMin_old = max(_tHeadMin_old, tHeadMin_old_can);
    }
  }

  ////
  //// if there is a NN interaction we cannot allow for the head to travel more than tMax
  //// otherwise the overlapping part of the interaction between site i and j wont be computed correctly
  ////
  if (has_U_nn) {
    if ( ! goingForwards) {
      // reaching current configuration by going forward
      const long originShift = tHead > _tHeadMin_old ? 0 : tMax;
      _tHeadMax_old = min(tHead + originShift, _tHeadMax_old);
    } else {
      // reaching current configuration by going backward
      const long originShift = tHead > _tHeadMax_old ? -tMax : 0;
      _tHeadMin_old = max(tHead + originShift, _tHeadMin_old);
    }
  }

  // convert to constants
  const long tHeadMin_old = _tHeadMin_old,
             tHeadMax_old = _tHeadMax_old;

  // in order to proceed the interval must have a nonzero length
  if (tHeadMin_old >= tHeadMax_old) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::overtake: tHeadMin_old >= tHeadMin_old  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // TEST: make sure that the bounds agree
  if (debug && test) {
    if (tHeadMin_old != this->tMin_prev) {
      cout << settings::cout::enterRed << "Worm::overtake: tHeadMin_old=" << tHeadMin_old << " vs tMin_prev=" << this->tMin_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tHeadMax_old != tMax_prev){
      cout << settings::cout::enterRed << "Worm::overtake: tHeadMax_old=" << tHeadMax_old << " vs tMax_prev=" << tMax_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }
  ////
  //// if the we are overtaking (or being overtaken by) a head or tail
  ////
  const auto boundaryNode = goingForwards ? contSeg->beg : contSeg->end;
  const bool overtakingDiscontinuity = boundaryNode->conn == nullptr;

  ////
  //// which head or tail are we overtaking
  ////
  bool overtakingHead = false;
  unsigned overtakenWormIndex = this->numWorms;
  // shared_ptr<Segment> passedHeadOrTail = nullptr;
  if (overtakingDiscontinuity) {
    // search in head segments
    for (unsigned wi = 0; wi < this->numWorms; wi++) {
      if (this->headSegments[wi]->end == boundaryNode) {
        overtakingHead = true;
        overtakenWormIndex = wi;
        break;
      }
    }

    // search in tail segments segments
    if (overtakenWormIndex == this->numWorms) {
      for (unsigned wi = 0; wi < this->numWorms; wi++) {
        if (this->tailSegments[wi]->beg == boundaryNode) {
          overtakingHead = false;
          overtakenWormIndex = wi;
          break;
        }
      }
    }
  }

  ////
  //// the part of the weight ratio of the proposed local
  //// update which will depend on the old (current) head time
  ////
  double localWeightRatioOld_base = 1,
         localWeightRatioOld_exponent = 0;
  if (allowMultiCompWorm) {
    cout << "Worm::overtake: TODO: implement me!" << endl;
    this->shutDown();
    // this->H.potenDiff(i,
    //                   this->wormPop,
    //                   goingForwards ? 1 : -1,
    //                   localWeightRatioOld_base,
    //                   localWeightRatioOld_exponent);
    // if (has_U) {
    //   localWeightRatioOld_exponent += this->H.interDiff(goingForwards ? segAfterHead->pop : headSegment->pop,
    //                                                     this->wormPop,
    //                                                     goingForwards ? 1 : -1);
    // }
    // this->H.discoDiffHead(contSeg->pop,
    //                       headSegment->pop,
    //                       this->wormPop,
    //                       goingForwards ? 1 : 0,
    //                       localWeightRatioOld_base,
    //                       localWeightRatioOld_exponent);

    // ////
    // //// whether or not we overtake the tail or simply a jump
    // ////
    // if (tailContributionModified) {
    //   // the tail contribution is modified
    //   this->H.discoDiffTail(tailSegment->beg->inco->pop,   // the population of the segment going to contain the tail head
    //                         tailSegment->beg->inco->pop,   // the current population of the segment
    //                         this->wormPop,                 // worm population
    //                         goingForwards ? 1 : -1,        // add or subtract the worm population of the "new population" (0 do nothing)
    //                         localWeightRatioOld_base,
    //                         localWeightRatioOld_exponent);
    // } else {
    //   // the contribution from a jump is modified

    //   // lattice site index at the other side of the jump
    //   const unsigned i2 = (goingForwards ? contSeg->beg : contSeg->end)->conn->inco->siteIndex;

    //   this->H.kinetDiff(i,                                                                // lattice site index 1
    //                     i2,                                                               // lattice site index 2
    //                     (goingForwards ? contSeg->beg->inco : contSeg)->pop,              // current population at lattice site indexed 1 before jump
    //                     (goingForwards ? contSeg : contSeg->end->outg)->pop,              // current population at lattice site indexed 1 after jump
    //                     (goingForwards ? contSeg->beg : contSeg->end)->conn->inco->pop,   // current population at lattice site indexed 2 before jump
    //                     (goingForwards ? contSeg->beg : contSeg->end)->conn->outg->pop,   // current population at lattice site indexed 2 after jump
    //                     this->wormPop,                                                    // worm population
    //                     goingForwards ? 1 : -1,                                           // add or subtract worm to/from site indexed 1 before jump
    //                     goingForwards ? 1 : -1,                                           // add or subtract worm to/from site after 1 before jump
    //                     0,
    //                     0,
    //                     localWeightRatioOld_base);
    // }
  } else {
    this->H.potenDiff(i,
                      actComp,
                      actPop,
                      goingForwards ? 1 : -1,
                      localWeightRatioOld_base,
                      localWeightRatioOld_exponent);
    if (has_U) {
      localWeightRatioOld_exponent += this->H.interDiff(goingForwards ? segAfterHead->pop : headSegment->pop,
                                                        actComp,
                                                        actPop,
                                                        goingForwards ? 1 : -1);
    }

    this->H.discoDiffHead(contSeg->pop[actComp],
                          headSegment->pop[actComp],
                          actPop,
                          goingForwards ? 1 : 0,
                          localWeightRatioOld_base,
                          localWeightRatioOld_exponent);

    ////
    //// whether or not we overtake a head/tail
    ////
    if (overtakingDiscontinuity) {
      if (actComp == this->actComps[overtakenWormIndex]) {
        // the other discontinuity is affected
        this->H.discoDiffHeadOrTail(boundaryNode->inco->pop[actComp],
                                    boundaryNode->outg->pop[actComp],
                                    actPop,
                                    goingForwards ? 1 : -1,
                                    localWeightRatioOld_base);
      }
    } else {
      // the contribution from a jump is modified

      // lattice site index at the other side of the jump
      const unsigned i2 = (goingForwards ? contSeg->beg : contSeg->end)->conn->inco->siteIndex;
      this->H.kinetDiff(i,
                        i2,
                        (goingForwards ? contSeg->beg->inco : contSeg)->pop[actComp],
                        (goingForwards ? contSeg : contSeg->end->outg)->pop[actComp],
                        (goingForwards ? contSeg->beg : contSeg->end)->conn->inco->pop[actComp],
                        (goingForwards ? contSeg->beg : contSeg->end)->conn->outg->pop[actComp],
                        actComp,
                        actPop,
                        goingForwards ? 1 : -1,
                        goingForwards ? 1 : -1,
                        0,
                        0,
                        localWeightRatioOld_base);
    }
  }

  ////
  //// parameters for the old (current) tHead distribution
  ////
  const long intLength_old = tHeadMax_old - tHeadMin_old;
  const double intervalLength_old = intLength_old * this->int2time(beta),
               lambda_old = -localWeightRatioOld_exponent;

  // TEST: make sure that the distribution parameters agree
  if (debug && test) {
    if (intLength_old != this->intLength_prev) {
      cout << settings::cout::enterRed << "Worm::overtake: intLength_old=" << intLength_old << " vs intLength_prev=" << this->intLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (lambda_old != this->lambda_prev) {
      cout << settings::cout::enterRed << "Worm::overtake: lambda_old=" << lambda_old << " vs lambda_prev=" << this->lambda_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  // the inverse probability (not exponential part) of having proposed the old (current) head time
  double PtHead_old_base;
  if (expDistrEnabled) {
    PtHead_old_base = abs(lambda_old) > pow(10., -10.) ?
                  (1 - exp(-lambda_old * intervalLength_old)) / lambda_old :
                  intervalLength_old;
  } else {
    PtHead_old_base = intervalLength_old;
  }

  ////
  //// the segment lengths being modified by this update
  ////
  // the segment length on the other side of the discontinuity
  const long modifiedLengthOld = goingForwards ?                          // going forwards in time?
                                 contSeg->beg->t - headSegment->end->t    //   y
                                 + (segAfterHead->end ? 0 : tMax) :       //   add one period of time if the segment is crossing t=beta boundary
                                 headSegment->end->t - contSeg->end->t    //   n
                                 + (headSegment->beg ? 0 : tMax);         //   add one period of time if the segment is crossing t=0 boundary

  // the segment length on the current side of the discontinuity
  const long modifiedLengthNew = goingForwards ?             // going forwards in time?
                                 tHeadRaw - tHeadMin + 1 :   //   y
                                 tHeadMax + 1 - tHeadRaw;    //   n

  ////
  //// the proposal distribution ratio of the proposed update: g(t) / g(t')
  //// OBS: cancellation of exponential dependence
  ////
  const double P = PtHead_base
                 / PtHead_old_base;


  ////
  //// other weight contributions
  ////
  // double base_other = 1,
  //        exp_other  = 0;

  // nearest neighbor interaction
  if (has_U_nn) {
    if (allowMultiCompWorm) {
      cout << "Worm::overtake: TODO: implement me!" << endl;
      this->shutDown();
    } else {
      cout << "Worm::overtake: TODO: implement me!" << endl;
      this->shutDown();
      // exp_other += computeNNinteractionDiff(i,
      //                                       contSeg->pop,
      //                                       propActComp,
      //                                       propActPop * (subtract ? -1 : 1),
      //                                       t1,
      //                                       t2);
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
  const double R = P
                 * abs(localWeightRatioOld_base / localWeightRatioNew_base)
                 * (expDistrEnabled ? 1 : exp( - (  localWeightRatioOld_exponent * modifiedLengthOld
                                                  + localWeightRatioNew_exponent * modifiedLengthNew )
                                               * this->int2time(beta)) ),
               A = min(1., R);

  // DEBUG: the acceptance ration must always be larger than zero
  if (debug && debugFrom <= this->currUpdateCount && A < 0) {
    cout << settings::cout::enterRed << "Worm::overtake: invalid value of the acceptance ration A=" << A << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  ////
  //// TEST: some final tests
  ////
  if (debug && test) {
    // modified lengths must agree
    if (modifiedLengthOld != this->modifiedLength2_prev) {
      cout << settings::cout::enterRed << "Worm::overtake: modifiedLengthOld=" << modifiedLengthOld << " vs modifiedLengthNew_prev=" << this->modifiedLength2_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (modifiedLengthNew != this->modifiedLength1_prev) {
      cout << settings::cout::enterRed << "Worm::overtake: modifiedLengthNew=" << modifiedLengthNew << " vs modifiedLength1_prev=" << this->modifiedLength1_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }

    // the acceptance ratios must agree
    if (abs(R - 1 / this->R_prev) > diffTresh) {
      cout << settings::cout::enterRed << "Worm::overtake: R=" << R << " vs 1/R_prev=" << 1 / this->R_prev << settings::cout::resetStyle << endl;
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
  this->sign *= sgn(localWeightRatioOld_base) * sgn(localWeightRatioNew_base);
  if (allowMultiCompWorm) {
    for (unsigned a = 0; a < numComps; a++) {
      // update the particle number
      this->numParticlesAtSite[i * numComps + a] += (goingForwards ? 1 : -1)
                                                   * (modifiedLengthOld + modifiedLengthNew)
                                                   * this->wormPop[a];

      this->numParticles[a] += (goingForwards ? 1 : -1)
                             * (modifiedLengthOld + modifiedLengthNew)
                             * this->wormPop[a];

      // update the particle number square
      for (unsigned b = a; b < numComps; b++) {
        this->numParticlesSquared[Worm::i_ab(a, b)] += (goingForwards ? 1 : -1) * modifiedLengthOld
                                                     * (  (int) (headSegment->pop[a] * headSegment->pop[b])
                                                        - (int) (segAfterHead->pop[a] * segAfterHead->pop[b]) )
                                                     + modifiedLengthNew
                                                     * (  (int) (  (contSeg->pop[a] + (goingForwards ? 1 : -1) * this->wormPop[a])
                                                                 * (contSeg->pop[b] + (goingForwards ? 1 : -1) * this->wormPop[b]) )
                                                        - (int) (contSeg->pop[a] * contSeg->pop[b]) );
      }
    }
  } else {
    // update the particle number
    this->numParticlesAtSite[i * numComps + actComp] += (goingForwards ? 1 : -1)
                                                             * (modifiedLengthOld + modifiedLengthNew)
                                                             * actPop;

    this->numParticles[actComp] += (goingForwards ? 1 : -1)
                                       * (modifiedLengthOld + modifiedLengthNew)
                                       * actPop;
    // update the particle number square
    for (unsigned b = 0; b < numComps; b++) {
      const unsigned i_ab = Worm::i_ab(actComp, b);

        this->numParticlesSquared[i_ab] += (goingForwards ? 1 : -1) * modifiedLengthOld
                                         * (  (int) (headSegment->pop[actComp] * headSegment->pop[b])
                                            - (int) (segAfterHead->pop[actComp] * segAfterHead->pop[b]) )
                                         + modifiedLengthNew
                                         * (  (int) (  (contSeg->pop[actComp] + (goingForwards ? 1 : -1) * actPop)
                                                     * (contSeg->pop[b] + (goingForwards ? 1 : -1) * (actComp == b ? actPop : 0)) )
                                            - (int) (contSeg->pop[actComp] * contSeg->pop[b]) );
    }
  }

  // // update nearest neighbor interaction
  // if (has_U_nn) this->U_nn += exp_U_nn / tMax;



  ////
  //// DEBUG: store the old (current) head time before being overwritten
  ////
  if (debug && debugFrom <= this->currUpdateCount) {
    // the old (current) head time relative t=0 for the new head time
    this->t_prev = headSegment->end->t + tHead - tHeadRaw;
  }

  ////
  //// update worm
  ////
  if (goingForwards) {
    // going forward

    if (segAfterHead->end) {
      // not going across t=beta boundary looking for containing segment

      if (tHeadRaw < tMax) {
        // not going across t=beta boundary to reach new head time

        // reconnect jump or fix discontinuities
        if (contSeg->beg->conn) {
          contSeg->beg->conn->conn = headSegment->end;
          headSegment->end->conn = contSeg->beg->conn;
          contSeg->beg->conn = nullptr;

          // fix dists
          headSegment->end->dist = contSeg->beg->dist;
        } else if (overtakingDiscontinuity) {
          if (overtakingHead) {
            this->headSegments[overtakenWormIndex] = this->headSegments[overtakenWormIndex]->beg->inco;
            this->headSegmentIndices[overtakenWormIndex]--;
          } else {
            this->tailSegments[overtakenWormIndex] = this->tailSegments[overtakenWormIndex]->beg->inco;
            this->tailSegmentIndices[overtakenWormIndex]--;
          }
        }

        // shift times
        headSegment->end->t = contSeg->beg->t;
        contSeg->beg->t = tHead;

        // set population
        if (allowMultiCompWorm) {
          contSeg->beg->inco->pop = Worm::sum<unsigned, unsigned, int>(contSeg->pop, this->wormPop);
        } else {
          contSeg->beg->inco->pop = contSeg->pop;
          contSeg->beg->inco->pop[actComp] += actPop;
        }

        // make head
        this->headSegments[wormIndex] = contSeg->beg->inco;
        this->headSegmentIndices[wormIndex] += 1;
      } else {
        // going across t=beta boundary to reach new head time

        // reconnect jump or head
        if (contSeg->beg->conn) {
          contSeg->beg->conn->conn = headSegment->end;
          headSegment->end->conn = contSeg->beg->conn;
          contSeg->beg->conn = nullptr;

          // fix dists
          headSegment->end->dist = contSeg->beg->dist;
        } else if (overtakingDiscontinuity) {
          if (overtakingHead) {
            this->headSegments[overtakenWormIndex] = headSegment;
            this->headSegmentIndices[overtakenWormIndex]--;
          }
        }

        // shift times
        headSegment->end->t = contSeg->beg->t;

        // remove segment
        this->removeSegment(i, this->headSegmentIndices[wormIndex] + 1);

        // insert new segment
        const auto newHeadSeg = this->splitSegment(i, 0, tHead);

        // set population
        if (allowMultiCompWorm) {
          Worm::add(newHeadSeg->pop, this->wormPop);
        } else {
          newHeadSeg->pop[actComp] += actPop;
        }
        this->sites[i].back()->pop = newHeadSeg->pop;

        // make head
        this->headSegments[wormIndex] = newHeadSeg;
        this->headSegmentIndices[wormIndex] = 0;

      }
    } else {
      // going across t=beta boundary looking for containing segment

      // remove segment
      this->removeSegment(i, this->headSegmentIndices[wormIndex] + 1);

      // insert new segment
      const auto newHeadSeg = this->splitSegment(i, 1, tHead);

      // set population
      this->sites[i][0]->pop = headSegment->pop;
      if (allowMultiCompWorm) {
        Worm::add(newHeadSeg->pop, this->wormPop);
      } else {
        newHeadSeg->pop[actComp] += actPop;
      }

      // make head
      this->headSegments[wormIndex] = newHeadSeg;
      this->headSegmentIndices[wormIndex] = 1;
    }
  } else {
    // going backwards

    if (headSegment->beg) {
      // not going across t=0 boundary looking for containing segment

      if (tHeadRaw > 0) {
        // not going across t=beta boundary to reach new head time

        // reconnect jump or fix discontinuities
        if (contSeg->end->conn) {
          contSeg->end->conn->conn = headSegment->end;
          headSegment->end->conn = contSeg->end->conn;
          contSeg->end->conn = nullptr;

          // fix dists
          headSegment->end->dist = contSeg->end->dist;
        } else if (overtakingDiscontinuity) {
          if (overtakingHead) {
            this->headSegments[overtakenWormIndex] = this->headSegments[overtakenWormIndex]->end->outg;
            this->headSegmentIndices[overtakenWormIndex]++;
          } else {
            this->tailSegments[overtakenWormIndex] = this->tailSegments[overtakenWormIndex]->end->outg;
            this->tailSegmentIndices[overtakenWormIndex]++;
          }
        }

        // shift times
        headSegment->end->t = contSeg->end->t;
        contSeg->end->t = tHead;

        // set population
        if (allowMultiCompWorm) {
          contSeg->end->outg->pop = Worm::diff<unsigned, unsigned, int>(contSeg->pop, this->wormPop);
        } else {
          contSeg->end->outg->pop = contSeg->pop;
          contSeg->end->outg->pop[actComp] -= actPop;
        }

        // make head
        this->headSegments[wormIndex] = contSeg;
        this->headSegmentIndices[wormIndex] -= 1;
      } else {
        // going across t=beta boundary to reach new head time

        // fix tail
        if (overtakingDiscontinuity) {
          if ( ! overtakingHead) {
            this->tailSegments[overtakenWormIndex] = this->tailSegments[overtakenWormIndex]->end->outg;
          }
        }

        // remove segment
        this->removeSegment(i, 1);

        // insert new segment
        const auto newHeadSeg = this->splitSegment(i,
                                                   this->sites[i].size() - 1,
                                                   tHead);

        // set population
        if (allowMultiCompWorm) {
          Worm::subtract(contSeg->pop, this->wormPop);
        } else {
          contSeg->pop[actComp] -= actPop;
        }
        newHeadSeg->end->outg->pop = contSeg->pop;

        // make head
        this->headSegments[wormIndex] = newHeadSeg;
        this->headSegmentIndices[wormIndex] = this->sites[i].size() - 2;
      }
    } else {
      // going across t=0 boundary looking for containing segment

      // remove segment
      this->removeSegment(i, 0);

      // insert new segment
      const auto newHeadSeg = this->splitSegment(i,
                                              this->sites[i].size() - 2,
                                              tHead);

      // set population
      this->sites[i].back()->pop = this->sites[i][0]->pop;
      if (allowMultiCompWorm) {
        Worm::subtract(newHeadSeg->end->outg->pop, this->wormPop);
      } else {
        newHeadSeg->end->outg->pop[actComp] -= actPop;
      }

      // make head
      this->headSegments[wormIndex] = newHeadSeg;
      this->headSegmentIndices[wormIndex] = this->sites[i].size() - 3;
    }
  }



  ////
  //// DEBUG: worm weight check
  ////
  if (debug_major && debugFrom <= this->currUpdateCount) {
    this->checkWormWeight(beta,
                          abs(localWeightRatioOld_base / localWeightRatioNew_base),
                          - (  localWeightRatioOld_exponent * modifiedLengthOld
                             + localWeightRatioNew_exponent * modifiedLengthNew) * this->int2time(beta) );
  }

  // VERBOSE: let the terminal know this update is complete
  if (verbose) cout << settings::cout::enterGreen << "...overtaken" << settings::cout::resetStyle << endl;

  ////
  //// DEBUG: test an anti-overtake
  ////
  if (debug && debugFrom <= this->currUpdateCount) {
    // which direction was taken
    this->goingForward_prev = goingForwards;


    // how much the origin will shift due to this update
    const auto originShift = tHeadRaw > 0 ?
                             (tHeadRaw / tMax) * tMax :
                             (tHeadRaw / tMax - 1) * tMax;

    // interval bounds for tHead
    this->tMin_prev = tHeadMin - originShift;
    tMax_prev = tHeadMax - originShift;

    // cout << "------------" << endl
    //      << "tHeadRaw = " << tHeadRaw << endl
    //      << "tHead = " << tHead << endl
    //      << "tHeadMin = " << tHeadMin << endl
    //      << "tHeadMax = " << tHeadMax << endl
    //      << "originShift = " << originShift << endl
    //      << "this->tMin_prev = " << this->tMin_prev << endl
    //      << "tMax_prev = " << tMax_prev << endl
    //      << "------------" << endl;

    // the previous tHead time
    // already stored

    // parameters for the exponential distribution
    this->intLength_prev = intLength;
    this->lambda_prev = lambda;

    // closely related to acceptance ratio
    this->modifiedLength1_prev = modifiedLengthOld;
    this->modifiedLength2_prev = modifiedLengthNew;
    this->R_prev = R;

    // try going back with the anti-update
    this->overtake(beta, wormIndex, P_insert, true);
  }
}