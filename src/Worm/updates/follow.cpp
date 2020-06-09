#include "../Worm.h"

using namespace std;

void Worm::follow (
  const double & beta,
  const unsigned wormIndex,
  const double & P_insert,
  const bool     test
) {
  // PROFILING
  // profiler("Worm::follow");

  // VERBOSE: let the terminal know that a jump has begun
  if (verbose && test) {
    cout << settings::cout::enterYellow << "testing to follow back (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  } else if (verbose) {
    cout << settings::cout::enterYellow << "about to follow (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  }

  // store the current update procedure
  constexpr unsigned updateNum = 8;
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
  const unsigned j = headSegment->siteIndex;
  const auto segAfterHead = headSegment->end->outg;
  const auto actComp = allowMultiCompWorm ? 0 : this->actComps[wormIndex];
  const auto actPop = allowMultiCompWorm ? 0 : this->actPops[wormIndex];

  ////
  //// determine the connection (if any) in the forward and backward direction
  ////
  const auto fwdConn = segAfterHead->end ?             // segment after head crossing t=beta boundary?
                       segAfterHead->end->conn :       //    n
                       this->sites[j][0]->end->conn,   //    y
             bwdConn = headSegment->beg ?                  // head segment crossing t=0 boundary?
                       headSegment->beg->conn :            //    n
                       this->sites[j].back()->beg->conn;   //    y

  // in order to proceed at least one possible direction must be available
  if ( ! fwdConn && ! bwdConn) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::follow: no directions to follow  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // DEBUG: make sure that the only possibility of not moving
  //        in one direction is if the worm tail is blocking
  if (debug && debugFrom <= this->currUpdateCount) {
    if ( ! fwdConn) {
      // not able to follow forward
      const auto blockingNode = segAfterHead->end ? segAfterHead->end : this->sites[j][0]->end;

      // two possible ways this blocking could happen
      const bool isHead = find(this->headSegments.begin(), this->headSegments.end(), blockingNode->inco) != this->headSegments.end();
      const bool isTail = find(this->tailSegments.begin(), this->tailSegments.end(), blockingNode->outg) != this->tailSegments.end();

      if ( ! isHead && ! isTail) {
        cout << settings::cout::enterRed << "Worm::follow: not able to follow fwd even though a tail is not blocking" << settings::cout::resetStyle << endl;
        if (shutItDown) this->shutDown();
      }
    }
    if ( ! bwdConn) {
      // unable to follow backwards
      const auto blockingNode = headSegment->beg ? headSegment->beg : this->sites[j].back()->beg;

      // two possible ways this blocking could happen
      const bool isHead = find(this->headSegments.begin(), this->headSegments.end(), blockingNode->inco) != this->headSegments.end();
      const bool isTail = find(this->tailSegments.begin(), this->tailSegments.end(), blockingNode->outg) != this->tailSegments.end();

      if ( ! isHead && ! isTail) {
        cout << settings::cout::enterRed << "Worm::follow: not able to follow bwd even though the tail is not blocking" << settings::cout::resetStyle << endl;
        if (shutItDown) this->shutDown();
      }
    }
  }


  ////
  //// propose a direction
  ////
  bool goingFwd;
  if (debug && test) {
    goingFwd = ! this->goingForward_prev;
  } else {
    goingFwd = (!! bwdConn && !! fwdConn) ?
               !! this->pseudoRandom.template Uint<unsigned>(0, 1) :
               !! fwdConn;
  }

  // TEST: check whether or not the proposed direction is allowed
  if (debug && test) {
    if ((goingFwd && ! fwdConn) || (! goingFwd && ! bwdConn)) {
      cout << settings::cout::enterRed << "Worm::follow: the proposed direction is not allowed" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  // inverse probability of proposing the particular direction
  const double PfwdOrBwd = (!! fwdConn && !! bwdConn) ? 2 : 1;

  // VERBOSE: let the terminal know what type of direction is proposed
  if (verbose) {
    array<int, numComps> wp = {};
    if (allowMultiCompWorm) {
      wp = this->wormPop;
    } else {
      wp[actComp] = actPop;
    }

    cout << settings::cout::enterYellow << "...by going "
         << (goingFwd ? "forwards " : "backwards ") << wp
         << "..." << settings::cout::resetStyle << endl;
  }

  ////
  //// containing the segment of the new head an the connecting node
  ////
  const auto contSeg = goingFwd ?
                       fwdConn->outg :
                       bwdConn->inco;
  const auto connectingNode = goingFwd ?
                              fwdConn :
                              bwdConn;

  // make sure the proposed update does not produce an invalid configuration
  const auto & nextPop = goingFwd ?
                         fwdConn->conn->outg->pop :
                         headSegment->end->outg->pop;
  const auto & currPop = goingFwd ?
                         headSegment->pop :
                         bwdConn->conn->inco->pop;
  if (Worm::hasInvalidPopDiff(nextPop, currPop)) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::follow: there would be an invalid population difference  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }


  ////
  //// the proposed lattice site index
  ////
  const unsigned i = contSeg->siteIndex;

  // TEST: make sure that the lattice indices agree
  if (debug && test) {
    if (i != this->j_prev || j != this->i_prev) {
      cout << settings::cout::enterRed
           << "Worm::follow: (i, j)=" << array<unsigned, 2>{{i, j}} << " vs (j, i)_prev="
           << array<unsigned, 2>{{this->j_prev, this->i_prev}}
           << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  // in order to proceed the population after the proposed update must not be negative
  if (
    allowMultiCompWorm ?
    Worm::wouldHaveInvalidPop(contSeg->pop, this->wormPop, goingFwd ? 1 : -1) :
    Worm::wouldHaveInvalidPop(contSeg->pop, actComp, actPop, goingFwd ? 1 : -1)
  ) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::follow: segment would have negative population  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }


  ////
  //// interval bounds for the new head time
  ////
  long _tHeadMin = goingFwd ?                                      // going forwards in time?
                   (segAfterHead->end ?                            //   y: crossed t=beta boundary locating the container segment?
                    contSeg->beg->t + 1 :                          //     n
                    contSeg->beg->t + 1 + tMax) :                  //     y
                   (  (contSeg->beg ?                              //   n: container segment is crossing t=0 boundary?
                       contSeg->beg->t + 1 :                       //     n
                       this->sites[i].back()->beg->t + 1 - tMax)   //     y
                    - (headSegment->beg ?                          //   additional t=0 crossing locating the container segment?
                       0 :                                         //     n
                       tMax)  ),                                   //     y
       _tHeadMax = goingFwd ?                                  // going forwards in time?
                   (  (contSeg->end ?                          //   y: container segment crossing t=beta boundary?
                       contSeg->end->t - 1 :                   //     n
                       this->sites[i][0]->end->t - 1 + tMax)   //     y
                    + (segAfterHead->end ?                     //   additional t=beta crossing locating the container segment?
                       0 :                                     //     n
                       tMax)  ) :                              //     y
                   (headSegment->beg ?                         //  n: crossed t=0 boundary locating the container segment?
                    contSeg->end->t - 1 :                      //    n
                    contSeg->end->t - 1 - tMax);               //    y

  // take into account canonical limitations
  if (allowMultiCompWorm) {
    cout << "Worm::follow: implement canonical limitations for \"allowMultiCompWorm\"!" << endl;
  } else if (this->isCanonical[actComp]) {
    // whether or not we are beyond the fixed particle number
    const auto beyond = (long) this->numParticles[actComp] - (long) tMax * this->Ns[actComp];

    if (goingFwd) {
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
    if (goingFwd) _tHeadMax = min(headSegment->end->t + tMax, _tHeadMax);
    else          _tHeadMin = max(headSegment->end->t - tMax, _tHeadMin);
  }

  // convert to constants
  const long tHeadMin = _tHeadMin,
             tHeadMax = _tHeadMax;

  // in order to proceed the interval must have a nonzero length
  if (tHeadMin >= tHeadMax) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::follow: tHeadMin >= tHeadMax  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // TEST: make sure that the bounds agree
  if (debug && test) {
    if (tHeadMin != this->tMin_old_prev) {
      cout << settings::cout::enterRed << "Worm::follow: tHeadMin=" << tHeadMin << " vs tMin_old_prev=" << this->tMin_old_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tHeadMax != this->tMax_old_prev) {
      cout << settings::cout::enterRed << "Worm::follow: tHeadMax=" << tHeadMax << " vs tMax_old_prev=" << this->tMax_old_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }


  ////
  //// the part of the weight ratio of the proposed local
  //// update which will depend on the new head time
  ////
  double localWeightRatioNew_base = 1,
         localWeightRatioNew_exponent = 0;
  if (allowMultiCompWorm) {
    cout << "worm:follow: implement this!" << endl;
  } else {
    ////
    //// diagonal matrix elements
    ////
    this->H.potenDiff(i,
                      actComp,
                      actPop,
                      goingFwd ? 1 : -1,
                      localWeightRatioNew_base,
                      localWeightRatioNew_exponent);
    if (has_U) {
      localWeightRatioNew_exponent += this->H.interDiff(contSeg->pop,
                                                       actComp,
                                                       actPop,
                                                       goingFwd ? 1 : -1);
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

    if (isnan(lambda)) {
      cout << "follow" << endl;
      exit(EXIT_FAILURE);
    }

    const double dt = expDistrEnabled ?
                      this->pseudoRandom.template Exp<double>(intervalLength, lambda) :
                      this->pseudoRandom.template U<double>(0, intervalLength);

    tHead = goingFwd ?
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
           << "Worm::follow: tHead == 0 || tHead == tMax  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // DEBUG: make sure that tHead is inside the bounds
  if (debug && debugFrom <= this->currUpdateCount) {
    if (tHead < tHeadMin) {
      cout << settings::cout::enterRed << "Worm::follow: tHead=" << tHead << " < " << tHeadMin << "=tHeadMin" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tHead > tHeadMax) {
      cout << settings::cout::enterRed << "Worm::follow: tHead=" << tHead << " > " << tHeadMax << "=tHeadMax" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// take into account the periodicity in time
  ////
  const long tHeadRaw = tHead;
  this->timeModulo(tHead);

  ////
  //// how many possible directions one may follow from the new configuration
  ////
  const double PfwdOrBwd_old = 1   // there is always the possibility of going back
                             + (goingFwd ?
                                (contSeg->end ?
                                 !! contSeg->end->conn :
                                 !! this->sites[i][0]->end->conn) :
                                (contSeg->beg ?
                                 !! contSeg->beg->conn :
                                 !! this->sites[i].back()->beg->conn));

  ////
  //// interval bounds for the old (current) head time distribution
  ////
  long _tHeadMin_old = headSegment->beg ?                          // is the head segment crossing the t=0 boundary?
                       headSegment->beg->t + 1 :                   //   n
                       this->sites[j].back()->beg->t + 1 - tMax,   //   y
       _tHeadMax_old = segAfterHead->end ?                     // is the segment after the head crossing the t=beta boundary?
                       segAfterHead->end->t - 1 :              //   n
                       this->sites[j][0]->end->t - 1 + tMax;   //   y

  // take into account canonical limitations
  if (allowMultiCompWorm) {
    cout << "Worm::timeShift: implement me!" << endl;
  } else if (this->isCanonical[actComp]) {
    // whether or not we are beyond the fixed particle number
    const auto beyond = (long) this->numParticles[actComp] - (long) tMax * this->Ns[actComp];

    if ( ! goingFwd) {
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
    if ( ! goingFwd) {
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
           << "Worm::follow: tHeadMin_old >= tHeadMax_old  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // TEST: make sure that the bounds agree
  if (debug && test) {
    if (tHeadMin_old != this->tMin_prev) {
      cout << settings::cout::enterRed << "Worm::follow: tHeadMin_old=" << tHeadMin_old << " vs tMin_prev=" << this->tMin_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tHeadMax_old != tMax_prev){
      cout << settings::cout::enterRed << "Worm::follow: tHeadMax_old=" << tHeadMax_old << " vs tMax_prev=" << tMax_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// the part of the weight ratio of the proposed local
  //// update which will depend on the old (current) head time
  ////
  double localWeightRatioOld_base = 1,
         localWeightRatioOld_exponent = 0;
  if (allowMultiCompWorm) {
    cout << "worm:follow: implement this!" << endl;
  } else {
    ////
    //// diagonal matrix elements
    ////
    this->H.potenDiff(j,
                      actComp,
                      actPop,
                      goingFwd ? 1 : -1,
                      localWeightRatioOld_base,
                      localWeightRatioOld_exponent);
    if (has_U) {
      localWeightRatioOld_exponent += this->H.interDiff(goingFwd ? segAfterHead->pop : headSegment->pop,
                                                        actComp,
                                                        actPop,
                                                        goingFwd ? 1 : -1);
    }


    ////
    //// off-diagonal matrix elements
    ////
    if (goingFwd) {
      this->H.jumpDiff(j,
                       connectingNode->conn->inco->pop,
                       connectingNode->conn->outg->pop,
                       i,
                       connectingNode->inco->pop,
                       connectingNode->outg->pop,
                       actComp,
                       actPop,
                       localWeightRatioOld_base);
    } else {
      this->H.jumpDiff(i,
                       connectingNode->inco->pop,
                       connectingNode->outg->pop,
                       j,
                       connectingNode->conn->inco->pop,
                       connectingNode->conn->outg->pop,
                       actComp,
                       -actPop,
                       localWeightRatioOld_base);
    }


    ////
    //// worm
    ////
    this->H.discoDiffHead(contSeg->pop[actComp],
                          headSegment->pop[actComp],
                          actPop,
                          goingFwd ? 1 : 0,
                          localWeightRatioOld_base,
                          localWeightRatioOld_exponent);
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
      cout << settings::cout::enterRed << "Worm::follow: intLength_old=" << intLength_old << " vs intLength_prev=" << this->intLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (lambda_old != this->lambda_prev) {
      cout << settings::cout::enterRed << "Worm::follow: lambda_old=" << lambda_old << " vs lambda_prev=" << this->lambda_prev << settings::cout::resetStyle << endl;
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
  const long modifiedLengthOld = goingFwd ?                              // going forwards in time?
                                 contSeg->beg->t - headSegment->end->t   //   y
                                 + (segAfterHead->end ? 0 : tMax) :      //   add one period of time if the segment is crossing t=beta boundary
                                 headSegment->end->t - contSeg->end->t   //   n
                                 + (headSegment->beg ? 0 : tMax);        //   add one period of time if the segment is crossing t=0 boundary

  // the segment length on the current side of the discontinuity
  const long modifiedLengthNew = goingFwd ?                  // going forwards in time?
                                 tHeadRaw - tHeadMin + 1 :   //   y
                                 tHeadMax + 1 - tHeadRaw;    //   n

  ////
  //// the proposal distribution ratio of the proposed update: g(t) / g(t')
  //// OBS: cancellation of exponential dependence
  ////
  const double P = PfwdOrBwd
                 * PtHead_base
                 / PfwdOrBwd_old
                 / PtHead_old_base;


  ////
  //// other weight contributions
  ////
  double exp_U_nn = 0;

  // nearest neighbor interaction
  if (has_U_nn) {
    if (allowMultiCompWorm) {
      cout << "Worm::reconnect: TODO: implement me!" << endl;
      this->shutDown();
    } else {
      if (goingFwd) {
        exp_U_nn += computeNNinteractionDiff(j,
                                             headSegment->pop,
                                             headSegment->end->outg->pop,
                                             headSegment->end->t,
                                             fwdConn->t);
        exp_U_nn += computeNNinteractionDiff(i,
                                             contSeg->pop,
                                             actComp,
                                             actPop,
                                             fwdConn->t,
                                             tHead);
      } else {
        exp_U_nn += computeNNinteractionDiff(j,
                                             headSegment->end->outg->pop,
                                             headSegment->pop,
                                             bwdConn->t,
                                             headSegment->end->t);
        exp_U_nn += computeNNinteractionDiff(i,
                                             contSeg->pop,
                                             actComp,
                                             - actPop,
                                             tHead,
                                             bwdConn->t);
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
  const double R = P
                 * abs(localWeightRatioOld_base / localWeightRatioNew_base)
                 * exp( - (expDistrEnabled ? 0 : (  localWeightRatioOld_exponent * modifiedLengthOld
                                                  + localWeightRatioNew_exponent * modifiedLengthNew ) * this->int2time(beta) )
                        - exp_U_nn * this->int2time(beta) ),
               A = min(1., R);

  // DEBUG: the acceptance ration must always be larger than zero
  if (debug && debugFrom <= this->currUpdateCount && A < 0) {
    cout << settings::cout::enterRed << "Worm::follow: invalid value of the acceptance ration A=" << A << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  ////
  //// TEST: some final tests
  ////
  if (debug && test) {
    // modified lengths must agree
    if (modifiedLengthOld != this->modifiedLength2_prev) {
      cout << settings::cout::enterRed << "Worm::follow: modifiedLengthOld=" << modifiedLengthOld << " vs modifiedLengthNew_prev=" << this->modifiedLength2_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (modifiedLengthNew != this->modifiedLength1_prev) {
      cout << settings::cout::enterRed << "Worm::follow: modifiedLengthNew=" << modifiedLengthNew << " vs modifiedLength1_prev=" << this->modifiedLength1_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }

    // the acceptance ratios must agree
    if (abs(R - 1 / this->R_prev) > diffTresh) {
      cout << settings::cout::enterRed << "Worm::follow: R=" << R << " vs 1/R_prev=" << 1 / this->R_prev << settings::cout::resetStyle << endl;
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
      // number of jumps
      if (this->wormPop[a]) {
        if ((int) connectingNode->inco->pop[a] - (int) connectingNode->outg->pop[a] == this->wormPop[a]) {
          // a jump is removed
          this->numJumps[a]--;
        } else {
          // a jump is created
          this->numJumps[a]++;
        }
      }

      // update the particle number
      this->numParticlesAtSite[j * numComps + a] += (goingFwd ? 1 : -1) * modifiedLengthOld * this->wormPop[a];
      this->numParticlesAtSite[i * numComps + a] += (goingFwd ? 1 : -1) * modifiedLengthNew * this->wormPop[a];

      // update the particle number square
      for (unsigned b = a; b < numComps; b++) {
        const unsigned i_ab = Worm::i_ab(actComp, b);

        this->numParticlesSquared[i_ab] += (goingFwd ? 1 : -1) * modifiedLengthOld
                                         * (  (int) (headSegment->pop[a] * headSegment->pop[b])
                                            - (int) (segAfterHead->pop[a] * segAfterHead->pop[b]) )
                                         + modifiedLengthNew
                                         * (  (int) (  (contSeg->pop[a] + (goingFwd ? 1 : -1) * this->wormPop[a])
                                                     * (contSeg->pop[b] + (goingFwd ? 1 : -1) * this->wormPop[b]) )
                                            - (int) (contSeg->pop[a] * contSeg->pop[b]));
      }

      // winding number
      auto W = this->lattice.boundaryCrossings(i, j);
      if (this->wormPop[a] > 0) {
        Worm::add(this->numWinds[a], W);
      } else if (this->wormPop[a] < 0) {
        Worm::subtract(this->numWinds[a], W);
      }
    }

    cout << "worm:follow: implement this update flow for \"allowMultiComp\"!" << endl;
    cout << "worm:follow: implement this update flow for \"numParticles\"!" << endl;

  } else {
    // update jump and exchange
    vector<unsigned> comps;
    for (unsigned a = 0; a < numComps; a++) {
      if (nextPop[a] != currPop[a]) comps.emplace_back(a);
    }
    if (comps.size() == 2) {
      // jump -> exchange
      this->numExchanges[i_ab_external(comps[0], comps[1])]++;
      this->numJumps[comps[0] == actComp ? comps[1] : comps[0]]--;
    } else if (comps.size() == 1) {
      // exchange -> jump
      this->numExchanges[i_ab_external(comps[0], actComp)]--;
      this->numJumps[comps[0] == actComp ? comps[1] : comps[0]]++;
    } else {
      cout << "Worm::follow: the heck!?" << endl;
      exit(EXIT_SUCCESS);
    }

    // update the particle number
    this->numParticlesAtSite[j * numComps + actComp] += (goingFwd ? 1 : -1) * modifiedLengthOld * actPop;
    this->numParticlesAtSite[i * numComps + actComp] += (goingFwd ? 1 : -1) * modifiedLengthNew * actPop;
    this->numParticles[actComp] += (goingFwd ? 1 : -1) * (modifiedLengthOld + modifiedLengthNew) * actPop;

    // update the particle number square
    for (unsigned b = 0; b < numComps; b++) {
      const unsigned i_ab = Worm::i_ab(actComp, b);

      this->numParticlesSquared[i_ab] += (goingFwd ? 1 : -1) * modifiedLengthOld
                                       * (  (int) (headSegment->pop[actComp] * headSegment->pop[b])
                                          - (int) (segAfterHead->pop[actComp]      * segAfterHead->pop[b]) )
                                       + modifiedLengthNew
                                       * (  (int) (  (contSeg->pop[actComp] + (goingFwd ? 1 : -1) * actPop)
                                                   * (contSeg->pop[b]             + (goingFwd ? 1 : -1) * (actComp == b ? actPop : 0)) )
                                          - (int) (contSeg->pop[actComp] * contSeg->pop[b]));
    }

    // winding number
    if (this->lattice.boundaryCrossed(i, j)) {
      const auto W = this->lattice.boundaryCrossings(i, j);
      Worm::addOrSubtract(this->numWinds[actComp], W, actPop > 0 ? 1 : -1);
    }

    // flow
    transform(connectingNode->dist.begin(),
              connectingNode->dist.end(),
              this->flow.begin() + (j * numComps + actComp) * this->numDims,
              this->flow.begin() + (j * numComps + actComp) * this->numDims,
              [&] (const double & a, const double & b) {
                return b - a * actPop;
              });
    transform(connectingNode->dist.begin(),
              connectingNode->dist.end(),
              this->flow.begin() + (i * numComps + actComp) * this->numDims,
              this->flow.begin() + (i * numComps + actComp) * this->numDims,
              [&] (const double & a, const double & b) {
                return b - a * actPop;
              });
  }

  // update nearest neighbor interaction
  if (has_U_nn)  this->U_nn += exp_U_nn / tMax;


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
  if (goingFwd) {
    // going forwards

    if (segAfterHead->end) {
      // not going across t=beta boundary looking for containing segment

      if (tHeadRaw < tMax) {
        // not going across t=beta boundary to reach new head time

        // retrieve segment index for segment to be splitted
        unsigned newHeadSegmentIndex;
        Worm::findSegmentIndex(i, tHead, newHeadSegmentIndex);

        // split segment
        const auto newHeadSeg = this->splitSegment(i, newHeadSegmentIndex, tHead);

        // set population
        headSegment->end->outg->pop = headSegment->pop;
        if (allowMultiCompWorm) add(newHeadSeg->pop, this->wormPop);
        else                    newHeadSeg->pop[actComp] += actPop;

        // update tail before removing
        if (headSegment == tailSegment) {
          this->tailSegments[wormIndex] = headSegment->end->outg;
          this->tailSegmentIndices[wormIndex] += 1;
        }

        // remove segment
        this->removeSegment(j, this->headSegmentIndices[wormIndex]);

        // make head
        this->headSegments[wormIndex] = newHeadSeg;
        this->headSegmentIndices[wormIndex] = newHeadSegmentIndex;
      } else {
        // going across t=beta boundary to reach new head time

        // split segment
        const auto newHeadSeg = this->splitSegment(i, 0, tHead);

        // set population
        headSegment->end->outg->pop = headSegment->pop;
        if (allowMultiCompWorm) add(contSeg->pop, this->wormPop);
        else                    contSeg->pop[actComp] += actPop;
        newHeadSeg->pop = contSeg->pop;

        // update tail before removing
        if (headSegment == tailSegment) {
          this->tailSegments[wormIndex] = headSegment->end->outg;
          this->tailSegmentIndices[wormIndex] += 1;
        }

        // remove segment
        this->removeSegment(j, this->headSegmentIndices[wormIndex]);

        // make head
        this->headSegments[wormIndex] = newHeadSeg;
        this->headSegmentIndices[wormIndex] = 0;
      }
    } else {
      // going across t=beta boundary looking for containing segment
      if (tHeadRaw > 2 * tMax) {
        // also going across t=beta once more to reach new head time

        // split segment
        const auto newHeadSeg = this->splitSegment(i, 0, tHead);

        // set population
        this->sites[j][0]->pop = headSegment->pop;
        if (allowMultiCompWorm) add(newHeadSeg->pop, this->wormPop);
        else                    newHeadSeg->pop[actComp] += actPop;
        this->sites[i].back()->pop = newHeadSeg->pop;

        // remove segment
        this->removeSegment(j, this->headSegmentIndices[wormIndex] + 1);

        // make head
        this->headSegments[wormIndex] = newHeadSeg;
        this->headSegmentIndices[wormIndex] = 0;

      } else {
        // not going across t=beta once more to reach new head time

        unsigned newHeadSegmentIndex;
        findSegmentIndex(i, tHead, newHeadSegmentIndex);

        // split segment
        const auto newHeadSeg = this->splitSegment(i, newHeadSegmentIndex, tHead);

        // set population
        this->sites[j][0]->pop = headSegment->pop;
        if (allowMultiCompWorm) add(newHeadSeg->pop, this->wormPop);
        else                    newHeadSeg->pop[actComp] += actPop;

        // remove segment
        this->removeSegment(j, this->headSegmentIndices[wormIndex] + 1);

        // make head
        this->headSegments[wormIndex] = newHeadSeg;
        this->headSegmentIndices[wormIndex] = newHeadSegmentIndex;
      }
    }
  } else {
    // going backwards

    if (headSegment->beg) {
      // not going across t=0 boundary looking for containing segment

      if (tHeadRaw > 0) {
        // not going across t=beta boundary to reach the new head time

        // retrieve segment index for segment to be splitted
        unsigned newHeadSegmentIndex;
        findSegmentIndex(i, tHead, newHeadSegmentIndex);

        // split segment
        const auto newHeadSeg = this->splitSegment(i, newHeadSegmentIndex, tHead);

        // set population
        if (allowMultiCompWorm) subtract(contSeg->pop, this->wormPop);
        else                    contSeg->pop[actComp] -= actPop;

        // remove segment
        this->removeSegment(j, this->headSegmentIndices[wormIndex]);

        // make head
        this->headSegments[wormIndex] = newHeadSeg;
        this->headSegmentIndices[wormIndex] = newHeadSegmentIndex;
      } else {
        // going across t=beta boundary to reach the new head time

        const unsigned newHeadSegmentIndex = this->sites[i].size() - 1;

        // split segment
        const auto newHeadSeg = this->splitSegment(i,
                                                   newHeadSegmentIndex,
                                                   tHead);

        // set population
        if (allowMultiCompWorm) subtract(contSeg->pop, this->wormPop);
        else                    contSeg->pop[actComp] -= actPop;
        newHeadSeg->end->outg->pop = contSeg->pop;

        // remove segment
        this->removeSegment(j, this->headSegmentIndices[wormIndex]);

        // make head
        this->headSegments[wormIndex] = newHeadSeg;
        this->headSegmentIndices[wormIndex] = newHeadSegmentIndex;
      }
    } else {
      // going across t=0 boundary looking for containing segment

      if (tHeadRaw < -tMax) {
        // going across t=0 once more to reach the new head time

        const unsigned newHeadSegmentIndex = this->sites[i].size() - 1;

        // split segment
        const auto newHeadSeg = this->splitSegment(i, newHeadSegmentIndex, tHead);

        // set population
        this->sites[j].back()->pop = headSegment->end->outg->pop;
        if (allowMultiCompWorm) subtract(contSeg->pop, this->wormPop);
        else                    contSeg->pop[actComp] -= actPop;
        newHeadSeg->end->outg->pop = contSeg->pop;

        // remove segment
        this->removeSegment(j, 0);

        // make head
        this->headSegments[wormIndex] = newHeadSeg;
        this->headSegmentIndices[wormIndex] = newHeadSegmentIndex;

      } else {
        // not going across t=0 once more to reach the new head time

        unsigned newHeadSegmentIndex;
        findSegmentIndex(i, tHead, newHeadSegmentIndex);

        // split segment
        const auto newHeadSeg = this->splitSegment(i, newHeadSegmentIndex, tHead);

        // set population
        this->sites[j].back()->pop = headSegment->end->outg->pop;
        if (allowMultiCompWorm) subtract(contSeg->pop, this->wormPop);
        else                    contSeg->pop[actComp] -= actPop;

        // remove segment
        this->removeSegment(j, 0);

        // make head
        this->headSegments[wormIndex] = newHeadSeg;
        this->headSegmentIndices[wormIndex] = newHeadSegmentIndex;
      }
    }
  }

  ////
  //// DEBUG: worm weight check
  ////
  if (debug_major && debugFrom <= this->currUpdateCount) {
    this->checkWormWeight(beta,
                          abs(localWeightRatioOld_base / localWeightRatioNew_base),
                          - (  localWeightRatioOld_exponent * modifiedLengthOld
                             + localWeightRatioNew_exponent * modifiedLengthNew) * this->int2time(beta)
                          - exp_U_nn * this->int2time(beta) );
  }

  // VERBOSE: let the terminal know this update is complete
  if (verbose) cout << settings::cout::enterGreen << "...followed" << settings::cout::resetStyle << endl;

  ////
  //// DEBUG: test an anti-follow
  ////
  if (debug && debugFrom <= this->currUpdateCount) {
    // which sites were involved in the jump
    this->i_prev = i;
    this->j_prev = j;

    // which direction was taken
    this->goingForward_prev = goingFwd;

    // how much the origin will shift due to this update
    const auto originShift = tHeadRaw > 0 ?
                             (tHeadRaw / tMax) * tMax :
                             (tHeadRaw / tMax - 1) * tMax;

    // interval bounds for tHead relative tHead
    this->tMin_prev = tHeadMin - originShift;
    this->tMax_prev = tHeadMax - originShift;

    // interval bounds for old tHead relative tHead
    this->tMin_old_prev = tHeadMin_old - originShift;
    this->tMax_old_prev = tHeadMax_old - originShift;

    // parameters for the exponential distribution
    this->intLength_prev = intLength;
    this->lambda_prev    = lambda;

    // closely related to acceptance ratio
    this->modifiedLength1_prev = modifiedLengthOld;
    this->modifiedLength2_prev = modifiedLengthNew;
    this->R_prev = R;

    // try going back with the anti-update
    this->follow(beta, wormIndex, P_insert, true);
  }
}