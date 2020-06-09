#include "../Worm.h"

using namespace std;

void Worm::antiJump (
  const double & beta,
  const unsigned wormIndex,
  const double & P_insert,
  const bool     test
) {
  // PROFILING
  // profiler("Worm::antiJump");

  // VERBOSE: let the terminal know that an anti-jump has begun
  if (verbose && test) {
    cout << settings::cout::enterYellow << "testing to anti-jump back (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  } else if (verbose) {
    cout << settings::cout::enterYellow << "about to anti-jump (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  }

  // store the current update procedure
  constexpr unsigned updateNum = 3;
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

  // segment before and after jump on jumped to site
  const auto segAfterJump_i  = headSegment->beg ?       // crossing t=0 boundary?
                               headSegment :            //   n
                               this->sites[i].back(),   //   y
             segBeforeJump_i = segAfterJump_i->beg->inco;

  // in order to proceed there must be a previous jump
  if ( ! segAfterJump_i->beg->conn) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::antiJump: the worm has no kinks  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // in order to proceed the previous jump must correspond the a worm of precisely the same population
  if (
    allowMultiCompWorm ?
    mismatchingPopDiff(segAfterJump_i->pop, segBeforeJump_i->pop, this->wormPop) :
    mismatchingPopDiff(segAfterJump_i->pop, segBeforeJump_i->pop, actComp, actPop)
  ) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::antiJump: this->wormPop != popDiff  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  ////
  //// store some more variables
  ////
  // the segment before the jump on the lattice site jumped from
  const auto segBeforeJump_j = segAfterJump_i->beg->conn->inco,
             segAfterJump_j = segBeforeJump_j->end->outg;

  // lattice site index jumped from
  const unsigned j = segBeforeJump_j->siteIndex;

  // TEST: make sure that the lattice indices match
  if (debug && test && ! (j == this->j_prev && i == this->i_prev)) {
    cout << settings::cout::enterRed
         << "Worm::antiJump: (i, j)=" << array<unsigned, 2>{{i, j}} << " vs (i, j)_prev="
         << array<unsigned, 2>{{this->i_prev, this->j_prev}}
         << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  ////
  //// calculate the maximum time allowed for the head segment on site j
  ////
  const long maxAllowedHeadTime = (segAfterJump_j->end ?                   // crossing the t=beta boundary?
                                   segAfterJump_j->end->t - 1 :            //   n
                                   this->sites[j][0]->end->t - 1 + tMax)   //   y
                                  - (headSegment->beg ? 0 : tMax);         // t=0 boundary crossed when reversing to latest jump

  // in order to proceed the current head time must not be greater than the maximum head time
  if (maxAllowedHeadTime < headSegment->end->t) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::antiJump: maxAllowedHeadTime < headSegment->end->t  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  ////
  //// the inverse probability of having proposing the NN site
  ////
  vector<unsigned> NNs;
  vector<vector<double> > dists;
  this->lattice.getNNsAndDists(j, NNs, dists);
  const double Psite = NNs.size();

  ////
  //// the segment population before the current jump was made
  ////
  const auto prevSegPop_i = segBeforeJump_i->pop,
             prevSegPop_j = segBeforeJump_j->pop;

  ////
  //// the weight ratio of the proposed local update: W(jump) / W(no jump)
  ////
  double localWeightRatio_base = 1,
         localWeightRatio_exponent = 0;
  if (allowMultiCompWorm) {
    this->H.potenDiff(i,
                      j,
                      this->wormPop,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
    if (has_U) {
      localWeightRatio_exponent += this->H.interDiff(prevSegPop_i,
                                                     prevSegPop_j,
                                                     this->wormPop);
    }
    this->H.discoDiffHead(prevSegPop_i,
                          prevSegPop_j,
                          this->wormPop,
                          1,
                          localWeightRatio_base,
                          localWeightRatio_exponent);
    this->H.kinetDiff(i,
                      j,
                      prevSegPop_i,
                      prevSegPop_j,
                      this->wormPop,
                      localWeightRatio_base);
  } else {
    this->H.potenDiff(i,
                      j,
                      actComp,
                      actPop,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
    if (has_U) {
      localWeightRatio_exponent += this->H.interDiff(prevSegPop_i,
                                                     prevSegPop_j,
                                                     actComp,
                                                     actPop);
    }
    this->H.discoDiffHead(prevSegPop_i[actComp],
                          prevSegPop_j[actComp],
                          actPop,
                          1,
                          localWeightRatio_base,
                          localWeightRatio_exponent);
    this->H.kinetDiff(i,
                      j,
                      prevSegPop_i[actComp],
                      prevSegPop_j[actComp],
                      actComp,
                      actPop,
                      localWeightRatio_base);
  }

  ////
  //// tJump and interval bounds
  ////
  const long tJumpMin_i = this->sites[i].size() == 3 ?                    // completely flat site?
                          1 - tMax :                                      //   y
                          (segAfterJump_i->end ?                          //   n: jump on the same side of t=0?
                           (segBeforeJump_i->beg ?                        //     y: end on the same side of t=0?
                            segBeforeJump_i->beg->t + 1 :                 //       y
                            this->sites[i].back()->beg->t + 1 - tMax) :   //       n
                           segBeforeJump_i->beg->t + 1 - tMax),           //     n
             tJumpMin_j = (segBeforeJump_j->beg ?                              // end on the same side of t=0?
                           segBeforeJump_j->beg->t + 1 :                       //   y
                           this->sites[j].back()->beg->t + 1 - tMax)           //   n
                          - (segBeforeJump_j->end->t > headSegment->end->t ?   // jump on other side of t=0?
                             tMax :                                            //   y
                             0),                                               //   n
             tJumpMin = max(tJumpMin_i, tJumpMin_j),
             tJumpMax = headSegment->end->t - 1,
             tJump = segAfterJump_i->end ?            // jumped on same side of t=0 as head?
                     segAfterJump_i->beg->t :         //    y
                     segAfterJump_i->beg->t - tMax;   //    n

  // in order to proceed the interval must have a nonzero length
  if (tJumpMin >= tJumpMax) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::jump: tJumpMin >= tJumpMax  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // TEST: make sure that the times agree
  if (debug && test) {
    if (tJumpMin_i != this->tMini_prev) {
      cout << settings::cout::enterRed << "Worm::antiJump: tJumpMin_i=" << tJumpMin_i << " vs tMini_prev=" << this->tMini_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tJumpMin_j != this->tMinj_prev) {
      cout << settings::cout::enterRed << "Worm::antiJump: tJumpMin_j=" << tJumpMin_j << " vs tMinj_prev=" << this->tMinj_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tJumpMin != this->tMin_prev) {
      cout << settings::cout::enterRed << "Worm::antiJump: tJumpMin=" << tJumpMin << " vs tMin_prev=" << this->tMin_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tJumpMax != tMax_prev) {
      cout << settings::cout::enterRed << "Worm::antiJump: tJumpMax=" << tJumpMax << " vs tMax_prev=" << tMax_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tJump != this->t_prev) {
      cout << settings::cout::enterRed << "Worm::antiJump: tJump=" << tJump << " vs tJump_prev=" << this->t_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// the segment length being modified by this update
  ////
  const long modifiedLength = headSegment->end->t - tJump;

  ////
  //// parameters for the tJump distribution
  ////
  const long intLength = tJumpMax - tJumpMin;
  const double lambda = localWeightRatio_exponent,
               intervalLength = intLength * this->int2time(beta);

  // TEST: make sure that the distribution parameters agree
  if (debug && test) {
    if (intLength != this->intLength_prev) {
      cout << settings::cout::enterRed << "Worm::antiJump: intLength=" << intLength << " vs intLength_prev=" << this->intLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (lambda != this->lambda_prev) {
      cout << settings::cout::enterRed << "Worm::antiJump: lambda=" << lambda << " vs lambda_prev=" << this->lambda_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  // the inverse probability (not exponential part) of having proposed this jump time
  double PtJump_base;
  if (expDistrEnabled) {
    PtJump_base = abs(lambda) > pow(10., -10.) ?
                  (1 - exp(-lambda * intervalLength)) / lambda :
                  intervalLength;
  } else {
    PtJump_base = intervalLength;
  }

  ////
  //// the weight ratio of picking this update procedure as compared to the reverse update procedure
  ////
  const double Pweight = (double) this->Wjump / (double) this->WantiJump;

  ////
  //// the proposal distribution ratio of the proposed update: g(jump -> no jump) / g(no jump -> jump)
  //// OBS: cancellation of exponential dependence
  ////
  const double P_base = Pweight
                      / Psite
                      / PtJump_base;


  ////
  //// other weight contributions
  ////
  double exp_U_nn = 0;

  // nearest neighbor interaction
  if (has_U_nn) {
    if (allowMultiCompWorm) {
      cout << "Worm::antiJump: TODO: implement me!" << endl;
      this->shutDown();
    } else {
      exp_U_nn += this->computeNNinteractionDiff(i,
                                                 segAfterJump_i->pop,
                                                 segBeforeJump_i->pop,
                                                 j,
                                                 segAfterJump_j->pop,
                                                 segBeforeJump_j->pop,
                                                 tJump,
                                                 headSegment->end->t);
    }
  }

  ////
  //// The acceptance ratio
  //// OBS: cancellation of exponential dependence
  ////
  ////                          ╭      W(no jump)    g(no jump -> jump)  ╮
  //// A(jump -> no jump) = min │ 1, -------------  -------------------- │
  ////                          ╰       W(jump)      g(jump -> no jump)  ╯
  ////
  const double R = P_base / abs(localWeightRatio_base)
                 * exp(   (expDistrEnabled ? 0 : localWeightRatio_exponent * this->int2time(beta) * modifiedLength)
                        + exp_U_nn * this->int2time(beta) ),
               A = min(1., R);

  // DEBUG: the acceptance ration must always be larger than zero
  if (debug && debugFrom <= this->currUpdateCount && A < 0) {
    cout << settings::cout::enterRed << "Worm::antiJump: invalid value of the acceptance ration A=" << A << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  ////
  //// TEST: some final tests
  ////
  if (debug && test) {
    // modified lengths must agree
    if (this->modifiedLength_prev != modifiedLength) {
      cout << settings::cout::enterRed << "Worm::antiJump: modifiedLength=" << modifiedLength << " vs modifiedLength_prev=" << this->modifiedLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }

    // the acceptance ratios must agree
    if (abs(R - 1 / this->R_prev) > diffTresh) {
      cout << settings::cout::enterRed << "Worm::antiJump: R=" << R << " vs 1/R_prev=" << 1 / this->R_prev << settings::cout::resetStyle << endl;
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
  this->numKinks--;
  if (allowMultiCompWorm) {
    for (unsigned a = 0; a < numComps; a++) {
      // jumps
      this->numJumps[a] -= !! this->wormPop[a];   // only 0 or 1 (hole anti-jump removes a particle jump)

      // flow
      transform(segBeforeJump_i->end->dist.begin(),
                segBeforeJump_i->end->dist.end(),
                this->flow.begin() + (j * numComps + a) * this->numDims,
                this->flow.begin() + (j * numComps + a) * this->numDims,
                [&] (const double & b, const double & c) {
                  return c + this->wormPop[a] * b;
                });
      transform(segBeforeJump_i->end->dist.begin(),
                segBeforeJump_i->end->dist.end(),
                this->flow.begin() + (i * numComps + a) * this->numDims,
                this->flow.begin() + (i * numComps + a) * this->numDims,
                [&] (const double & b, const double & c) {
                  return c + this->wormPop[a] * b;
                });
      // particle number
      this->numParticlesAtSite[j * numComps + a] += modifiedLength * this->wormPop[a];
      this->numParticlesAtSite[i * numComps + a] -= modifiedLength * this->wormPop[a];

      // update the particle number square
      for (unsigned b = a; b < numComps; b++) {
        this->numParticlesSquared[Worm::i_ab(a, b)] += modifiedLength
                                                     * (  this->wormPop[b] * ((int) segAfterJump_j->pop[a] - (int) headSegment->pop[a])
                                                        + this->wormPop[a] * ((int) segAfterJump_j->pop[b] - (int) headSegment->pop[b])
                                                        + 2 * this->wormPop[a] * this->wormPop[b]);
      }

      // winding number
      if (this->lattice.boundaryCrossed(i, j)) {
        const auto W = this->lattice.boundaryCrossings(i, j);
        if (this->wormPop[a] > 0) {
          Worm::subtract(this->numWinds[a], W);
        } else if (this->wormPop[a] < 0) {
          Worm::add(this->numWinds[a], W);
        }
      }
    }

  } else {
    // jumps
    this->numJumps[actComp] -= !! actPop;   // only 0 or 1 (hole anti-jump creates a particle jump)

    // particle number
    this->numParticlesAtSite[j * numComps + actComp] += modifiedLength * actPop;
    this->numParticlesAtSite[i * numComps + actComp] -= modifiedLength * actPop;

    // update the particle number square
    for (unsigned b = 0; b < numComps; b++) {
      const unsigned i_ab = Worm::i_ab(actComp, b);

      if (b == actComp) {
        this->numParticlesSquared[i_ab] += modifiedLength
                                         * 2 * actPop
                                         * (  (int) segAfterJump_j->pop[actComp]
                                            - (int) headSegment->pop[actComp]
                                            + actPop );

      } else {
        this->numParticlesSquared[i_ab] += modifiedLength
                                         * actPop
                                         * ((int) segAfterJump_j->pop[b] - (int) headSegment->pop[b]);
      }
    }

    // winding number
    if (this->lattice.boundaryCrossed(i, j)) {
      const auto W = this->lattice.boundaryCrossings(i, j);
      Worm::addOrSubtract(this->numWinds[actComp],
                          W,
                          actPop > 0 ? -1 : 1);
    }

    // flow
    transform(segBeforeJump_i->end->dist.begin(),
              segBeforeJump_i->end->dist.end(),
              this->flow.begin() + (j * numComps + actComp) * this->numDims,
              this->flow.begin() + (j * numComps + actComp) * this->numDims,
              [&] (const double & a, const double & b) {
                return actPop > 0 ? b + a : b - a;
              });
    transform(segBeforeJump_i->end->dist.begin(),
              segBeforeJump_i->end->dist.end(),
              this->flow.begin() + (i * numComps + actComp) * this->numDims,
              this->flow.begin() + (i * numComps + actComp) * this->numDims,
              [&] (const double & a, const double & b) {
                return actPop > 0 ? b + a : b - a;
              });
  }

  // update nearest neighbor interaction
  if (has_U_nn) this->U_nn -= exp_U_nn / tMax;



  ////
  //// update worm
  ////
  if (headSegment->beg) {
    // jump occurring at the same side of t=0 boundary

    // shift times
    segBeforeJump_j->end->t = headSegment->end->t;

    // unlink
    segBeforeJump_j->end->conn = nullptr;
    headSegment->beg->conn = nullptr;

    // if this is true, the tailSegment will be removed and thus
    // the segment to become the new tail must be assigned to this->tailSegments
    if (headSegment->beg->inco == tailSegment) {
      this->tailSegments[wormIndex] = headSegment->end->outg;
    }

    // remove former head and the previous segment
    this->removeSegment(i, this->headSegmentIndices[wormIndex]);
    this->removeSegment(i, this->headSegmentIndices[wormIndex] - 1);

    // update head
    this->findSegment(j,
                      tJump,
                      this->headSegments[wormIndex],
                      this->headSegmentIndices[wormIndex]);
  } else {
    // jump occurring on the other side of t=0 boundary

    // unlink
    segBeforeJump_j->end->conn->conn = nullptr;
    segBeforeJump_j->end->conn = nullptr;

    // insert new head
    this->splitSegment(j, 0, headSegment->end->t);

    // remove segments
    this->removeSegment(i, this->sites[i].size() - 1);
    this->removeSegment(i, 0);
    this->removeSegment(j, this->sites[j].size() - 1);

    // update head (headSegmentIndex already correct)
    this->headSegments[wormIndex] = this->sites[j][0];

    // increase population
    if (allowMultiCompWorm) {
      this->add(this->headSegments[wormIndex]->pop, this->wormPop);
    } else {
      this->headSegments[wormIndex]->pop[actComp] += actPop;
    }
  }

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
  if (verbose) cout << settings::cout::enterGreen << "...anti-jumped" << settings::cout::resetStyle << endl;

  ////
  //// DEBUG: test a jump
  ////
  if (debug && debugFrom <= this->currUpdateCount) {
    // which sites were involved in the jump
    this->i_prev = i;
    this->j_prev = j;

    // interval bounds for tJump
    this->tMini_prev = tJumpMin_i;
    this->tMinj_prev = tJumpMin_j;
    this->tMin_prev = tJumpMin;
    tMax_prev = tJumpMax;

    // the jump time
    this->t_prev = tJump;

    // parameters for the exponential distribution
    this->intLength_prev = intLength;
    this->lambda_prev = lambda;

    // closely related to acceptance ratio
    this->modifiedLength_prev = modifiedLength;
    this->R_prev = R;

    // try going back with the anti-update
    this->jump(beta, wormIndex, P_insert, true);
  }
}