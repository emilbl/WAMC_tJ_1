#include "../Worm.h"

using namespace std;

void Worm::jump (
  const double & beta,
  const unsigned wormIndex,
  const double & P_insert,
  const bool     test
) {
  // PROFILING
  // profiler("Worm::jump");

  // VERBOSE: let the terminal know that a jump has begun
  if (verbose && test) {
    cout << settings::cout::enterYellow << "testing to jump back (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  } else if (verbose) {
    cout << settings::cout::enterYellow << "about to jump (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  }

  // store the current update procedure
  constexpr unsigned updateNum = 2;
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
  const unsigned j = headSegment->siteIndex;
  const auto actComp = allowMultiCompWorm ? 0 : this->actComps[wormIndex];
  const auto actPop = allowMultiCompWorm ? 0 : this->actPops[wormIndex];

  ////
  //// fetch NN indices list
  ////
  vector<unsigned> NNs;
  vector<vector<double> > dists;
  this->lattice.getNNsAndDists(j, NNs, dists);

  // TEST: make sure that the previous jumped to site is in the list
  if (debug && test) {
    if (find(NNs.begin(), NNs.end(), this->i_prev) == NNs.end()) {
      cout << settings::cout::enterRed << "Worm::jump: i_prev=" << this->i_prev << " not in NNs=" << NNs << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// propose a NN site to jump to
  ////
  unsigned i;
  vector<double> dist;
  if (debug && test) {
    i = this->i_prev;
    // no need to assign a distance here since it will never be used anyways
  } else {
    const unsigned NNindex = this->pseudoRandom.template Uint<unsigned>(0, NNs.size() - 1);
    i = NNs[NNindex];
    dist = dists[NNindex];
  }

  // the inverse probability of having proposing the NN site
  const double Psite = NNs.size();

  // TEST: make sure that the lattice indices agree
  if (debug && test) {
    if (j != this->j_prev || i != this->i_prev) {
      cout << settings::cout::enterRed
           << "Worm::jump: (i, j)=" << array<unsigned, 2>{{i, j}} << " vs (i, j)_prev="
           << array<unsigned, 2>{{this->i_prev, this->j_prev}}
           << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// find the segment at the proposed site which contains the head time
  ////
  unsigned NNsegmentIndex;
  shared_ptr<Segment> NNsegment;
  this->findSegment(i, headSegment->end->t, NNsegment, NNsegmentIndex);

  // in order to proceed there must be a positive population on the new segment
  if (
    allowMultiCompWorm ?
    Worm::wouldHaveInvalidPop(NNsegment->pop, this->wormPop, 1) :
    Worm::wouldHaveInvalidPop(NNsegment->pop, actComp, actPop, 1)
  ) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::jump: segment would have invalid population  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // in order to proceed the the segment ends are not allowed to coincide in time
  // since this would imply the head to be located on top of a jump
  if (NNsegment->end && NNsegment->end->t == headSegment->end->t) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::jump: NNsegment->end->t == headSegment->end->t  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  ////
  //// interval bounds for the jump time
  ////
  const long tJumpMin_i = NNsegment->beg ?                              // crossing t=0 boundary_i?
                          NNsegment->beg->t + 1 :                       //   n
                          (this->sites[i].size() == 1 ?                 //   y: completely flat site?
                           1 - tMax :                                   //     y
                           this->sites[i].back()->beg->t + 1 - tMax),   //     n
             tJumpMin_j = headSegment->beg ?                          // crossing t=0 boundary_j?
                          headSegment->beg->t + 1 :                   //    n
                          this->sites[j].back()->beg->t + 1 - tMax,   //    y
             tJumpMin = max(tJumpMin_i, tJumpMin_j),
             tJumpMax = headSegment->end->t - 1;

  // in order to proceed the interval must have a nonzero length
  if (tJumpMin >= tJumpMax) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::jump: tJumpMin >= tJumpMax  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // TEST: make sure that the bounds agree
  if (debug && test) {
    if (tJumpMin_i != this->tMini_prev) {
      cout << settings::cout::enterRed << "Worm::jump: tJumpMin_i=" << tJumpMin_i << " vs tMini_prev=" << this->tMini_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tJumpMin_j != this->tMinj_prev){
      cout << settings::cout::enterRed << "Worm::jump: tJumpMin_j=" << tJumpMin_j << " vs tMinj_prev=" << this->tMinj_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tJumpMin != this->tMin_prev) {
      cout << settings::cout::enterRed << "Worm::jump: tJumpMin=" << tJumpMin << " vs tMin_prev=" << this->tMin_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tJumpMax != tMax_prev) {
      cout << settings::cout::enterRed << "Worm::jump: tJumpMax=" << tJumpMax << " vs tMax_prev=" << tMax_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

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
      localWeightRatio_exponent += this->H.interDiff(NNsegment->pop,
                                                     headSegment->pop,
                                                     this->wormPop);
    }
    this->H.discoDiffHead(NNsegment->pop,
                          headSegment->pop,
                          this->wormPop,
                          1,
                          localWeightRatio_base,
                          localWeightRatio_exponent);
    this->H.kinetDiff(i,
                      j,
                      NNsegment->pop,
                      headSegment->pop,
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
      localWeightRatio_exponent += this->H.interDiff(NNsegment->pop,
                                                     headSegment->pop,
                                                     actComp,
                                                     actPop);
    }
    this->H.discoDiffHead(NNsegment->pop[actComp],
                          headSegment->pop[actComp],
                          actPop,
                          1,
                          localWeightRatio_base,
                          localWeightRatio_exponent);
    this->H.kinetDiff(i,
                      j,
                      NNsegment->pop[actComp],
                      headSegment->pop[actComp],
                      actComp,
                      actPop,
                      localWeightRatio_base);
  }

  ////
  //// parameters for the tJump distribution
  ////
  const long intLength = tJumpMax - tJumpMin;
  const double lambda = localWeightRatio_exponent,
               intervalLength = intLength * this->int2time(beta);

  if (isnan(lambda)) {
    cout << settings::cout::enterRed << "Worm::jump: lambda=" << lambda << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  // TEST: make sure that the distribution parameters agree
  if (debug && test) {
    if (intLength != this->intLength_prev) {
      cout << settings::cout::enterRed << "Worm::jump: intLength=" << intLength << " vs intLength_prev=" << this->intLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (lambda != this->lambda_prev) {
      cout << localWeightRatio_exponent << endl;
      cout << settings::cout::enterRed << "Worm::jump: lambda=" << lambda << " vs lambda_prev=" << this->lambda_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// propose a jump time
  ////
  long tJump;
  if (debug && test) {
    tJump = this->t_prev;
  } else {
    ////
    //// what type of distribution should be used
    ////
    if (expDistrEnabled) {
      const double dt = this->pseudoRandom.template Exp<double>(intervalLength, lambda);
      tJump = tJumpMax - round(dt / this->int2time(beta));
    } else {
      const double dt = this->pseudoRandom.template U<double>(0, intervalLength);
      tJump = tJumpMin + round(dt / this->int2time(beta));
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

  // DEBUG: make sure that tJump is inside the bounds
  if (debug && debugFrom <= this->currUpdateCount) {
    if (tJump < tJumpMin) {
      cout << settings::cout::enterRed << "Worm::jump: tJump=" << tJump << " < " << tJumpMin << "=tJumpMin" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tJump > tJumpMax) {
      cout << settings::cout::enterRed << "Worm::jump: tJump=" << tJump << " > " << tJumpMax << "=tJumpMax" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// the segment length being modified by this update
  ////
  const long modifiedLength = headSegment->end->t - tJump;

  ////
  //// take into account the periodicity in time
  ////
  const long tJumpRaw = tJump;
  this->timeModulo(tJump);

  // in order to proceed the jump time must not occur on the time boundary
  if (tJump == 0) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::jump: tJump = 0  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  ////
  //// the weight ratio of picking this update procedure as compared to the reverse update procedure
  ////
  const double Pweight = (double) this->WantiJump / (double) this->Wjump;

  ////
  //// the proposal distribution ratio of the proposed update: g(jump -> no jump) / g(no jump -> jump)
  //// OBS: cancellation of exponential dependence
  ////
  const double P_base = Pweight
                      * Psite
                      * PtJump_base;



  ////
  //// other weight contributions
  ////
  double exp_U_nn = 0;

  // nearest neighbor interaction
  if (has_U_nn) {
    if (allowMultiCompWorm) {
      cout << "Worm::jump: TODO: implement me!" << endl;
      this->shutDown();
    } else {
      exp_U_nn += this->computeNNinteractionDiff(i,
                                                 NNsegment->pop,
                                                 actComp,
                                                 actPop,
                                                 j,
                                                 headSegment->pop,
                                                 actComp,
                                                 - actPop,
                                                 tJump,
                                                 headSegment->end->t);
    }
  }


  ////
  //// The acceptance ratio
  //// OBS: cancellation of exponential dependence
  ////
  ////                          ╭       W(jump)      g(jump -> no jump)  ╮
  //// A(no jump -> jump) = min │ 1, -------------  -------------------- │
  ////                          ╰      W(no jump)    g(no jump -> jump)  ╯
  ////
  const double R = abs(localWeightRatio_base) * P_base
                 * exp( - (expDistrEnabled ? 0 : localWeightRatio_exponent * this->int2time(beta) * modifiedLength)
                        - exp_U_nn * this->int2time(beta) ),
               A = min(1., R);

  // DEBUG: the acceptance ration must always be larger than zero
  if (debug && debugFrom <= this->currUpdateCount && A < 0) {
    cout << settings::cout::enterRed << "Worm::jump: invalid value of the acceptance ration A=" << A << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  ////
  //// TEST: some final tests
  ////
  if (debug && test) {
    // modified lengths must agree
    if (this->modifiedLength_prev != modifiedLength) {
      cout << settings::cout::enterRed << "Worm::jump: modifiedLength=" << modifiedLength << " vs modifiedLength_prev=" << this->modifiedLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }

    // the acceptance ratios must agree
    if (abs(R - 1 / this->R_prev) > diffTresh) {
      cout << settings::cout::enterRed << "Worm::jump: R=" << R << " vs 1/R_prev=" << 1 / this->R_prev << settings::cout::resetStyle << endl;
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
  this->numKinks++;
  if (allowMultiCompWorm) {
    for (unsigned a = 0; a < numComps; a++) {
      // jumps
      this->numJumps[a] += !! this->wormPop[a];   // only 0 or 1 (hole anti-jump creates a particle jump)

      // flow
      transform(dist.begin(),
                dist.end(),
                this->flow.begin() + (j * numComps + a) * this->numDims,
                this->flow.begin() + (j * numComps + a) * this->numDims,
                [&] (const double & b, const double & c) {
                  return c + (this->wormPop[a] == 0 ? 0 : (this->wormPop[a] * b));
                });
      transform(dist.begin(),
                dist.end(),
                this->flow.begin() + (i * numComps + a) * this->numDims,
                this->flow.begin() + (i * numComps + a) * this->numDims,
                [&] (const double & b, const double & c) {
                  return c + (this->wormPop[a] == 0 ? 0 : (this->wormPop[a] * b ));
                });

      // particle number
      this->numParticlesAtSite[i * numComps + a] += modifiedLength * this->wormPop[a];
      this->numParticlesAtSite[j * numComps + a] -= modifiedLength * this->wormPop[a];

      // update the particle number square
      for (unsigned b = a; b < numComps; b++) {
        this->numParticlesSquared[Worm::i_ab(a, b)] += modifiedLength
                                                     * (  this->wormPop[b] * ((int) NNsegment->pop[a] - (int) headSegment->pop[a])
                                                        + this->wormPop[a] * ((int) NNsegment->pop[b] - (int) headSegment->pop[b])
                                                        + 2 * this->wormPop[a] * this->wormPop[b]);
      }

      // winding number
      if (this->lattice.boundaryCrossed(i, j)) {
        const auto W = this->lattice.boundaryCrossings(i, j);
        if (this->wormPop[a] > 0) {
          Worm::add(this->numWinds[a], W);
        } else if (this->wormPop[a] < 0) {
          Worm::subtract(this->numWinds[a], W);
        }
      }
    }

  } else {
    // jumps
    this->numJumps[actComp] += !! actPop;   // only 0 or 1 (hole anti-jump creates a particle jump)

    // particle number
    this->numParticlesAtSite[i * numComps + actComp] += modifiedLength * actPop;
    this->numParticlesAtSite[j * numComps + actComp] -= modifiedLength * actPop;

    // update the particle number square
    for (unsigned b = 0; b < numComps; b++) {
      const unsigned i_ab = Worm::i_ab(actComp, b);

      if (b == actComp) {
        this->numParticlesSquared[i_ab] += modifiedLength
                                         * 2 * actPop
                                         * (  (int) NNsegment->pop[actComp]
                                            - (int) headSegment->pop[actComp]
                                            + actPop );

      } else {
        this->numParticlesSquared[i_ab] += modifiedLength
                                         * actPop
                                         * ((int) NNsegment->pop[b] - (int) headSegment->pop[b]);
      }
    }

    // winding number
    if (this->lattice.boundaryCrossed(i, j)) {
      const auto W = this->lattice.boundaryCrossings(i, j);
      Worm::addOrSubtract(this->numWinds[actComp],
                          W,
                          actPop > 0 ? 1 : -1);
    }

    // flow
    transform(dist.begin(),
              dist.end(),
              this->flow.begin() + (j * numComps + actComp) * this->numDims,
              this->flow.begin() + (j * numComps + actComp) * this->numDims,
              [&] (const double & a, const double & b) {
                return actPop > 0 ? b + a : b - a;
              });
    transform(dist.begin(),
              dist.end(),
              this->flow.begin() + (i * numComps + actComp) * this->numDims,
              this->flow.begin() + (i * numComps + actComp) * this->numDims,
              [&] (const double & a, const double & b) {
                return actPop > 0 ? b + a : b - a;
              });
  }

  // update nearest neighbor interaction
  if (has_U_nn) this->U_nn += exp_U_nn / tMax;


  ////
  //// update worm
  ////
  if (tJump < headSegment->end->t) {
    // jump occurring at the same side of t=0 boundary

    // insert new segments
    const auto segBeforeJump_i = this->splitSegment(i,
                                                    NNsegmentIndex,
                                                    tJump),
               newHeadSegment = this->splitSegment(i,
                                                   NNsegmentIndex + 1,
                                                   headSegment->end->t);

    // shift time
    headSegment->end->t = tJump;

    // link nodes
    headSegment->end->conn = segBeforeJump_i->end;
    segBeforeJump_i->end->conn = headSegment->end;

    // fix population imbalance
    if (allowMultiCompWorm) {
      this->add(newHeadSegment->pop, this->wormPop);
    } else {
      newHeadSegment->pop[actComp] += actPop;
    }

    // add distances to the new nodes
    headSegment->end->dist = dist;
    transform(dist.begin(),
              dist.end(),
              segBeforeJump_i->end->dist.begin(),
              negate<double>());

    // update head
    this->headSegments[wormIndex] = newHeadSegment;
    this->headSegmentIndices[wormIndex] = NNsegmentIndex + 1;
  } else {
    // jump occurring on the other side of t=0 boundary

    // insert new segments
    const auto newHeadSegment = this->splitSegment(i, 0, headSegment->end->t),
               segBeforeJump_j = this->splitSegment(j, this->sites[j].size() - 1, tJump),
               segBeforeJump_i = this->splitSegment(i, this->sites[i].size() - 1, tJump);

    // link nodes
    segBeforeJump_j->end->conn = segBeforeJump_i->end;
    segBeforeJump_i->end->conn = segBeforeJump_j->end;

    // fix imbalance population
    if (allowMultiCompWorm) {
      this->subtract(segBeforeJump_j->end->outg->pop, this->wormPop);
      this->add(segBeforeJump_i->end->outg->pop, this->wormPop);
      this->add(newHeadSegment->pop, this->wormPop);
    } else {
      segBeforeJump_j->end->outg->pop[actComp] -= actPop;
      segBeforeJump_i->end->outg->pop[actComp] += actPop;
      newHeadSegment->pop[actComp] += actPop;
    }

    // remove former head segment
    this->removeSegment(j, this->headSegmentIndices[wormIndex]);

    // add distances to the new nodes
    segBeforeJump_j->end->dist = dist;
    transform(dist.begin(),
              dist.end(),
              segBeforeJump_i->end->dist.begin(),
              negate<double>());

    // update head
    this->headSegments[wormIndex] = newHeadSegment;
  }

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
  if (verbose) cout << settings::cout::enterGreen << "...jumped" << settings::cout::resetStyle << endl;

  ////
  //// DEBUG: test an anti-jump
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
    this->t_prev = tJumpRaw;

    // parameters for the exponential distribution
    this->intLength_prev = intLength;
    this->lambda_prev = lambda;

    // closely related to acceptance ratio
    this->modifiedLength_prev = modifiedLength;
    this->R_prev = R;

    // try going back with the anti-update
    this->antiJump(beta, wormIndex, P_insert, true);
  }
}