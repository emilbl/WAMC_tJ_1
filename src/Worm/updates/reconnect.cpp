#include "../Worm.h"

using namespace std;

void Worm::reconnect (
  const double & beta,
  const unsigned wormIndex,
  const double & P_insert,
  const bool     test
) {
  // PROFILING
  // profiler("Worm::reconnect");

  // VERBOSE: let the terminal know that a reconnect has begun
  if (verbose && test) {
    cout << settings::cout::enterYellow << "testing to reconnect back (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  } else if (verbose) {
    cout << settings::cout::enterYellow << "about to reconnect (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  }

  // store the current update procedure
  constexpr unsigned updateNum = 4;
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

  // TEST: make sure that the previous reconnected to site is in the list
  if (debug && test) {
    if (find(NNs.begin(), NNs.end(), this->i_prev) == NNs.end()) {
      cout << settings::cout::enterRed << "Worm::reconnect: i_prev=" << this->i_prev << " not in NNs=" << NNs << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// propose a NN site to reconnect to
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
           << "Worm::reconnect: (i, j)=" << array<unsigned, 2>{{i, j}} << " vs (i, j)_prev="
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

  if ( Worm::wouldBreakFermiPrinc(NNsegment->pop, this->wormPop) ) {
    if (verbose || (debug && test))
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::reconnect: About to insert fermion on a fermion occupied site (site: " << i << ") ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // in order to proceed there must be a positive population on the new segment
  if (
    allowMultiCompWorm ?
    Worm::wouldHaveInvalidPop(NNsegment->pop, this->wormPop, -1) :
    Worm::wouldHaveInvalidPop(NNsegment->pop, actComp, actPop, -1)
  ) {
    if (verbose)
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::reconnect: segment would have invalid population  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // in order to proceed the the segment ends are not allowed to coincide in time
  // since this would imply the head to be located on top of a jump
  if (NNsegment->end && NNsegment->end->t == headSegment->end->t) {
    if (verbose)
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::reconnect: NNsegment->end->t == headSegment->end->t  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  ////
  //// store some more variables
  ////
  // part of NN segment connected to an end (if any) and segment
  const auto NNsegmentHead = NNsegment->end ?     // crossing t=beta?
                             NNsegment :          //   n
                             this->sites[i][0];   //   y

  // segment after worm head
  const auto segAfterHead = headSegment->end->outg;

  ////
  //// interval bounds for the reconnection time
  ////
  const long tReconMax_i = NNsegment->end ?                      // crossing t=beta boundary?
                           NNsegmentHead->end->t - 1 :           //   n
                           (this->sites[i].size() == 1 ?         //   y: completely flat site?
                            2 * tMax - 1 :                       //     y
                            NNsegmentHead->end->t - 1 + tMax),   //     n
             tReconMax_j = segAfterHead->end ?                     // crossing t=beta boundary?
                           segAfterHead->end->t - 1 :              //   n
                           this->sites[j][0]->end->t - 1 + tMax,   //   y
             tReconMax = min(tReconMax_i, tReconMax_j),
             tReconMin = headSegment->end->t + 1;

  // in order to proceed the interval must have a nonzero length
  if (tReconMax <= tReconMin) {
    if (verbose)
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::reconnect: tReconMax <= tReconMin  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // TEST: make sure that the bounds agree
  if (debug && test) {
    if (tReconMin != this->tMin_prev) {
      cout << settings::cout::enterRed << "Worm::reconnect: tReconMin=" << tReconMin << " vs tMin_prev=" << this->tMin_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tReconMax_i != tMaxi_prev) {
      cout << settings::cout::enterRed << "Worm::reconnect: tReconMax_i=" << tReconMax_i << " vs tMaxi_prev=" << tMaxi_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tReconMax_j != tMaxj_prev){
      cout << settings::cout::enterRed << "Worm::reconnect: tReconMax_j=" << tReconMax_j << " vs tMaxj_prev=" << tMaxj_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tReconMax != tMax_prev) {
      cout << settings::cout::enterRed << "Worm::reconnect: tReconMax=" << tReconMax << " vs tMax_prev=" << tMax_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// the weight ratio of the proposed local update: W(reconnected) / W(not reconnected)
  ////
  double localWeightRatio_base = 1,
         localWeightRatio_exponent = 0;
  if (allowMultiCompWorm) {
    this->H.potenDiff(j,
                      i,
                      this->wormPop,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
    if (has_U) {
      localWeightRatio_exponent += this->H.interDiff(segAfterHead->pop,
                                                     NNsegmentHead->pop,
                                                     this->wormPop);
    }
    this->H.discoDiffHead(NNsegment->pop,
                          headSegment->pop,
                          this->wormPop,
                          0,
                          localWeightRatio_base,
                          localWeightRatio_exponent);
    this->H.kinetDiff(i,
                      j,
                      segAfterHead->pop,
                      NNsegmentHead->pop,
                      this->wormPop,
                      localWeightRatio_base);
  } else {
    this->H.potenDiff(j,
                      i,
                      actComp,
                      actPop,
                      localWeightRatio_base,
                      localWeightRatio_exponent);
    if (has_U) {
      localWeightRatio_exponent += this->H.interDiff(segAfterHead->pop,
                                                     NNsegmentHead->pop,
                                                     actComp,
                                                     actPop);
    }
    this->H.discoDiffHead(NNsegment->pop[actComp],
                          headSegment->pop[actComp],
                          actPop,
                          0,
                          localWeightRatio_base,
                          localWeightRatio_exponent);
    this->H.kinetDiff(i,
                      j,
                      segAfterHead->pop[actComp],
                      NNsegmentHead->pop[actComp],
                      actComp,
                      actPop,
                      localWeightRatio_base);
  }

  ////
  //// parameters for the tJump distribution
  ////
  const long intLength = tReconMax - tReconMin;
  const double lambda = localWeightRatio_exponent,
               intervalLength = intLength * this->int2time(beta);

  // TEST: make sure that the distribution parameters agree
  if (debug && test) {
    if (intLength != this->intLength_prev) {
      cout << settings::cout::enterRed << "Worm::reconnect: intLength=" << intLength << " vs intLength_prev=" << this->intLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (lambda != this->lambda_prev) {
      cout << settings::cout::enterRed << "Worm::reconnect: lambda=" << lambda << " vs lambda_prev=" << this->lambda_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// propose a reconnection time
  ////
  long tRecon;
  if (debug && test) {
    tRecon = this->t_prev;
  } else {
    ////
    //// what type of distribution should be used
    ////
    if (expDistrEnabled) {
      const double dt = this->pseudoRandom.template Exp<double>(intervalLength, lambda);
      tRecon = tReconMin + round(dt / this->int2time(beta));
    } else {
      const double dt = this->pseudoRandom.template U<double>(0, intervalLength);
      tRecon = tReconMin + round(dt / this->int2time(beta));
    }
  }

  // the inverse probability (not exponential part) of having proposed this reconnection time
  double PtRecon_base;
  if (expDistrEnabled) {
    PtRecon_base = abs(lambda) > pow(10., -10.) ?
                   (1 - exp(-lambda * intervalLength)) / lambda :
                   intervalLength;
  } else {
    PtRecon_base = intervalLength;
  }

  // DEBUG: make sure that tRecon is inside the bounds
  if (debug && debugFrom <= this->currUpdateCount) {
    if (tRecon < tReconMin) {
      cout << settings::cout::enterRed << "Worm::reconnect: tRecon=" << tRecon << " < " << tReconMin << "=tReconMin" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tRecon > tReconMax) {
      cout << settings::cout::enterRed << "Worm::reconnect: tRecon=" << tRecon << " > " << tReconMax << "=tReconMax" << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// the segment length being modified by this update
  ////
  long modifiedLength = tRecon - headSegment->end->t;

  ////
  //// take into account the periodicity in time
  ////
  const long tReconRaw = tRecon;
  this->timeModulo(tRecon);

  // in order to proceed the reconnection time must not occur on the time boundary
  if (tRecon == 0) {
    if (verbose)
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::reconnect: tRecon == 0  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  ////
  //// the weight ratio of picking this update procedure as compared to the reverse update procedure
  ////
  const double Pweight = (double) this->WantiReconnect / (double) this->Wreconnect;

  ////
  //// the proposal distribution ratio of the proposed update: g(recon -> not recon) / g(not recon -> recon)
  //// OBS: cancellation of exponential dependence
  ////
  const double P_base = Pweight
                      * Psite
                      * PtRecon_base;


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
      exp_U_nn += this->computeNNinteractionDiff(i,
                                                 NNsegment->pop,
                                                 actComp,
                                                 - actPop,
                                                 j,
                                                 segAfterHead->pop,
                                                 actComp,
                                                 actPop,
                                                 headSegment->end->t,
                                                 tRecon);
    }
  }


  ////
  //// The acceptance ratio
  //// OBS: cancellation of exponential dependence
  ////
  ////                             ╭        W(recon)      g(recon -> not recon)  ╮
  //// A(not recon -> recon) = min │ 1, ---------------  ----------------------- │
  ////                             ╰      W(not recon)    g(not recon -> recon)  ╯
  ////
  const double R = abs(localWeightRatio_base) * P_base
                 // * (expDistrEnabled ? 1 : exp(-localWeightRatio_exponent * this->int2time(beta) * modifiedLength)),
                 * exp( - (expDistrEnabled ? 0 : localWeightRatio_exponent * this->int2time(beta) * modifiedLength)
                        - exp_U_nn * this->int2time(beta) ),
               A = min(1., R);

  // DEBUG: the acceptance ration must always be larger than zero
  if (debug && debugFrom <= this->currUpdateCount && A < 0) {
    cout << settings::cout::enterRed << "Worm::reconnect: invalid value of the acceptance ration A=" << A << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  ////
  //// TEST: some final tests
  ////
  if (debug && test) {
    // modified lengths must agree
    if (this->modifiedLength_prev != modifiedLength) {
      cout << settings::cout::enterRed << "Worm::reconnect: modifiedLength=" << modifiedLength << " vs modifiedLength_prev=" << this->modifiedLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }

    // the acceptance ratios must agree
    if (abs(R - 1 / this->R_prev) > diffTresh) {
      cout << settings::cout::enterRed << "Worm::reconnect: R=" << R << " vs 1/R_prev=" << 1 / this->R_prev << settings::cout::resetStyle << endl;
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
      this->numJumps[a] += !! this->wormPop[a];   // only 0 or 1 (hole anti-reconnect creates a particle jump)

      // flow
      transform(dist.begin(),
                dist.end(),
                this->flow.begin() + (j * numComps + a) * this->numDims,
                this->flow.begin() + (j * numComps + a) * this->numDims,
                [&] (const double & b, const double & c) {
                  return c + (this->wormPop[a] == 0 ? 0 : (this->wormPop[a] * b ));
                });
      transform(dist.begin(),
                dist.end(),
                this->flow.begin() + (i * numComps + a) * this->numDims,
                this->flow.begin() + (i * numComps + a) * this->numDims,
                [&] (const double & b, const double & c) {
                  return c + (this->wormPop[a] == 0 ? 0 : (this->wormPop[a] * b));
                });
      // particle number
      this->numParticlesAtSite[i * numComps + a] -= modifiedLength * this->wormPop[a];
      this->numParticlesAtSite[j * numComps + a] += modifiedLength * this->wormPop[a];

      // update the particle number square
      for (unsigned b = a; b < numComps; b++) {
        this->numParticlesSquared[Worm::i_ab(a, b)] += modifiedLength
                                                     * (  this->wormPop[b] * ((int) segAfterHead->pop[a] - (int) NNsegment->pop[a])
                                                        + this->wormPop[a] * ((int) segAfterHead->pop[b] - (int) NNsegment->pop[b])
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

    cout << "Worm::reconnect: implement \"flow\" for multicomponent worms" << endl;
  } else {
    // jumps
    this->numJumps[actComp] += !! actPop;   // only 0 or 1 (hole anti-jump creates a particle jump)

    // particle number
    this->numParticlesAtSite[i * numComps + actComp] -= modifiedLength * actPop;
    this->numParticlesAtSite[j * numComps + actComp] += modifiedLength * actPop;

    // update the particle number square
    for (unsigned b = 0; b < numComps; b++) {
      const unsigned i_ab = Worm::i_ab(actComp, b);

      if (b == actComp) {
        this->numParticlesSquared[i_ab] += modifiedLength
                                         * 2 * actPop
                                         * (  (int) segAfterHead->pop[actComp]
                                            - (int) NNsegment->pop[actComp]
                                            + actPop );

      } else {
        this->numParticlesSquared[i_ab] += modifiedLength
                                         * actPop
                                         * ((int) segAfterHead->pop[b] - (int) NNsegment->pop[b]);
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
  //// try force recomputation of fermionic exchange sign
  ////
  if (allowMultiCompWorm) {
    cout << "Worm::antiReconnect: implement me!" << endl;
  } else {
    if ( ! isBosonic[actComp] && ( ! isCanonical[actComp] || Ns[actComp] > 1 )) {
      this->recomputeFermionicExchangeSign = true;
    }
  }



  ////
  //// update worm
  ////
  if (headSegment->end->t < tRecon) {
    // reconnection is occurring at the same side of t=beta

    // split segments at site i
    auto segBeforeJump_i = this->splitSegment(i,
                                              NNsegmentIndex,
                                              tRecon),
         newHeadSegment = this->splitSegment(i,
                                             NNsegmentIndex,
                                             headSegment->end->t);

    // shift time
    headSegment->end->t = tRecon;

    // link
    headSegment->end->conn = segBeforeJump_i->end;
    segBeforeJump_i->end->conn = headSegment->end;

    // fix population imbalance
    if (allowMultiCompWorm) {
      this->subtract(segBeforeJump_i->pop, this->wormPop);
    } else {
      segBeforeJump_i->pop[actComp] -= actPop;
    }

    // add distances to the new nodes
    headSegment->end->dist = dist;
    transform(dist.begin(),
              dist.end(),
              segBeforeJump_i->end->dist.begin(),
              negate<double>());

    // update head
    this->headSegments[wormIndex] = newHeadSegment;
    this->headSegmentIndices[wormIndex] = NNsegmentIndex;
  } else {
    // reconnection is occurring at the other side of t=beta

    // split
    auto segBeforeJump_i = this->splitSegment(i, 0, tRecon),
         newHeadSegment  = this->splitSegment(i, this->sites[i].size() - 1, headSegment->end->t),
         segBeforeJump_j = this->splitSegment(j, 0, tRecon);
    // remove
    this->removeSegment(j, this->headSegmentIndices[wormIndex] + 1);

    // link
    segBeforeJump_i->end->conn = segBeforeJump_j->end;
    segBeforeJump_j->end->conn = segBeforeJump_i->end;

    // fix population imbalance
    if (allowMultiCompWorm) {
      this->subtract(newHeadSegment->end->outg->pop, this->wormPop);
      this->subtract(segBeforeJump_i->pop, this->wormPop);
      this->add(segBeforeJump_j->pop, this->wormPop);
    } else {
      newHeadSegment->end->outg->pop[actComp] -= actPop;
      segBeforeJump_i->pop[actComp] -= actPop;
      segBeforeJump_j->pop[actComp] += actPop;
    }

    // add distances to the new nodes
    segBeforeJump_j->end->dist = dist;
    transform(dist.begin(),
              dist.end(),
              segBeforeJump_i->end->dist.begin(),
              negate<double>());

    // update head
    this->headSegments[wormIndex] = newHeadSegment;
    this->headSegmentIndices[wormIndex] = NNsegmentIndex + 1;
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
  if (verbose) cout << settings::cout::enterGreen << "...reconnected" << settings::cout::resetStyle << endl;

  ////
  //// DEBUG: test an anti-jump
  ////
  if (debug && debugFrom <= this->currUpdateCount) {
    // which sites were involved in the jump
    this->i_prev = i;
    this->j_prev = j;

    // interval bounds for tJump
    this->tMin_prev = tReconMin;
    tMaxi_prev = tReconMax_i;
    tMaxj_prev = tReconMax_j;
    tMax_prev = tReconMax;

    // the jump time
    this->t_prev = tReconRaw;

    // parameters for the exponential distribution
    this->intLength_prev = intLength;
    this->lambda_prev = lambda;

    // closely related to acceptance ratio
    this->modifiedLength_prev = modifiedLength;
    this->R_prev = R;

    // try going back with the anti-update
    this->antiReconnect(beta, wormIndex, P_insert, true);
  }
}