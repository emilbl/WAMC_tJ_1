#include "../Worm.h"

using namespace std;

void Worm::antiReconnect (
  const double & beta,
  const unsigned wormIndex,
  const double & P_insert,
  const bool     test
) {
  // PROFILING
  // profiler("Worm::antiReconnect");

  // VERBOSE: let the terminal know that an anti-reconnect has begun
  if (verbose && test) {
    cout << settings::cout::enterYellow << "testing to anti-reconnect back (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  } else if (verbose) {
    cout << settings::cout::enterYellow << "about to anti-reconnect (" << wormIndex << ")..." << settings::cout::resetStyle << endl;
  }

  // store the current update procedure
  constexpr unsigned updateNum = 5;
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

  // segments on site "i" after and before reconnection
  const auto segAfterRecon_i  = headSegment->end->outg->end ?         // crossing t=beta boundary?
                                headSegment->end->outg->end->outg :   //   n
                                this->sites[i][1],                    //   y
             segBeforeRecon_i = segAfterRecon_i->beg->inco;

  // in order to proceed there must be reconnection to anti-reconnect
  if ( ! segAfterRecon_i->beg->conn) {
    if (verbose)
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::antiReconnect: ! segAfterRecon_i->beg->conn  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  // in order to proceed this reconnection must further more exactly match the population of the current worm
  if (
    allowMultiCompWorm ?
    mismatchingPopDiff(segAfterRecon_i->pop, segBeforeRecon_i->pop, this->wormPop) :
    mismatchingPopDiff(segAfterRecon_i->pop, segBeforeRecon_i->pop, actComp, actPop)
  ) {
    if (verbose)
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::antiReconnect: this->wormPop != popDiff  ->  RETURN" << settings::cout::resetStyle << endl;
    if (debug && test && shutItDown) this->shutDown(); else return;
  }

  ////
  //// store some more variables
  ////
  // lattice site index of site reconnected from
  const unsigned j = segAfterRecon_i->beg->conn->inco->siteIndex;

  // segments on site "j" after and before reconnection
  const auto segAfterRecon_j = segAfterRecon_i->beg->conn->outg,
             segBeforeRecon_j = segAfterRecon_j->beg->inco;

  // TEST: make sure that the lattice indices match
  if (debug && test && ! (j == this->j_prev && i == this->i_prev)) {
    cout << settings::cout::enterRed
         << "Worm::antiReconnect: (i, j)=" << array<unsigned, 2>{{i, j}} << " vs (i, j)_prev="
         << array<unsigned, 2>{{this->i_prev, this->j_prev}}
         << settings::cout::resetStyle << endl;
    if (shutItDown) this->shutDown();
  }

  ////
  //// calculate the minimum allowed time on "j" site which the head time must be larger than
  ////
  const long minAllowedHeadTime = (segBeforeRecon_j->beg ?                     // crossing t=0 boundary?
                                   segBeforeRecon_j->beg->t + 1 :              //   n
                                   this->sites[j].back()->beg->t + 1 - tMax)   //   y
                                + (headSegment->end->outg->end ? 0 : tMax);    // might also have passed t=beta boundary looking for reconnection

  // in order to proceed the current head time must be larger than this min time
  if (headSegment->end->t < minAllowedHeadTime) {
    if (verbose)
      cout << (test ? settings::cout::enterRed : settings::cout::enterOrange)
           << "Worm::antiReconnect: headSegment->end->t < minAllowedHeadTime  ->  RETURN" << settings::cout::resetStyle << endl;
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
  //// the segment population before the current reconnection was made
  ////
  const auto prevSegPop_i = segAfterRecon_i->pop,
                            prevSegPop_j = segAfterRecon_j->pop,
                            prevHeadSegPop = segAfterRecon_i->beg->conn->inco->pop;

  ////
  //// the weight ratio of the proposed local update: W(reconnected) / W(not reconnected)
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
      localWeightRatio_exponent += this->H.interDiff(prevSegPop_j,
                                                     prevSegPop_i,
                                                     this->wormPop);
    }
    this->H.discoDiffHead(prevSegPop_i,
                          prevHeadSegPop,
                          this->wormPop,
                          0,
                          localWeightRatio_base,
                          localWeightRatio_exponent);
    this->H.kinetDiff(i,
                      j,
                      prevSegPop_j,
                      prevSegPop_i,
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
      localWeightRatio_exponent += this->H.interDiff(prevSegPop_j,
                                                     prevSegPop_i,
                                                     actComp,
                                                     actPop);
    }
    this->H.discoDiffHead(prevSegPop_i[actComp],
                          prevHeadSegPop[actComp],
                          actPop,
                          0,
                          localWeightRatio_base,
                          localWeightRatio_exponent);
    this->H.kinetDiff(i,
                      j,
                      prevSegPop_j[actComp],
                      prevSegPop_i[actComp],
                      actComp,
                      actPop,
                      localWeightRatio_base);
  }

  ////
  //// calculate the reconnection time bounds
  ////
  const long tReconMin = headSegment->end->t + 1,
             tReconMax_i = this->sites[i].size() == 3 ?             // was a completely flat site?
                           2 * tMax - 1 :                           //   y
                           segAfterRecon_i->end ?                   //   n: segAfterRecon crossing t=beta boundary?
                           (headSegment->end->outg->end ?           //     n: segBeforeRecon crossing t=beta boundary?
                            segAfterRecon_i->end->t - 1 :           //       n
                            segAfterRecon_i->end->t - 1 + tMax) :   //       y
                           this->sites[i][0]->end->t - 1 + tMax,    //     y
             tReconMax_j = segAfterRecon_j->end ?                   // segAfterRecon crossing t=beta boundary?
                           (headSegment->end->outg->end ?           //   n: segBeforeRecon crossing t=beta boundary?
                            segAfterRecon_j->end->t - 1 :           //     n
                            segAfterRecon_j->end->t - 1 + tMax) :   //     y
                           this->sites[j][0]->end->t - 1 + tMax,    //   y
             tReconMax = min(tReconMax_i, tReconMax_j),
             tRecon = segBeforeRecon_i->beg ?           // crossing t=beta boundary?
                      segAfterRecon_i->beg->t :         //   n
                      segAfterRecon_i->beg->t + tMax;   //   y

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
      cout << settings::cout::enterRed << "Worm::antiReconnect: tReconMin=" << tReconMin << " vs tMin_prev=" << this->tMin_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tReconMax_i != tMaxi_prev) {
      cout << settings::cout::enterRed << "Worm::antiReconnect: tReconMax_i=" << tReconMax_i << " vs tMaxi_prev=" << tMaxi_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tReconMax_j != tMaxj_prev){
      cout << settings::cout::enterRed << "Worm::antiReconnect: tReconMax_j=" << tReconMax_j << " vs tMaxj_prev=" << tMaxj_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tReconMax != tMax_prev) {
      cout << settings::cout::enterRed << "Worm::antiReconnect: tReconMax=" << tReconMax << " vs tMax_prev=" << tMax_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (tRecon != this->t_prev) {
      cout << settings::cout::enterRed << "Worm::antiJump: tRecon=" << tRecon << " vs tRecon_prev=" << this->t_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// the segment length being modified by this update
  ////
  const long modifiedLength = tRecon - headSegment->end->t;

  ////
  //// parameters for the tRecon distribution
  ////
  const long intLength = tReconMax - tReconMin;
  const double lambda = localWeightRatio_exponent,
               intervalLength = intLength * this->int2time(beta);

  // TEST: make sure that the distribution parameters agree
  if (debug && test) {
    if (intLength != this->intLength_prev) {
      cout << settings::cout::enterRed << "Worm::antiReconnect: intLength=" << intLength << " vs intLength_prev=" << this->intLength_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
    if (lambda != this->lambda_prev) {
      cout << settings::cout::enterRed << "Worm::antiReconnect: lambda=" << lambda << " vs lambda_prev=" << this->lambda_prev << settings::cout::resetStyle << endl;
      if (shutItDown) this->shutDown();
    }
  }

  // the inverse probability (not exponential part) of having proposed this jump time
  double PtRecon_base;
  if (expDistrEnabled) {
    PtRecon_base = abs(lambda) > pow(10., -10.) ?
                   (1 - exp(-lambda * intervalLength)) / lambda :
                   intervalLength;
  } else {
    PtRecon_base = intervalLength;
  }

  ////
  //// the weight ratio of picking this update procedure as compared to the reverse update procedure
  ////
  const double Pweight = (double) this->Wreconnect / (double) this->WantiReconnect;

  ////
  //// the proposal distribution ratio of the proposed update: g(recon -> not recon) / g(not recon -> recon)
  //// OBS: cancellation of exponential dependence
  ////
  const double P_base = Pweight
                      / Psite
                      / PtRecon_base;


  ////
  //// other weight contributions
  ////
  double exp_U_nn = 0;

  // nearest neighbor interaction
  if (has_U_nn) {
    if (allowMultiCompWorm) {
      cout << "Worm::antiReconnect: TODO: implement me!" << endl;
      this->shutDown();
    } else {
      exp_U_nn += this->computeNNinteractionDiff(i,
                                                 segBeforeRecon_i->pop,
                                                 segAfterRecon_i->pop,
                                                 j,
                                                 segBeforeRecon_j->pop,
                                                 segAfterRecon_j->pop,
                                                 headSegment->end->t,
                                                 tRecon);
    }
  }


  ////
  //// The acceptance ratio
  //// OBS: cancellation of exponential dependence
  ////
  ////                             ╭      W(not recon)    g(not recon -> recon)  ╮
  //// A(recon -> not recon) = min │ 1, ---------------  ----------------------- │
  ////                             ╰        W(recon)      g(recon -> not recon)  ╯
  ////
  const double R = P_base / abs(localWeightRatio_base)
                 * exp(   (expDistrEnabled ? 0 : localWeightRatio_exponent * this->int2time(beta) * modifiedLength)
                        + exp_U_nn * this->int2time(beta) ),
               A = min(1., R);

  // DEBUG: the acceptance ration must always be larger than zero
  if (debug && debugFrom <= this->currUpdateCount && A < 0) {
    cout << settings::cout::enterRed << "Worm::antiReconnect: invalid value of the acceptance ration A=" << A << settings::cout::resetStyle << endl;
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
  this->numKinks--;
  if (allowMultiCompWorm) {
    for (unsigned a = 0; a < numComps; a++) {
      // jumps
      this->numJumps[a] -= !! this->wormPop[a];   // only 0 or 1 (hole anti-reconnect removes a particle jump)

      // flow
      transform(segAfterRecon_i->beg->dist.begin(),
                segAfterRecon_i->beg->dist.end(),
                this->flow.begin() + (j * numComps + a) * this->numDims,
                this->flow.begin() + (j * numComps + a) * this->numDims,
                [&] (const double & b, const double & c) {
                  return c + this->wormPop[a] * b;
                });
      transform(segAfterRecon_i->beg->dist.begin(),
                segAfterRecon_i->beg->dist.end(),
                this->flow.begin() + (i * numComps + a) * this->numDims,
                this->flow.begin() + (i * numComps + a) * this->numDims,
                [&] (const double & b, const double & c) {
                  return c + this->wormPop[a] * b;
                });
      // particle number per site
      this->numParticlesAtSite[j * numComps + a] -= modifiedLength * this->wormPop[a];
      this->numParticlesAtSite[i * numComps + a] += modifiedLength * this->wormPop[a];

      // update the particle number square
      for (unsigned b = a; b < numComps; b++) {
        this->numParticlesSquared[Worm::i_ab(a, b)] += modifiedLength
                                                     * (  this->wormPop[b] * ((int) headSegment->end->outg->pop[a] - (int) segAfterRecon_i->beg->conn->inco->pop[a])
                                                        + this->wormPop[a] * ((int) headSegment->end->outg->pop[b] - (int) segAfterRecon_i->beg->conn->inco->pop[b])
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

    cout << "Worm::anitReconnect: implement \"flow\" for multicomponent worms" << endl;
  } else {
    // jumps
    this->numJumps[actComp] -= !! actPop;   // only 0 or 1 (hole anti-jump creates a particle jump)

    // particle number per site
    this->numParticlesAtSite[j * numComps + actComp] -= modifiedLength * actPop;
    this->numParticlesAtSite[i * numComps + actComp] += modifiedLength * actPop;

    // update the particle number square
    for (unsigned b = 0; b < numComps; b++) {
      const unsigned i_ab = Worm::i_ab(actComp, b);

      if (b == actComp) {
        this->numParticlesSquared[i_ab] += modifiedLength
                                         * 2 * actPop
                                         * (  (int) headSegment->end->outg->pop[actComp]
                                            - (int) segAfterRecon_i->beg->conn->inco->pop[actComp]
                                            + actPop );

      } else {
        this->numParticlesSquared[i_ab] += modifiedLength
                                         * actPop
                                         * (  (int) headSegment->end->outg->pop[b]
                                            - (int) segAfterRecon_i->beg->conn->inco->pop[b] );
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
    transform(segAfterRecon_i->beg->dist.begin(),
              segAfterRecon_i->beg->dist.end(),
              this->flow.begin() + (j * numComps + actComp) * this->numDims,
              this->flow.begin() + (j * numComps + actComp) * this->numDims,
              [&] (const double & a, const double & b) {
                return actPop > 0 ? b + a : b - a;
              });
    transform(segAfterRecon_i->beg->dist.begin(),
              segAfterRecon_i->beg->dist.end(),
              this->flow.begin() + (i * numComps + actComp) * this->numDims,
              this->flow.begin() + (i * numComps + actComp) * this->numDims,
              [&] (const double & a, const double & b) {
                return actPop > 0 ? b + a : b - a;
              });
  }

  // update nearest neighbor interaction
  if (has_U_nn)  this->U_nn -= exp_U_nn / tMax;


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
  if (headSegment->end->outg->end) {
    // reconnection occurs on the same side of t=beta boundary

    // shift time
    segBeforeRecon_j->end->t = headSegment->end->t;

    // unlink
    segBeforeRecon_i->end->conn = nullptr;
    segBeforeRecon_j->end->conn = nullptr;

    // if this is true, the tailSegment will be removed and thus
    // the segment to become the new tail must be assigned to tailSegment
    if (headSegment == tailSegment) {
      this->tailSegments[wormIndex] = segAfterRecon_i;
    }

    // remove
    this->removeSegment(i, this->headSegmentIndices[wormIndex] + 1);
    this->removeSegment(i, this->headSegmentIndices[wormIndex] + 0);

    // update head
    this->headSegments[wormIndex] = segBeforeRecon_j;
    this->findSegmentIndex(j, this->headSegments[wormIndex]->end->t, this->headSegmentIndices[wormIndex]);
  } else {
    // reconnection occurs on the other side of t=beta boundary

    // split
    auto newHeadSegment = this->splitSegment(j, this->sites[j].size() - 1, headSegment->end->t);

    // unlink
    segAfterRecon_i->beg->conn = nullptr;
    segAfterRecon_j->beg->conn = nullptr;

    // remove
    this->removeSegment(i, this->headSegmentIndices[wormIndex] + 1);
    this->removeSegment(i, 0);
    this->removeSegment(j, 0);

    // fix population imbalance
    if (allowMultiCompWorm) {
      this->subtract(newHeadSegment->end->outg->pop, this->wormPop);
    } else {
      newHeadSegment->end->outg->pop[actComp] -= actPop;
    }

    // update head
    this->headSegments[wormIndex] = newHeadSegment;
    this->headSegmentIndices[wormIndex] = this->sites[j].size() - 2;
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
  if (verbose) cout << settings::cout::enterGreen << "...anti-reconnected" << settings::cout::resetStyle << endl;

  ////
  //// DEBUG: test a jump
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
    this->t_prev = tRecon;

    // parameters for the exponential distribution
    this->intLength_prev = intLength;
    this->lambda_prev = lambda;

    // closely related to acceptance ratio
    this->modifiedLength_prev = modifiedLength;
    this->R_prev = R;

    // try going back with the anti-update
    this->reconnect(beta, wormIndex, P_insert, true);
  }
}