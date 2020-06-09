#include "Worm.h"

using namespace std;


void Worm::computeWormQuantities (
  unsigned &                             numKinks,
  double &                               U_nn,
  array<unsigned, numComps> &            numJumps,
  array<unsigned, numExternalInters> &   numExchanges,
  array<unsigned long long, numComps> &  numParticles,
  array<unsigned long long, numInters> & numParticlesSquared,
  array<vector<int>, numComps> &         numWinds,
  vector<double> &                       flow,
  vector<unsigned long long> &           numParticlesAtSite
) const {

  ////
  //// initialize / reset
  ////
  numKinks            = 0;
  U_nn                = 0;
  numJumps            = {};
  numExchanges        = {};
  numParticles        = {};
  numParticlesSquared = {};

  fill(numWinds.begin(), numWinds.end(), vector<int>(this->numDims));

  numParticlesAtSite = vector<unsigned long long>(numSites * numComps);
  flow               = vector<double>(numSites * numComps * numDims);


  ////
  //// compute nearest neighbor interaction
  ////
  if (has_U_nn) U_nn = this->compute_U_nn();


  for (unsigned i = 0; i < this->numSites; i++) {
    for (unsigned si = 0; si < this->sites[i].size(); si++) {
      const auto segment = this->sites[i][si];


      ////
      //// particle number and particle number square
      ////
      for (unsigned a = 0; a < numComps; a++) {

        unsigned long timeDiff = (segment->end ? segment->end->t : tMax)
                               - (segment->beg ? segment->beg->t : 0);

        numParticlesAtSite[i * numComps + a] += timeDiff * segment->pop[a];

        numParticles[a] += timeDiff * segment->pop[a];

        for (unsigned b = a; b < numComps; b++) {
          numParticlesSquared[Worm::i_ab(a, b)] += timeDiff * segment->pop[a] * segment->pop[b];
        }
      }

      ////
      //// jumps, exchange, winding and flow number
      ////

      // is there a jump?
      if (segment->end && segment->end->conn) {

        // calculate the population difference in the jump
        array<int, numComps> popDiff = Worm::diff<int, unsigned, unsigned>(segment->end->outg->pop, segment->pop);

        // to avoid double counting
        if (segment->siteIndex > segment->end->conn->outg->siteIndex) {

          // kinks
          numKinks++;

          // jump or exchange
          vector<unsigned> comps;
          for (unsigned a = 0; a < numComps; a++) {
            if (popDiff[a] != 0) comps.emplace_back(a);
          }
          if (comps.size() == 1) {
            // is jump
            numJumps[comps[0]]++;
          } else if (comps.size() == 2) {
            // is exchange
            numExchanges[i_ab_external(comps[0], comps[1])]++;
          } else {
            //this should never happen
            cout << "Worm::checkWormQuantities: ERROR: super exchange" << endl;
            this->shutDown();
          }
        }


        // the winding number that goes with the jump j -> i
        auto W = this->lattice.boundaryCrossings(i, segment->end->conn->outg->siteIndex);

        for (unsigned a = 0; a < numComps; a++) {
          if (popDiff[a] > 0) {
            // add the winding contribution
            Worm::add(numWinds[a], W);
          }

          // for flow we need to count all nonzero contributions
          if (popDiff[a]) {
            transform(segment->end->dist.begin(),
                      segment->end->dist.end(),
                      flow.begin() + (i * numComps + a) * this->numDims,
                      flow.begin() + (i * numComps + a) * this->numDims,
                      [&] (const double & A, const double & B) {
                        return (popDiff[a] > 0 ? B - A : B + A);
                      });
          }
        }
      }
    }
  }


}