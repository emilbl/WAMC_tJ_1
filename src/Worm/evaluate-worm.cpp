#include "Worm.h"

using namespace std;

void Worm::evaluateWorm (
  const double & beta,
  int &          sign,
  double &       base,
  double &       exponent
) const {
  // reset
  base     = 1;
  exponent = 0;



  ////
  //// discontinuities
  ////
  for (const auto & headSeg : this->headSegments) {
    base *= this->H.disco(headSeg->pop, headSeg->end->outg->pop, beta);
  }
  for (const auto & tailSeg : this->tailSegments) {
    base *= this->H.disco(tailSeg->beg->inco->pop, tailSeg->pop, beta);
  }



  ////
  //// kinetic and exchange
  ////
  for (unsigned i = 0; i < this->numSites; i++) {
    for (auto seg : this->sites[i]) {
      if (
        seg->end &&
        seg->end->conn &&
        seg->siteIndex > seg->end->conn->outg->siteIndex   // <- no double counting
      ) {
        base *= this->H.jump(seg->siteIndex,                    // j
                             seg->pop,
                             seg->end->outg->pop,
                             seg->end->conn->inco->siteIndex,   // i
                             seg->end->conn->inco->pop,
                             seg->end->conn->outg->pop);
      }
    }
  }



  ////
  //// potential
  ////
  {
    double exp_U_mu = 0;

    for (unsigned i = 0; i < this->numSites; i++) {
      for (auto seg : this->sites[i]) {

        const double dt = (  (seg->end ? seg->end->t : tMax)
                           - (seg->beg ? seg->beg->t : 0) );

        exp_U_mu += this->H.poten(i, seg->pop) * dt;
      }
    }

    // convert to time (obs plus sign since chemical potential comes with negative sign in Hamiltonian definition)
    exponent += exp_U_mu * this->int2time(beta);
  }



  ////
  //// on-site interaction
  ////
  if (has_U) {
    double exp_U_os = 0;

    for (unsigned i = 0; i < this->numSites; i++) {
      for (auto seg : this->sites[i]) {

        const long dt = (  (seg->end ? seg->end->t : tMax)
                         - (seg->beg ? seg->beg->t : 0) );

        exp_U_os += this->H.inter(seg->pop) * dt;
      }
    }

    // convert to time (obs minus sign)
    exponent -= exp_U_os * this->int2time(beta);
  }



  ////
  //// nearest neighbor interaction
  ////
  if (has_U_nn) {
    double exp_U_nn = this->compute_U_nn();

    // convert to action (obs minus sign)
    exponent -= exp_U_nn * tMax * this->int2time(beta);
  }


  ////
  //// the fermionic sign is computed separately
  ////


  ////
  //// extract sign
  ////
  if (base < 0) sign = -1;
  else          sign = 1;


  ////
  //// remove sign from base
  ////
  if (base < 0) base = -base;
}
