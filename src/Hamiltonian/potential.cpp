#include "Hamiltonian.h"

using namespace std;


double Hamiltonian::poten (
  const unsigned i,
  const array<unsigned, numComps> & pop
) const {
  double U_mu = 0;

  for (unsigned a = 0; a < numComps; a++) {
    U_mu += pop[a] * this->V_ia[this->i_ia(i, a)];
  }

  return U_mu;
}


////
////
////
void Hamiltonian::potenDiff (
  const unsigned i,                       // site index
  const array<int, numComps> & wormPop,   // worm population
  const int sign,                         // add (+1) or subtract (-1) worm pop
  double & base,                          // base
  double & exponent                       // exponent
) const {
  for (unsigned a = 0; a < numComps; a++) {
    exponent -= sign * wormPop[a] * this->V_ia[this->i_ia(i, a)];
  }
}


void Hamiltonian::potenDiff (
  const unsigned siteIndex,   // site index
  const unsigned a,           // active component
  const int wormPop,          // worm population of active component
  const int sign,             // add (+1) or subtract (-1) worm pop
  double & base,              // base
  double & exponent           // exponent
) const {
  exponent -= sign * wormPop * this->V_ia[this->i_ia(siteIndex, a)];
}


////
////
////
void Hamiltonian::potenDiff (
  const unsigned i,                       // proposed site index
  const unsigned j,                       // current site index
  const array<int, numComps> & wormPop,   // worm population
  double & base,                          // base
  double & exponent                       // exponent
) const {
  for (unsigned a = 0; a < numComps; a++) {
   exponent -= wormPop[a] * (this->V_ia[this->i_ia(i, a)] - this->V_ia[this->i_ia(j, a)]);
  }
}

void Hamiltonian::potenDiff (
  const unsigned i,    // proposed site index
  const unsigned j,    // current site index
  const unsigned a,    // active component
  const int wormPop,   // worm population of active component
  double & base,       // base
  double & exponent    // exponent
) const {
  exponent -= wormPop * (this->V_ia[this->i_ia(i, a)] - this->V_ia[this->i_ia(j, a)]);
}