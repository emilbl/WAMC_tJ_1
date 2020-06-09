#include "Hamiltonian.h"

using namespace std;


double Hamiltonian::inter (
  const array<unsigned, numComps> & pop
) const {
  double U_os = 0;

  for (unsigned a = 0; a < numComps; a++) {
    for (unsigned b = a; b < numComps; b++) {
      U_os += this->U_ab[a][b] * pop[a] * pop[b];
    }
  }

  return U_os;
}


////
////
////
double Hamiltonian::interDiff (
  const array<unsigned, numComps> & segPop,    // site population before update
  const array<int, numComps> &      wormPop,   // worm population
  const int                         sign       // add worm pop (+1) or remove worm pop (-1)
) const {
  double U_os = 0;

  for (unsigned a = 0; a < numComps; a++) {
    for (unsigned b = a; b < numComps; b++) {
      U_os += this->U_ab[a][b]
            * (  sign * (int) segPop[a] * wormPop[b]
               + sign * (int) segPop[b] * wormPop[a]
               + wormPop[a] * wormPop[b] );
    }
  }

  return U_os;
}

double Hamiltonian::interDiff (
  const array<unsigned, numComps> & segPop,    // site population before update
  const unsigned                    a,         // active component
  const int                         wormPop,   // worm population of active component
  const int                         sign       // add worm pop (+1) or remove worm pop (-1)
) const {
  double U_os = 0;

  for (unsigned b = 0; b < numComps; b++) {
    if (b == a) {
      U_os += this->U_ab[a][a]
            * (2 * sign * (int) segPop[a] + wormPop) * wormPop;
    } else {
      U_os += this->U_ab[a][b]
            * sign * (int) segPop[b] * wormPop;
    }
  }

  return U_os;
}


////
////
////
double Hamiltonian::interDiff (
  const array<unsigned, numComps> & segPop_i,   // proposed site population before update
  const array<unsigned, numComps> & segPop_j,   // current site population before update
  const array<int, numComps> &      wormPop     // worm population
) const {
  double U_os = 0;

  for (unsigned a = 0; a < numComps; a++) {
    for (unsigned b = a; b < numComps; b++) {
      U_os += this->U_ab[a][b]
            * (  2 * wormPop[a] * wormPop[b]
               + wormPop[a] * ((int) segPop_i[b] - (int) segPop_j[b])
               + wormPop[b] * ((int) segPop_i[a] - (int) segPop_j[a]) );
    }
  }

  return U_os;
}

double Hamiltonian::interDiff (
  const array<unsigned, numComps> & segPop_i,   // proposed site population before update
  const array<unsigned, numComps> & segPop_j,   // current site population before update
  const unsigned                    a,          // active component
  const int                         wormPop     // worm population of active component
) const {
  double U_os = 0;

  for (unsigned b = 0; b < numComps; b++) {
    if (b == a) {
      U_os += this->U_ab[a][a]
            * 2 * wormPop
            * (wormPop +  (int) segPop_i[a] - (int) segPop_j[a]);
    } else {
      U_os += this->U_ab[a][b]
            * wormPop
            * ((int) segPop_i[b] - (int) segPop_j[b]);
    }
  }

  return U_os;
}