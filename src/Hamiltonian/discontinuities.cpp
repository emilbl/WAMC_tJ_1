#include "Hamiltonian.h"

using namespace std;


double Hamiltonian::disco (
  array<unsigned, numComps> & popBefore,   // before in time
  array<unsigned, numComps> & popAfter,    // after in time
  const double &              beta         // the current temperature
) const {

  double popContr = 1,
         etaContr = 1;

  for (unsigned compIndex = 0; compIndex < numComps; compIndex++) {
    if (popBefore[compIndex] > popAfter[compIndex]) {
      // worm head
      popContr *= popBefore[compIndex];
      etaContr *= this->eta_a[compIndex];
    } else if (popBefore[compIndex] < popAfter[compIndex]) {
      // worm tail
      popContr *= popAfter[compIndex];
      etaContr *= this->eta_a[compIndex];
    }
  }

  return etaContr * sqrt(popContr / (beta * this->numSites));
}


////
////
////
void Hamiltonian::discoDiff (
  const array<unsigned, numComps> & segPop,   // site population before update
  const array<int, numComps> & wormPop,       // worm population
  const int sign,                             // add (+1) or subtract (-1) worm pop
  const double & beta,                        // the current temperature
  double & base,                              // base
  double & exponent                           // exponent
) const {
  double popContr = 1,
         etaContr = 1;

  for (unsigned a = 0; a < numComps; a++) {
    if (sign * wormPop[a] < 0) {
      // hole worm
      popContr *= segPop[a];
      etaContr *= this->eta_a[a];
    } else if (sign * wormPop[a] > 0) {
      // particle worm
      popContr *= segPop[a] + 1;
      etaContr *= this->eta_a[a];
    }
  }

  base *= pow(etaContr, 2)
        * popContr
        / (beta * this->numSites);
}

void Hamiltonian::discoDiff (
  const unsigned segPop,    // site population before update
  const unsigned a,         // active component
  const int      wormPop,   // worm population of active component
  const int      sign,      // add (+1) or subtract (-1) worm pop
  const double & beta,      // the current temperature
  double &       base,      // base
  double &       exponent   // exponent
) const {
  base *= pow(this->eta_a[a], 2)
        * (sign * wormPop > 0 ?
           segPop + 1 :                    // particle worm
           segPop)                         // hole worm
        / (beta * this->numSites);
}



////
////
////
void Hamiltonian::discoDiffHead (
  const array<unsigned, numComps> & newPop,    // the new population of the head segment
  const array<unsigned, numComps> & currPop,   // the current population of the head segment
  const array<int, numComps> & wormPop,        // worm population
  const int addOrSubtract,                     // add (+1) or subtract (-1) worm pop
  double & base,                               // base
  double & exponent                            // exponent
) const {

  double popRatio = 1;

  for (unsigned a = 0; a < numComps; a++) {
    if (wormPop[a] > 0) {
      // particle worm
      popRatio *= (newPop[a] + addOrSubtract * wormPop[a]) / (double) currPop[a];
    } else if (wormPop[a] < 0) {
      // hole worm
      popRatio *= (newPop[a] + 1 + addOrSubtract * wormPop[a]) / (double) (currPop[a] + 1.);
    }
  }

  base *= sqrt(popRatio);
}

void Hamiltonian::discoDiffTail (
  const array<unsigned, numComps> & newPop,    // the new population of the segment before the tail segment
  const array<unsigned, numComps> & currPop,   // the current population of the segment before the tail segment
  const array<int, numComps> & wormPop,        // worm population
  const int addOrSubtract,                     // add or subtract the worm population of the "new population" (0 do nothing)
  double & base,
  double & exponent
) const {

  double popRatio = 1;

  for (unsigned a = 0; a < numComps; a++) {
    if (wormPop[a] > 0) {
      // particle worm
      popRatio *= (newPop[a] + 1 + addOrSubtract * wormPop[a]) / (double) (currPop[a] + 1.);
    } else if (wormPop[a] < 0) {
      // hole worm
      popRatio *= (newPop[a] + addOrSubtract * wormPop[a]) / (double) currPop[a];
    }
  }

  base *= sqrt(popRatio);
}

////
////
////
void Hamiltonian::discoDiffHead (
  const unsigned newPop,          // the new population of the head segment
  const unsigned currPop,         // the current population of the head segment
  const int      wormPop,         // worm population of active component
  const int      addOrSubtract,   // add (+1) or subtract (-1) worm pop
  double &       base,            // base
  double &       exponent         // exponent
) const {
  if (wormPop > 0) {
    // particle worm
    base *= sqrt((newPop + addOrSubtract * wormPop) / (double) currPop);
  } else if (wormPop < 0) {
    // hole worm
    base *= sqrt((newPop + 1 + addOrSubtract * wormPop) / (double) (currPop + 1.));
  }
}

void Hamiltonian::discoDiffTail (
  const unsigned newPop,          // the new population of the segment before the tail segment
  const unsigned currPop,         // the current population of the segment before the tail segment
  const int      wormPop,         // worm population
  const int      addOrSubtract,   // add or subtract the worm population of the "new population" (0 do nothing)
  double &       base,
  double &       exponent
) const {
  if (wormPop > 0) {
    // particle worm
    base *= sqrt((newPop + 1 + addOrSubtract * wormPop) / (double) (currPop + 1.));
  } else if (wormPop < 0) {
    // hole worm
    base *= sqrt((newPop + addOrSubtract * wormPop) / (double) currPop);
  }
}

void Hamiltonian::discoDiffHeadOrTail (
  const unsigned segPop_l,        // current population on the left hand side of discontinuity
  const unsigned segPop_r,        // current population on the right hand side of discontinuity
  const int      wormPop,         // worm population
  const int      addOrSubtract,   // add or subtract the worm population to the current populations
  double &       base
) const {
  if (segPop_l > segPop_r) {
    // head discontinuity
    base *= sqrt((segPop_l + addOrSubtract * wormPop) / (double) segPop_l);
  } else if (segPop_l < segPop_r) {
    // tail discontinuity
    base *= sqrt((segPop_r + addOrSubtract * wormPop) / (double) segPop_r);
  }
}