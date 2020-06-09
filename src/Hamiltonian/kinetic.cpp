#include "Hamiltonian.h"

using namespace std;


void Hamiltonian::kinet (
  unsigned i1,                                // lattice site index 1
  unsigned i2,                                // lattice site index 2
  array<unsigned, numComps> & popBefJump_1,   // population at lattice site indexed 1 before jump
  array<unsigned, numComps> & popAftJump_1,   // population at lattice site indexed 1 after jump
  array<unsigned, numComps> & popBefJump_2,   // population at lattice site indexed 2 before jump
  array<unsigned, numComps> & popAftJump_2,   // population at lattice site indexed 2 after jump
  double & base                               // base
) const {

  double popContr = 1,
         tContr = 1;

  for (unsigned a = 0; a < numComps; a++) {
    // calculate worm pop if the worm was going from i1 to i2
    int wormPopA = popBefJump_1[a] - popAftJump_1[a];

    if (wormPopA > 0) {
      // particle jumps from i1 -> i2
      popContr *= popBefJump_1[a] * popAftJump_2[a];
      tContr *= this->t_ija[this->i_ija(i2, i1, a)];
    } else if (wormPopA < 0) {
      // particle jumps from i2 -> i1
      popContr *= popBefJump_2[a] * popAftJump_1[a];
      tContr *= this->t_ija[this->i_ija(i1, i2, a)];
    }
  }

  base *= tContr * sqrt(popContr)
        * -1;   // the minus sign (-1)^n
}



////
////
////
void Hamiltonian::kinetDiff (
  const unsigned i,                             // proposed site index
  const unsigned j,                             // current site index
  const array<unsigned, numComps> & segPop_i,   // proposed site population before update
  const array<unsigned, numComps> & segPop_j,   // current site population before update
  const array<int, numComps> & wormPop,         // worm population
  double & base                                 // base
) const {
  int popContr = 1;

  for (unsigned a = 0; a < numComps; a++) {
    if (wormPop[a] > 0) {
      // particle worm
      popContr *= (segPop_i[a] + 1) * segPop_j[a];
      base *= this->t_ija[this->i_ija(i, j, a)];
    } else if (wormPop[a] < 0) {
      // hole worm
      popContr *= segPop_i[a] * (segPop_j[a] + 1);
      base *= this->t_ija[this->i_ija(j, i, a)];
    }
  }

  base *= sqrt((double) popContr)
        * -1;   // the minus sign (-1)^n
}


void Hamiltonian::kinetDiff (
  const unsigned i,          // proposed site index
  const unsigned j,          // current site index
  const unsigned segPop_i,   // proposed site population before update
  const unsigned segPop_j,   // current site population before update
  const unsigned a,          // active component
  const int wormPop,         // worm population of active component
  double & base              // base
) const {
  if (wormPop > 0) {
    // particle worm
    base *= this->t_ija[this->i_ija(i, j, a)]
          * sqrt((segPop_i + 1) * segPop_j)
          * -1;   // the minus sign (-1)^n
  } else if (wormPop < 0) {
    // hole worm
    base *= this->t_ija[this->i_ija(j, i, a)]
          * sqrt(segPop_i * (segPop_j + 1))
          * -1;   // the minus sign (-1)^n
  }
}




////
////
////
void Hamiltonian::kinetDiff (
  const unsigned i1,                                // lattice site index 1
  const unsigned i2,                                // lattice site index 2
  const array<unsigned, numComps> & popBefJump_1,   // current population at lattice site indexed 1 before jump
  const array<unsigned, numComps> & popAftJump_1,   // current population at lattice site indexed 1 after jump
  const array<unsigned, numComps> & popBefJump_2,   // current population at lattice site indexed 2 before jump
  const array<unsigned, numComps> & popAftJump_2,   // current population at lattice site indexed 2 after jump
  const array<int, numComps> & wormPop,             // worm population
  const int addPopBefJump_1,                        // add or subtract worm to/from site indexed 1 before jump
  const int addPopAftJump_1,                        // add or subtract worm to/from site indexed 1 after jump
  const int addPopBefJump_2,                        // add or subtract worm to/from site indexed 2 before jump
  const int addPopAftJump_2,                        // add or subtract worm to/from site indexed 2 before jump
  double & base                                     // base
) const {
  cout << "Hamiltonian::kinetDiff: implement me! 11" << endl;
  exit(EXIT_SUCCESS);

  // cout << "Hamiltonian::kinetDiff: this would rather be an exchange! 1" << endl;

  // double jumpAmplitudeContr = 1,
  //        popContr = 1;

  // for (unsigned a = 0; a < numComps; a++) {
  //   ////
  //   //// current contribution
  //   ////
  //   const int popDiffBefore = popBefJump_1[a] - popAftJump_1[a];
  //   if (popDiffBefore > 0) {
  //     // particle worm from i1 -> i2
  //     popContr /= popBefJump_1[a] * popAftJump_2[a];
  //     jumpAmplitudeContr /= this->t_ija[this->i_ija(i2, i1, a)];
  //   } else if (popDiffBefore < 0) {
  //     // particle worm from i2 -> i1
  //     popContr /= popBefJump_2[a] * popAftJump_1[a];
  //     jumpAmplitudeContr /= this->t_ija[this->i_ija(i1, i2, a)];
  //   }

  //   ////
  //   //// contribution after update
  //   ////
  //   const int popDiffAfter = popDiffBefore + (addPopBefJump_1 - addPopAftJump_1) * wormPop[a];
  //   if (popDiffAfter > 0) {
  //     // particle worm from i1 -> i2
  //     popContr *= (popBefJump_1[a] + addPopBefJump_1 * wormPop[a])
  //               * (popAftJump_2[a] + addPopAftJump_2 * wormPop[a]);
  //     jumpAmplitudeContr *= this->t_ija[this->i_ija(i2, i1, a)];
  //   } else if (popDiffAfter  < 0) {
  //     // particle worm from i2 -> i1
  //     popContr *= (popBefJump_2[a] + addPopBefJump_2 * wormPop[a])
  //               * (popAftJump_1[a] + addPopAftJump_1 * wormPop[a]);
  //     jumpAmplitudeContr *= this->t_ija[this->i_ija(i1, i2, a)];
  //   }
  // }

  // base *= jumpAmplitudeContr * sqrt(popContr);
}


void Hamiltonian::kinetDiff (
  const unsigned i1,             // lattice site index 1
  const unsigned i2,             // lattice site index 2
  const unsigned popBefJump_1,   // current population at lattice site indexed 1 before jump
  const unsigned popAftJump_1,   // current population at lattice site indexed 1 after jump
  const unsigned popBefJump_2,   // current population at lattice site indexed 2 before jump
  const unsigned popAftJump_2,   // current population at lattice site indexed 2 after jump
  const unsigned a,              // active component
  const int wormPop,             // worm population of active component
  const int addPopBefJump_1,     // add or subtract worm to/from site indexed 1 before jump
  const int addPopAftJump_1,     // add or subtract worm to/from site indexed 1 after jump
  const int addPopBefJump_2,     // add or subtract worm to/from site indexed 2 before jump
  const int addPopAftJump_2,     // add or subtract worm to/from site indexed 2 before jump
  double & base                  // base
) const {

  double jumpAmplitudeContr = 1,
         popContr = 1;

  ////
  //// current contribution
  ////
  const int popDiffBefore = popBefJump_1 - popAftJump_1;
  if (popDiffBefore > 0) {
    // particle worm from i1 -> i2
    popContr /= popBefJump_1 * popAftJump_2;
    jumpAmplitudeContr /= this->t_ija[this->i_ija(i2, i1, a)];
  } else if (popDiffBefore < 0) {
    // particle worm from i2 -> i1
    popContr /= popBefJump_2 * popAftJump_1;
    jumpAmplitudeContr /= this->t_ija[this->i_ija(i1, i2, a)];
  }

  ////
  //// contribution after update
  ////
  const int popDiffAfter = popDiffBefore + (addPopBefJump_1 - addPopAftJump_1) * wormPop;
  if (popDiffAfter > 0) {
    // particle worm from i1 -> i2
    popContr *= (popBefJump_1 + addPopBefJump_1 * wormPop)
              * (popAftJump_2 + addPopAftJump_2 * wormPop);
    jumpAmplitudeContr *= this->t_ija[this->i_ija(i2, i1, a)];
  } else if (popDiffAfter  < 0) {
    // particle worm from i2 -> i1
    popContr *= (popBefJump_2 + addPopBefJump_2 * wormPop)
              * (popAftJump_1 + addPopAftJump_1 * wormPop);
    jumpAmplitudeContr *= this->t_ija[this->i_ija(i1, i2, a)];
  }

  base *= jumpAmplitudeContr * sqrt(popContr);
}

