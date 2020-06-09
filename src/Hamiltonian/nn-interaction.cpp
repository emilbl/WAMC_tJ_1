#include "Hamiltonian.h"

using namespace std;


double Hamiltonian::inter_nn (
  const array<unsigned, numComps> & pop1,
  const array<unsigned, numComps> & pop2
) const {

  ////
  //// tJ
  ////
  if (modelType == tJ) {
    // the only contribution will come from neighboring sites with opposite spin
    const int dPop0 = (int) pop1[0] - (int) pop2[0];

    if (abs(dPop0) == 1 && pop1[1] == 0 && pop2[1] == 0) return this->U_nn_ab[0][1];
    else                                                 return 0;
  }


  cout << "Hamiltonian::inter_NN: ERROR: no implementation found!" << endl;
  exit(EXIT_SUCCESS);
}


double Hamiltonian::interDiff_nn (
  const array<unsigned, numComps> & pop1_aft,
  const array<unsigned, numComps> & pop1_bef,
  const array<unsigned, numComps> & pop2
) const {

  ////
  //// tJ
  ////
  if (modelType == tJ) {
    // return difference
    return this->inter_nn(pop1_aft, pop2) - this->inter_nn(pop1_bef, pop2);
  }


  cout << "Hamiltonian::interDiff_NN: ERROR: no implementation found!" << endl;
  exit(EXIT_SUCCESS);
}