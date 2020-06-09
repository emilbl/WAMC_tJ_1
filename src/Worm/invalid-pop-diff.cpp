#include "Worm.h"

using namespace std;


#pragma message("move to Hamiltonian")
bool Worm::hasInvalidPopDiff (
  const array<unsigned, numComps> & nextPop,
  const array<unsigned, numComps> & currPop
) {

  ////
  //// tJ
  ////
  if (modelType == tJ) {
    const int dPop0 = nextPop[0] - currPop[0],
              dPop1 = nextPop[1] - currPop[1];

    ////
    //// Allowed: ±[1, 0], ±[0, 1], ±[1, -1]
    ////
    return dPop0 == dPop1 || abs(dPop0) > 1 || abs(dPop1) > 1;
  }


  ////
  //// unknown model
  ////
  cout << settings::cout::enterRed
       << "Worm::hasInvalidPopDiff: ERROR: unknown model \"" << modelName << "\"  ->  EXIT"<< endl
       << settings::cout::resetStyle;
  exit(EXIT_SUCCESS);
}