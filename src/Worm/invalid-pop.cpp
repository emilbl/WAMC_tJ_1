#include "Worm.h"

using namespace std;



#pragma message("move to Hamiltonian")
bool Worm::hasInvalidPop (
  const array<unsigned, numComps> & pop
) {

  ////
  //// tJ
  ////
  if (modelType == tJ) {
    if (
      pop[0] == UINT_MAX || pop[1] == UINT_MAX || // pop[0] < 0 || pop[1] < 0 ||
      pop[0] > 1 || pop[1] > 1 || (pop[0] == 1 && pop[1] == 1)
    ) {
      return true;
    } else {
      return false;
    }
  }

  ////
  //// unknown model
  ////
  cout << settings::cout::enterRed
       << "Worm::hasInvalidPop: ERROR: unknown model \"" << modelName << "\"  ->  EXIT"<< endl
       << settings::cout::resetStyle;
  exit(EXIT_SUCCESS);
}


bool Worm::wouldHaveInvalidPop (
  const array<unsigned, numComps> & segmentPop,
  const unsigned                    actComp,
  const int                         actPop,
  const int                         sumOrDiff
) {
  ////
  //// tJ
  ////
  if (modelType == tJ) {
    const int pop0 = segmentPop[actComp] + sumOrDiff * actPop,
              pop1 = segmentPop[1 - actComp];

    if (pop0 < 0 || pop0 > 1 || (pop0 == 1 && pop1 == 1)) {
      return true;
    } else {
      return false;
    }
  }


  ////
  //// unknown model
  ////
  cout << settings::cout::enterRed
       << "Worm::wouldHaveInvalidPop: ERROR: unknown model \"" << modelName << "\"  ->  EXIT"<< endl
       << settings::cout::resetStyle;
  exit(EXIT_SUCCESS);
}


bool Worm::wouldHaveInvalidPop (
  const array<unsigned, numComps> & segmentPop,
  const array<int, numComps> &      wormPop,
  const int                         sumOrDiff
) {

  ////
  //// tJ
  ////
  if (modelType == tJ) {
    const int pop0 = segmentPop[0] + sumOrDiff * wormPop[0],
              pop1 = segmentPop[1] + sumOrDiff * wormPop[1];

    if (
      pop0 < 0 || pop1 < 0 ||
      pop0 > 1 || pop1 > 1 || (pop0 == 1 && pop1 == 1)
    ) {
      return true;
    } else {
      return false;
    }
  }


  ////
  //// unknown model
  ////
  cout << settings::cout::enterRed
       << "Worm::wouldHaveInvalidPop: ERROR: unknown model \"" << modelName << "\"  ->  EXIT"<< endl
       << settings::cout::resetStyle;
  exit(EXIT_SUCCESS);
}