#include "Worm.h"

using namespace std;

void Worm::preCheckWormWeight (
  const double & beta
) {
  // PROFILING
  // profiler("Worm::preCheckWormWeight");

  // evaluate worm weight and store the weight to be used after the update
  this->evaluateWorm(beta,
                     this->preWormWeightSign,
                     this->preWormWeightBase,
                     this->preWormWeightExponent);
}


void Worm::checkWormWeight (
  const double & beta,
  const double & localWeightRatioBase,
  const double & localWeightRatioExponent
) const {
  // PROFILING
  // profiler("Worm::checkWormWeight");

  // non-unity acceptance ratio

  ////
  //// OBS: preCheckWormWeight must have been run before!
  ////

  int postWormWeightSign;
  double postWormWeightBase,
         postWormWeightExponent;
  this->evaluateWorm(beta,
                     postWormWeightSign,
                     postWormWeightBase,
                     postWormWeightExponent);


  if (settings::mode::verbose) {
    cout << "Worm::checkWormWeight: W = " << postWormWeightBase << " * exp(" << postWormWeightExponent << ")" << endl;
  }


  double globalWeightRatioBase     = postWormWeightBase        / this->preWormWeightBase,
         globalWeightRatioExponent = postWormWeightExponent    - this->preWormWeightExponent,
         WeightRatioRatioBase      = globalWeightRatioBase     / localWeightRatioBase,
         WeightRatioRatioExponent  = globalWeightRatioExponent - localWeightRatioExponent,
         WeightRatioRatio          = WeightRatioRatioBase      * exp(WeightRatioRatioExponent);


  // check for sign
  if (localWeightRatioBase < 0) {
    cout << settings::cout::enterRed
         << "Worm::checkWormWeight: ERROR: negative worm weight ratio: " << localWeightRatioBase << endl
         << settings::cout::resetStyle;
    if (settings::mode::shutItDown) this->shutDown();
  }



  // check the sign
  if (this->sign != postWormWeightSign) {
    cout << settings::cout::enterRed
         << "Worm::checkWormWeight: ERROR: the signs do not agree" << endl
         << "current: " << this->sign << endl
         << "actual:  " << postWormWeightSign << endl
         << settings::cout::resetStyle;
    if (settings::mode::shutItDown) this->shutDown();
  }


  #pragma message("increase detailedBalanceMaxDiff for each kink since they lead to computational inaccuracies")


  // the check
  if (abs(1 - WeightRatioRatio) > detailedBalanceMaxDiff) {
    // the base should not be to small -> rounding errors
    if (postWormWeightBase < pow(10., -300) && this->preWormWeightBase < pow(10., -300)) {
      // if (settings::mode::verbose) {
        cout << settings::cout::enterRed
             << "Worm::checkWormWeight: the detailed balance maximum difference surpassed"
             << " but the global bases to small for an accurate ratio ("
             << postWormWeightBase << ", " << this->preWormWeightBase << ")"
             << settings::cout::resetStyle << endl;
      // }
      return;
    }

    cout << settings::cout::enterRed
         << "Worm::checkWormWeight: the detailed balance maximum difference surpassed"
         << " (" << abs(1 - WeightRatioRatio)
         <<  " = | 1 - " << WeightRatioRatioBase << " * exp(" << WeightRatioRatioExponent << " ) |"
         << " > " << detailedBalanceMaxDiff << ")"
         << settings::cout::resetStyle << endl;
    if (settings::mode::verbose) {
      cout << "globalWormWeightRatio = " << globalWeightRatioBase << " * exp( " << globalWeightRatioExponent << " ) = "
           << globalWeightRatioBase * exp(globalWeightRatioExponent) << endl
           << "localWormWeightRatio  = " << localWeightRatioBase << " * exp( " << localWeightRatioExponent << " ) = "
           << localWeightRatioBase * exp(localWeightRatioExponent) << endl
           << "ratio                 = " << WeightRatioRatioBase << " * exp( " << WeightRatioRatioExponent << " ) = "
           << WeightRatioRatio << endl;
    }
    if (settings::mode::shutItDown) this->shutDown();
  }


}