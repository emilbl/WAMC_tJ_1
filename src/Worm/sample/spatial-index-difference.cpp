#include "../Worm.h"

using namespace std;


unsigned Worm::spatialIndexDifference (
  const unsigned I1,
  const unsigned I2
) const {
  ////
  //// I1 == I2 -> [size_1 / 2, size_2 / 2, ...]
  //// (centered, not the origin!)
  ////

  // unravel index
  const auto is1 = this->lattice.i2is(I1);
  const auto is2 = this->lattice.i2is(I2);

  unsigned I = 0;
  if (this->numDims == 1) {
    if (this->isPeriodic[0]) {
      I = positiveModulo<int, unsigned>( (int) is1[0] - (int) is2[0] + (int) this->latticeSize[0]/2, this->latticeSize[0]);
    } else {
      cout << "Worm::spatialIndexDifference: Implement one-dimensional non-periodic lattice!" << endl;
      if (shutItDown) this->shutDown();
    }
  } else if (this->numDims == 2) {
    if (this->isPeriodic[0] && this->isPeriodic[1]) {
      I = positiveModulo<int, unsigned>( (int) is1[0] - (int) is2[0] + (int) this->latticeSize[0]/2, this->latticeSize[0]) * this->latticeSize[1]
        + positiveModulo<int, unsigned>( (int) is1[1] - (int) is2[1] + (int) this->latticeSize[1]/2, this->latticeSize[1]);
    } else if (this->isPeriodic[0]) {
      I = positiveModulo<int, unsigned>( (int) is1[0] - (int) is2[0] + (int) this->latticeSize[0]/2, this->latticeSize[0]);
    } else {
      I = positiveModulo<int, unsigned>( (int) is1[1] - (int) is2[1] + (int) this->latticeSize[1]/2, this->latticeSize[1]);
    }
  } else {
    cout << "Worm::spatialIndexDifference: Implement dimensions higher than 2!" << endl;
    if (shutItDown) this->shutDown();
  }

  if (I < 0 || I >= this->numSites){
    cout << "Worm::spatialIndexDifference: I < 0 || I >= this->numSites" << endl;
    if (shutItDown) this->shutDown();
  }

  return I;
}