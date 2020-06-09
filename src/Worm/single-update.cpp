#include "Worm.h"

using namespace std;


void Worm::singleUpdate () {


  // while (5758 - 10 > this->currUpdateCount) {
  //   this->update(this->numGsCount);
  //   this->trySample(this->beta_target);
  // }
  // this->trySample(this->beta_target);
  // int i = 0;
  // do {
  //   this->update();
  // } while ( ! this->sampleDDC() && i++ < 1000000);

  this->update();
  // this->sampleDDC();

}