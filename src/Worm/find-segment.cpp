#include "Worm.h"

using namespace std;

void Worm::findSegmentIndex (
  const unsigned i,
  const long     t,
  unsigned &     segmentIndex
) const {
  // to be the segment on site "i" which contains the time "t"
  // and the corresponding location i that lattice vector

  unsigned lo = 0,
           up = this->sites[i].size() - 1;

  shared_ptr<Segment> curSeg;

  while (true) {
    unsigned j = (up + lo) / 2;
    curSeg = this->sites[i][j];

    if (
      ( ! curSeg->end || curSeg->end->t >= t) &&   // segment end time larger than or equal to t
      ( ! curSeg->beg || curSeg->beg->t < t)       // segment beg time smaller than t
    ) {
      // is found
      segmentIndex = j;
      return;
    } else if (debug && debugFrom <= this->currUpdateCount && lo == up) {
      // is not found
      cout << settings::cout::enterRed
           << "Worm::findSegmentIndex ERROR (0): no segment containing t=" << t << " in site " << i <<  " was found  ->  EXIT" << endl
           << settings::cout::resetStyle;
      this->shutDown();
    }

    // reduce interval
    if ( ! curSeg->end || curSeg->end->t > t) {
      // is in the lower half of the interval
      up = j - 1;
    } else {
      // is in the upper half of the interval
      lo = j + 1;
    }
  }

  // no segment was found
  cout << settings::cout::enterRed
       << "Worm::findSegmentIndex ERROR (1): no segment containing t=" << t << " in site " << i <<  " was found  ->  EXIT" << endl
       << settings::cout::resetStyle;
  this->shutDown();

  return;
}



void Worm::findSegment (
  const unsigned        i,
  const long            t,
  shared_ptr<Segment> & segment,       // to be the segment on site "i" which contains the time "t"
  unsigned &            segmentIndex   // and the corresponding location i that lattice vector
) const {
  this->findSegmentIndex(i, t, segmentIndex);
  segment = this->sites[i][segmentIndex];
}