#include "Worm.h"

using namespace std;

shared_ptr<Segment> Worm::splitSegment (
  const unsigned      siteIndex,
  const unsigned      segmentIndex,
  const unsigned long splittingTime
) {

  //
  // ╭───╮
  // ┤   ├───────────────────────────────────
  // ╰───╯
  //
  //                     ↓
  //
  // ╭───╮       new   ╔═══╗
  // ┤   ╞═════════════╣   ╟─────────────────
  // ╰───╯             ╚═══╝
  //


  auto segment = this->sites[siteIndex][segmentIndex];


  auto newSegment = make_shared<Segment>(siteIndex,
                                         segment->pop,
                                         segment->beg,
                                         nullptr);

  auto newNode = make_shared<Node>(splittingTime,
                                   newSegment,
                                   segment,
                                   this->numDims);

  // link
  newSegment->end = newNode;

  // relink
  if (segment->beg) segment->beg->outg = newSegment;
  segment->beg = newNode;

  // insert the new segments into lattice vector
  this->sites[siteIndex].insert(this->sites[siteIndex].begin() + segmentIndex, newSegment);

  // update head segment index accordingly
  for (unsigned w = 0; w < this->numWorms; w++) {
    if (this->headSegments[w]->siteIndex == siteIndex && this->headSegmentIndices[w] >= segmentIndex) {
      this->headSegmentIndices[w]++;
    }
    if (this->tailSegments[w]->siteIndex == siteIndex && this->tailSegmentIndices[w] > segmentIndex) {
      this->tailSegmentIndices[w]++;
    }
  }

  // update tail
  for (unsigned w = 0; w < this->numWorms; w++) {
    if (this->tailSegments[w] == segment) {
      this->tailSegments[w] = newSegment;
      break;
    }
  }

  return newSegment;
}
