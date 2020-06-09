#include "Worm.h"

using namespace std;

void Worm::removeSegment (
  const unsigned siteIndex,
  const unsigned segmentIndex
) {

  //             NOT THE LAST SEGMENT                                            LAST SEGMENT
  //
  //           ╭───╮    remove   ╔═══╗                               ╭───╮           ╔═══╗   remove
  // ┈ ────────┤   ╞═════════════╣   ╟───────── ┈        or        ┈ ┤   ├───────────╢   ╠═════════════
  //           ╰───╯             ╚═══╝                               ╰───╯           ╚═══╝
  //
  //                       ↓                                                           ↓
  //
  //           ╭───╮                                                 ╭───╮
  // ┈ ────────┤   ├─────────────────────────── ┈        or        ┈ ┤   ├─────────────────────────────
  //           ╰───╯                                                 ╰───╯
  //

  auto segment = this->sites[siteIndex][segmentIndex];

  // update head and tail quantities
  for (unsigned w = 0; w < this->numWorms; w++) {
    if (this->headSegments[w]->siteIndex == siteIndex && this->headSegmentIndices[w] > segmentIndex) {
      this->headSegmentIndices[w]--;
    }
    if (this->tailSegments[w]->siteIndex == siteIndex) {
      if (this->tailSegmentIndices[w] > segmentIndex) this->tailSegmentIndices[w]--;
      if (segment == this->tailSegments[w] && segment->end) {
        // if the segment is crossing the t=beta boundary and being removed it cannot belong
        // to a worm other than one being removed and in that case it does not matter
        this->tailSegments[w] = segment->end->outg;
      }
    }

  }

  shared_ptr<Node> node;
  if (segment->end) {
    node = segment->end;
  } else {
    node = segment->beg;
  }

  // make sure that the node to be removed is not connected to another node
  if (debug && debugFrom <= this->currUpdateCount && node->conn) {
    cout << settings::cout::enterRed
         << "Worm::removeSegment ERROR: node to be removed still connected  ->  EXIT" << endl
         << settings::cout::resetStyle;
    this->shutDown();
  }

  // relink
  if (segment->beg && segment->end) {
    segment->beg->outg = node->outg;
    node->outg->beg    = segment->beg;
  } else if (segment->end) {
    node->outg->beg = nullptr;
    node->outg      = nullptr;
  } else {
    node->inco->end = nullptr;
    node->inco      = nullptr;
  }

  // remove from vector
  this->sites[siteIndex].erase(this->sites[siteIndex].begin() + segmentIndex);

  // unlink
  if (segment->end) {
    segment->end = nullptr;
    node->inco = nullptr;
  } else {
    segment->beg = nullptr;
    node->outg = nullptr;
  }

}