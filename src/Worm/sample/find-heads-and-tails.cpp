#include "../Worm.h"

using namespace std;

tuple<shared_ptr<Segment>, shared_ptr<Segment>, shared_ptr<Segment>, shared_ptr<Segment> >
Worm::findHeadsAndTails (
  const unsigned actComp
) const {
  shared_ptr<Segment> hs1 = nullptr;
  shared_ptr<Segment> hs2 = nullptr;
  shared_ptr<Segment> ts1 = nullptr;
  shared_ptr<Segment> ts2 = nullptr;

  if (this->actPops[0] > 0) {
    hs1 = this->headSegments[0];

    // true head segment of second worm
    // (the tails are not necessarily in order with the heads)
    if (this->actPops[1] > 0) {
      hs2 = this->headSegments[1];

      ts1 = this->tailSegments[1];
      ts2 = this->tailSegments[0];
    } else if (this->tailSegments[0]->beg->inco->pop[actComp] > this->tailSegments[0]->pop[actComp]) {
      hs2 = this->tailSegments[0]->beg->inco;

      ts1 = this->tailSegments[1];
      ts2 = this->headSegments[1]->end->outg;
    } else {
      hs2 = this->tailSegments[1]->beg->inco;

      ts1 = this->tailSegments[0];
      ts2 = this->headSegments[1]->end->outg;
    }
  } else if (this->actPops[1] > 0) {
    hs1 = this->headSegments[1];

    // true head segment of second worm
    // (the tails are not necessarily in order with the heads)
    if (this->tailSegments[0]->beg->inco->pop[actComp] > this->tailSegments[0]->pop[actComp]) {
      hs2 = this->tailSegments[0]->beg->inco;

      ts1 = this->tailSegments[1];
      ts2 = this->headSegments[0]->end->outg;
    } else {
      hs2 = this->tailSegments[1]->beg->inco;

      ts1 = this->tailSegments[0];
      ts2 = this->headSegments[0]->end->outg;
    }
  } else {
    // both heads carry negative population
    hs1 = this->tailSegments[0]->beg->inco;

    // true head segment of second worm
    hs2 = this->tailSegments[1]->beg->inco;

    ts1 = this->headSegments[0]->end->outg;
    ts2 = this->headSegments[1]->end->outg;
  }
  return {hs1, hs2, ts1, ts2};
}