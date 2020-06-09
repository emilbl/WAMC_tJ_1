#include "Segment.h"

using namespace std;

Segment::Segment (
  unsigned siteIndex,
  array<unsigned, numComps> pop,
  shared_ptr<Node> beg,
  shared_ptr<Node> end
) : siteIndex{siteIndex}
{
  this->pop = pop;
  this->beg = beg;
  this->end = end;
}

void Segment::print () {
  cout << "Segment: " << this << endl
       << "• siteIndex=" << this->siteIndex << endl
       << "• pop=" << this->pop << endl
       << "• beg=" << this->beg << endl
       << "• end=" << this->end << endl;
}