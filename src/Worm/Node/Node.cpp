#include "Node.h"

using namespace std;

Node::Node (
  const long & t,
  shared_ptr<Segment> inco,
  shared_ptr<Segment> outg,
  const unsigned numDims
) :
  t{t},
  inco{inco},
  outg{outg},
  conn{nullptr},
  dist(numDims)
{ }

void Node::print () {
  cout << "Node: " << this << endl
       << "• t=" << this->t << endl
       << "• inco=" << this->inco << endl
       << "• outg=" << this->outg << endl
       << "• conn=" << this->conn << endl;
}