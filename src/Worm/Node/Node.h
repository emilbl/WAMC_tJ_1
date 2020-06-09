#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <memory>
#include <vector>

#include "../../etc/cout.h"

#include "../../settings/settings.h"

using namespace settings::model;

struct Segment;

struct Node : public std::enable_shared_from_this<Node> {

  long t;   // the imaginary time position of the node

  std::shared_ptr<Segment> inco,   // incoming segment
                           outg;   // outgoing segment

  std::shared_ptr<Node> conn;   // connecting node

  // the distance to the connecting node
  std::vector<double> dist;


  Node (const long &,               // time
        std::shared_ptr<Segment>,   // incoming segment
        std::shared_ptr<Segment>,   // outgoing segment
        const unsigned);            // number of dimensions in order to prepare the dist vector

  void print ();
};

#endif