#ifndef SEGMENT_H
#define SEGMENT_H

#include <iostream>
#include <array>
#include <memory>

#include "../../etc/cout.h"

#include "../../settings/settings.h"

using namespace settings::model;

struct Node;

struct Segment : public std::enable_shared_from_this<Segment> {

  const unsigned siteIndex;

  std::array<unsigned, numComps> pop;

  std::shared_ptr<Node> beg,   // node at beginning
                        end;   // node at end

  Segment (unsigned,                         // lattice site index
           std::array<unsigned, numComps>,   // population
           std::shared_ptr<Node>,            // node at the beginning
           std::shared_ptr<Node>);           // node at the end

  void print ();
};

#endif