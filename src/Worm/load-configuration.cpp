#include "Worm.h"

using namespace std;


void Worm::loadConfiguration (
  const nlohmann::json & config
) {
  // timer
  Tic timer{};

  // allowMultiCompWorm must agree
  if (config["allowMultiCompWorm"] != allowMultiCompWorm) {
    cout << settings::cout::enterRed
         << "Worm::loadConfiguration: ERROR: mismatch of allowMultiCompWorm: " << config["allowMultiCompWorm"] << " vs " << allowMultiCompWorm
         << settings::cout::resetStyle
         << endl;
  }

  // the lattice type must agree
  if (config["lattice"] != this->lattice.name) {
    cout << settings::cout::enterRed
         << "Worm::loadConfiguration: ERROR: mismatch of lattice type: " << config["lattice"] << " vs " << this->lattice.name
         << settings::cout::resetStyle
         << endl;
  }

  // system sizes must agree
  if (config["sites"].size() != this->numSites) {
    cout << settings::cout::enterRed
         << "Worm::loadConfiguration: ERROR: mismatch of the number of lattice sites: " << config["sites"].size() << " vs " << this->numSites
         << settings::cout::resetStyle
         << endl;
  }

  // the tMax must agree
  if (config["tMax"] != tMax) {
    cout << settings::cout::enterRed
         << "Worm::loadConfiguration: ERROR: mismatch of \"tMax\": " << config["tMax"] << " vs " << tMax
         << settings::cout::resetStyle
         << endl;
  }


  ////
  //// prepare / clear system
  ////
  this->numWorms = config["numWorms"].get<unsigned>();
  this->headSegments.resize(this->numWorms);
  this->tailSegments.resize(this->numWorms);
  this->headSegmentIndices.resize(this->numWorms);
  this->tailSegmentIndices.resize(this->numWorms);
  if (allowMultiCompWorm) {
    cout << "Worm::loadConfiguration: implement me!" << endl;
  } else {
    this->actComps.resize(this->numWorms);
    this->actPops.resize(this->numWorms);
  }
  for (auto & site : this->sites) site.clear();


  ////
  //// load first only the segments
  ////
  for (unsigned i = 0; i < numSites; i++) {
    for (unsigned seg_i = 0; seg_i < config["sites"][i].size(); seg_i++) {

      // create new segment
      auto segment = make_shared<Segment>(i, config["sites"][i][seg_i]["pop"], nullptr, nullptr);

      // link with previous segment
      if (seg_i > 0) {
        // create node
        auto node = make_shared<Node>(config["sites"][i][seg_i - 1]["end"],
                                      this->sites[i].back(),
                                      segment,
                                      this->numDims);

        // connect segments to node
        segment->beg = node;
        this->sites[i].back()->end = node;
      }

      // if head or tail
      if (config["sites"][i][seg_i].find("isHead") != config["sites"][i][seg_i].end()) {
        const auto w = config["sites"][i][seg_i]["isHead"].get<unsigned>();
        this->headSegments[w]       = segment;
        this->headSegmentIndices[w] = seg_i;
        if (allowMultiCompWorm) {
          cout << "Worm::loadConfiguration: implement me!" << endl;
        } else {
          this->actComps[w] = config["sites"][i][seg_i]["actComp"].get<unsigned>();
          this->actPops[w]  = config["sites"][i][seg_i]["actPop"].get<int>();
        }
      }
      if (config["sites"][i][seg_i].find("isTail") != config["sites"][i][seg_i].end()) {
        const auto w = config["sites"][i][seg_i]["isTail"].get<unsigned>();
        this->tailSegments[w]       = segment;
        this->tailSegmentIndices[w] = seg_i;
      }

      // add to site
      this->sites[i].push_back(segment);
    }
  }


  ////
  //// connect nodes
  ////
  for (auto & site : this->sites) {
    for (unsigned seg_i = 0; seg_i < site.size(); seg_i++) {
      auto segment = site[seg_i];

      // check if connected to node
      if (segment->end) {
        const unsigned i = segment->siteIndex;
        const auto seg_json = config["sites"][i][seg_i];

        // look up the connection, if any
        if (seg_json.find("conn") != seg_json.end()) {
          const unsigned j = seg_json["conn"];

          // no double connections
          if (i > j) {

            // find corresponding node
            auto connectedNode = this->sites[j][0]->end;
            while (connectedNode->t != segment->end->t) connectedNode = connectedNode->outg->end;

            // connect
            segment->end->conn = connectedNode;
            connectedNode->conn = segment->end;

            // set the distance
            vector<unsigned> NNs;
            vector<vector<double> > dists;
            this->lattice.getNNsAndDists(i, NNs, dists);
            for (unsigned k = 0; k < NNs.size(); k++) {
              if (NNs[k] == j) {
                segment->end->dist = dists[k];
                for (unsigned d = 0; d < this->numDims; d++) connectedNode->dist[d] = -dists[k][d];
                break;
              }
            }
          }
        }
      }
    }
  }


  ////
  //// compute quantities
  ////
  this->computeWormQuantities(this->numKinks,
                              this->U_nn,
                              this->numJumps,
                              this->numExchanges,
                              this->numParticles,
                              this->numParticlesSquared,
                              this->numWinds,
                              this->flow,
                              this->numParticlesAtSite);


  // ensure everything is ok
  this->checkWormConfiguration();
  this->checkWormQuantities();

  // output notification to user
  printf("[Configuration loaded in %.2fs]\n", timer.toc());
  fflush(stdout);
}