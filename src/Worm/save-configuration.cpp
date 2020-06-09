#include "Worm.h"

using namespace std;


void Worm::saveConfiguration (
  const double & beta
) const {
  // timer
  Tic timer{};

  // container
  nlohmann::json config;


  ////
  //// relevant parameters
  ////
  config["lattice"]            = this->lattice.name;
  config["tMax"]               = tMax;
  config["allowMultiCompWorm"] = allowMultiCompWorm;
  config["numWorms"]           = this->numWorms;
  config["beta"]               = beta;


  ////
  //// the worm configuration
  ////
  config["sites"] = nlohmann::json::array();
  for (unsigned i = 0; i < this->numSites; i++) {
    const auto & site = this->sites[i];
    config["sites"].push_back(nlohmann::json::array());
    for (const auto & segment : site) {
      nlohmann::json seg;

      // population
      seg["pop"] = segment->pop;

      // end time, if not wrapping around t=beta
      if (segment->end) seg["end"] = segment->end->t;

      // outgoing connection, if any
      if (segment->end && segment->end->conn) seg["conn"] = segment->end->conn->outg->siteIndex;

      // if head and/or tail
      if (this->numWorms) {
        const auto headItr = find(this->headSegments.begin(), this->headSegments.end(), segment);
        const auto tailItr = find(this->tailSegments.begin(), this->tailSegments.end(), segment);

        if (headItr != this->headSegments.end()) {
          const unsigned w = distance(this->headSegments.begin(), headItr);
          seg["isHead"] = w;

          // worm population
          if (allowMultiCompWorm) {
            cout << "Worm::saveConfiguration: implement me!" << endl;
          } else {
            seg["actComp"] = this->actComps[w];
            seg["actPop"]  = this->actPops[w];
          }
        }
        if (tailItr != this->tailSegments.end()) {
          const unsigned w = distance(this->tailSegments.begin(), tailItr);
          seg["isTail"] = w;
        }
      }

      // append
      config["sites"].back().push_back(seg);
    }
  }


  ////
  //// write prettified json object to file
  ////
  stringstream ss;
  ss << setw(2) << config << endl;
  this->writeToFile.aString(ss.str(), "saved-configuration.json");


  // output notification to user
  printf("[Configuration saved in %.2fs]\n", timer.toc());
  fflush(stdout);
}