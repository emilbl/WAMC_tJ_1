#include "Worm.h"

using namespace std;

void Worm::calculateAverageSegmentLengths (
  const double & beta,
  vector<vector<double> > & container
) const {
  // reset
  container = vector<vector<double> >(2 * numComps);

  for (unsigned a = 0; a < numComps; a++) {

    // segment length container
    vector<double> segLengs;
    segLengs.reserve(10000);

    for (unsigned i = 0; i < this->numSites; i++) {
      // will hold the times
      long tFirst = 0,
           tStart = 0;

      for (unsigned s = 0; s < this->sites[i].size() - 1; s++) {
        // temporarily store the segment
        const auto seg = this->sites[i][s];

        // look for a population difference
        if (seg->end->outg->pop[a] != seg->pop[a]) {
          // calculate the difference in time in case of the segment being populated
          if (tStart > 0 && seg->pop[a] > 0) {
            segLengs.push_back(Worm::int2time(beta) * (seg->end->t - tStart));
          }

          // store the time of the first found population imbalance
          if ( ! tFirst) tFirst = seg->end->t;

          // increment start time
          tStart = seg->end->t;
        }
      }

      // add the segment crossing t=beta if populated
      if (this->sites[i][0]->pop[a] > 0) {
        segLengs.push_back(Worm::int2time(beta) * (tFirst + tMax - tStart));
      }
    }
    vector<double> vals, counts;
    Analytics::histogram<double>({segLengs.begin(), 1},
                                 segLengs.size(),
                                 this->numSites,
                                 vals,
                                 counts);

    container[2 * a] = vals;
    container[2 * a + 1] = counts;
  }
}
