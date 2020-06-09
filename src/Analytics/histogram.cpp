#include "Analytics.h"

using namespace std;

template<typename T>
void Analytics::histogram (
  const DataSet<T> & dataSet,
  const typename vector<T>::size_type & size,
  const unsigned numBins,
  vector<double> & vals,
  vector<double> & counts
) {
  // reset
  vals = vector<double>(numBins);
  counts = vector<double>(numBins);

  ////
  //// find min and max
  ////
  auto itr = dataSet.beg;
  double minVal = *itr,
         maxVal = *itr;
  while (itr < dataSet.beg + size * dataSet.incr) {
    minVal = min(minVal, (double) *itr);
    maxVal = max(maxVal, (double) *itr);

    // increment iterator
    itr += dataSet.incr;
  }

  ////
  //// discretized values spectrum
  ////
  const double binSize = (maxVal - minVal) / numBins;
  for (unsigned i = 0; i < numBins; i++) {
    vals[i] = minVal + (i + 0.5) * binSize;
  }

  ////
  //// increment counters accordingly
  ////
  itr = dataSet.beg;
  while (itr < dataSet.beg + size * dataSet.incr) {
    const unsigned binNum = max((unsigned) min((unsigned) ((*itr - minVal) / binSize), numBins - 1), (unsigned) 0);
    counts[binNum]++;

    // increment iterator
    itr += dataSet.incr;
  }
}
template void Analytics::histogram<double> (
  const DataSet<double> &,
  const vector<double>::size_type &,
  const unsigned,
  vector<double> &,
  vector<double> &
);
