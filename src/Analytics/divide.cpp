#include "Analytics.h"

using namespace std;


////
////
////
template<typename Ti, typename To>
bool Analytics::divide (
  const function<To(const vector<DataSet<Ti> > &,
                    const typename vector<Ti>::size_type &,
                    const typename vector<Ti>::size_type &)> & func,
  const vector<DataSet<Ti> > &                                 dataSets,
  const typename vector<Ti>::size_type &                       size,
  const typename vector<To>::size_type &                       minNumSegs,
  const typename vector<To>::size_type &                       minSegSize,
  vector<To> & output
) {
  bool succDiv = true;

  // work out the numbers of segments and their sizes
  const typename vector<To>::size_type segSize = size / minNumSegs,               // size of each subvector
                                       remain  = size % minNumSegs,               // remainder after integer division
                                       numBins = segSize ? minNumSegs : remain;   // actual number of bins

  // prepare
  output.clear();
  output.reserve(numBins);

  // to much correlation between adjacent segments
  if (minSegSize > segSize) {
    cout << "Analytics::divide: minSegSize = " << minSegSize << " > " << segSize << " = segSize" << endl;
    succDiv = false;
  }

  // the data set was to short
  if (minNumSegs > numBins) {
    cout << "Analytics::divide: minNumSegs = " << minNumSegs << " > " << numBins << " = numBins" << endl;
    succDiv = false;
  }


  // evaluate the function on the segments
  bool nanOrInfEncountered = false;
  typename vector<To>::size_type beg = 0,
                                 end = 0;
  for (unsigned i = 0; i < numBins; i++) {
    end = beg + (numBins - remain <= i ? segSize + 1 : segSize);

    const auto val = func(dataSets, beg, end);;
    if ( ! isnan(val) && ! isinf(val)) {
      output.emplace_back(val);
    } else if ( ! nanOrInfEncountered) {
      nanOrInfEncountered = true;
    }

    beg = end;
  }

  if (nanOrInfEncountered) {
    cout << "Analytics::divide: \"nan\" or \"inf\" encountered in " <<  numBins - output.size() << " out of " << numBins << " segments" << endl;
    succDiv = false;
  }

  // return whether the division was successful or not
  return succDiv;
}
template bool Analytics::divide<double, double> (
  const function<double(const vector<DataSet<double> > &,
                        const vector<double>::size_type &,
                        const vector<double>::size_type &)> &,
  const vector<DataSet<double> > &,
  const vector<double>::size_type &,
  const vector<double>::size_type &,
  const vector<double>::size_type &,
  vector<double> &
);
template bool Analytics::divide<int, double> (
  const function<double(const vector<DataSet<int> > &,
                        const vector<int>::size_type &,
                        const vector<int>::size_type &)> &,
  const vector<DataSet<int> > &,
  const vector<int>::size_type &,
  const vector<double>::size_type &,
  const vector<double>::size_type &,
  vector<double> &
);
template bool Analytics::divide<unsigned, double> (
  const function<double(const vector<DataSet<unsigned> > &,
                        const vector<unsigned>::size_type &,
                        const vector<unsigned>::size_type &)> &,
  const vector<DataSet<unsigned> > &,
  const vector<unsigned>::size_type &,
  const vector<double>::size_type &,
  const vector<double>::size_type &,
  vector<double> &
);