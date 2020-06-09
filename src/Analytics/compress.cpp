#include "Analytics.h"

using namespace std;

template<typename T>
Analytics::Compressed Analytics::compress (
  const DataSet<T> & dataSet,                   // the data
  const typename vector<T>::size_type & size,   // size of list, padding excluded
  const string & quantityName                   // the name of quantity being compressed which might be used in dialog outputs
) {
  ////
  //// calculate mean
  ////
  double mean = Analytics::statAvg(dataSet, 0, size);

  ////
  //// calculate symmetric std
  ////
  const auto squaredData = Analytics::elemProd<T, T, T>(dataSet, dataSet, size);
  DataSet<T> squaredDataSet = {squaredData.begin(), 1};
  const auto squaredMean = Analytics::statAvg(squaredDataSet, 0, size);
  const auto std = sqrt(squaredMean - mean * mean);

  ////
  //// calculate correlation length
  ////
  const auto tau = Analytics::correlationLength<T>(dataSet,
                                                   min(size, (decltype(size)) 10000),   // cutoff: no reason waste any more computations
                                                   exp(-1));

  ////
  //// calculate mean std via bootstrap
  ////
  double mean_std = -1;

  auto func = [](
    const vector<Analytics::DataSet<T> > & dataSets,
    const typename vector<T>::size_type & beg,
    const typename vector<T>::size_type & end
  ) {
    return Analytics::statAvg(dataSets[0], beg, end);
  };

  mean_std = Analytics::bootstrap<T>(func,
                                     {dataSet},
                                     size,
                                     mean,
                                     tau,
                                     quantityName);




  ////
  //// return
  ////
  Analytics::Compressed compressed{mean,
                                   std,
                                   mean_std,
                                   (unsigned) tau};

  return compressed;
}
template Analytics::Compressed Analytics::compress<double> (
  const DataSet<double> & dataSet,
  const typename vector<double>::size_type &,
  const string &
);
template Analytics::Compressed Analytics::compress<int> (
  const DataSet<int> & dataSet,
  const typename vector<int>::size_type &,
  const string &
);
// template Analytics::Compressed Analytics::compress<unsigned> (
//   const DataSet<unsigned> & dataSet,
//   const typename vector<unsigned>::size_type &,
//   const string &
// );
