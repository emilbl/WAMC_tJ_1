#include "Analytics.h"

using namespace std;

template<typename T>
Analytics::Compressed Analytics::compress_signed (
  const DataSet<T> &                    dataSet,       // the data
  const DataSet<T> &                    signs,         // the signs
  const typename vector<T>::size_type & size,          // size of list, padding excluded
  const string &                        quantityName   // the name of quantity being compressed which might be used in dialog outputs
) {
  ////
  //// multiply data with sign and store in new data set
  ////
  const auto signedData = Analytics::elemProd<T, T, T>(dataSet, signs, size);
  DataSet<T> signedDataSet = {signedData.begin(), 1};

  ////
  //// calculate mean
  ////
  const auto meanSign = Analytics::statAvg(signs, 0, size);
  const auto mean     = Analytics::statAvg(signedDataSet, 0, size) / meanSign;

  ////
  //// calculate symmetric std
  ////
  const auto signedSquaredData = Analytics::elemProd<T, T, T>(dataSet, signedDataSet, size);
  DataSet<T> signedSquaredDataSet = {signedSquaredData.begin(), 1};
  const auto squaredMean = Analytics::statAvg<T>(signedSquaredDataSet, 0, size) / meanSign;
  const auto std = sqrt(squaredMean - mean * mean);


  ////
  //// calculate correlation length
  ////
  const auto tau = Analytics::correlationLength<T>(dataSet,
                                                   min(size, (decltype(size)) 10000),   // cutoff: no reason to waste any more computations
                                                   exp(-1));

  ////
  //// calculate mean std via bootstrap
  ////
  double mean_std = -1;

  auto func = [] (
    const vector<Analytics::DataSet<T> > & dataSets,
    const typename vector<T>::size_type & beg,
    const typename vector<T>::size_type & end
  ) {
    const auto nominator   = Analytics::statAvg(dataSets[0], beg, end),
               denominator = Analytics::statAvg(dataSets[1], beg, end);
    return nominator / denominator;
  };

  mean_std = Analytics::bootstrap<T>(func,
                                     {signedDataSet, signs},
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
template Analytics::Compressed Analytics::compress_signed<double> (
  const DataSet<double> &,
  const DataSet<double> &,
  const typename vector<double>::size_type &,
  const string &
);
template Analytics::Compressed Analytics::compress_signed<int> (
  const DataSet<int> &,
  const DataSet<int> &,
  const typename vector<int>::size_type &,
  const string &
);
// the quantity may not be of unsigned type!

