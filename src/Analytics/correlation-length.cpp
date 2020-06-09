#include "Analytics.h"

using namespace std;

template<typename T>
typename vector<T>::size_type Analytics::correlationLength (
  const DataSet<T> & dataSet,                       // the data
  const typename vector<T>::size_type & size,       // size of list, padding excluded
  const double & treshold                           // threshold
) {
  // the last element
  const auto fin = dataSet.beg + dataSet.incr * size;

  // quick check if all elements are the same
  bool identicalData = true;
  for (auto citr = dataSet.beg + dataSet.incr; citr < fin; citr += dataSet.incr) {
    if (*citr != *dataSet.beg) {
      identicalData = false;
      break;
    }
  }

  if (identicalData) {
    return size;
  }

  // useful variable
  double summed = 0;
  for (auto citr = dataSet.beg; citr < fin; citr += dataSet.incr) {
    summed += *citr;
  }

  // autocorrelation with periodic boundaries
  double crit;
  for (typename vector<T>::size_type i = 0; i < size/2 + 1; i++) {
    double sum = 0;
    for (typename vector<T>::size_type j = 0; j < size; j++) {
      sum += *(dataSet.beg + dataSet.incr * j) * *(dataSet.beg + dataSet.incr * ((i + j) % size));
    }

    // subtract the "mean" value
    sum -= summed * summed / (double) size;

    // define the critical value and check for correlation length
    if ( ! i) {
      crit = sum * treshold;
    } else if (sum < crit) {
      return i;
    }
  }

  // in case the values have not yet fallen off quickly enough
  return size;
}
template typename vector<double>::size_type Analytics::correlationLength<double> (
  const DataSet<double> &,
  const vector<double>::size_type &,
  const double &
);
template typename vector<int>::size_type Analytics::correlationLength<int> (
  const DataSet<int> &,
  const vector<int>::size_type &,
  const double &
);
template typename vector<unsigned>::size_type Analytics::correlationLength<unsigned> (
  const DataSet<unsigned> &,
  const vector<unsigned>::size_type &,
  const double &
);
