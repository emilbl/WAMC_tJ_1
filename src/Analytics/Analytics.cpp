#include "Analytics.h"

using namespace std;


////
////
////
template<typename T>
double Analytics::statProdAvg (
  const DataSet<T> & dataSet1,
  const DataSet<T> & dataSet2,
  const typename vector<T>::size_type & beg,   // elements in the range
  const typename vector<T>::size_type & end    // [beg, end[
) {
  // the inner product between the two vectors
  double innerProd = 0;
  auto citr1 = dataSet1.beg + beg * dataSet1.incr,
       citr2 = dataSet2.beg + beg * dataSet2.incr;
  while (citr1 < dataSet1.beg + end * dataSet1.incr) {
    innerProd += *citr1 * *citr2;

    // increment iterators
    citr1 += dataSet1.incr;
    citr2 += dataSet2.incr;
  }

  // statistically average
  return innerProd / (end - beg);
}
template double Analytics::statProdAvg<int> (
  const DataSet<int> &,
  const DataSet<int> &,
  const vector<int>::size_type &,
  const vector<int>::size_type &
);

////
////
////
template<typename T>
double Analytics::statAvg (
  const DataSet<T> & dataSet,
  const typename vector<T>::size_type & beg,   // elements in the range
  const typename vector<T>::size_type & end    // [beg, end[
) {
  // the inner product between the two vectors
  double sum = 0;
  auto citr = dataSet.beg + beg * dataSet.incr;

  while (citr < dataSet.beg + end * dataSet.incr) {
    sum += *citr;

    // increment iterator
    citr += dataSet.incr;
  }

  // statistically average
  return sum / (end - beg);
}
template double Analytics::statAvg<double> (
  const DataSet<double> &,
  const vector<double>::size_type &,
  const vector<double>::size_type &
);
template double Analytics::statAvg<int> (
  const DataSet<int> &,
  const vector<int>::size_type &,
  const vector<int>::size_type &
);
template double Analytics::statAvg<unsigned> (
  const DataSet<unsigned> &,
  const vector<unsigned>::size_type &,
  const vector<unsigned>::size_type &
);



template<typename To, typename Ti1, typename Ti2>
vector<To> Analytics::elemProd (
  const DataSet<Ti1> &                    dataSet1,
  const DataSet<Ti2> &                    dataSet2,
  const typename vector<To>::size_type & size
  // vector<T> & result
) {
  // resize
  vector<To> result(size);
  // result.resize(size);

  auto itr1 = dataSet1.beg,
       itr2 = dataSet2.beg;
  auto itr3 = result.begin();
  while (itr1 < dataSet1.beg + size * dataSet1.incr) {
    *(itr3++) = *itr1 * *itr2;

    // increment iterators
    itr1 += dataSet1.incr;
    itr2 += dataSet2.incr;
  }

  return result;
}
template vector<int> Analytics::elemProd<int, int, int> (
  const DataSet<int> &,
  const DataSet<int> &,
  const vector<int>::size_type &
);
template vector<double> Analytics::elemProd<double, double, double> (
  const DataSet<double> &,
  const DataSet<double> &,
  const vector<double>::size_type &
);




template<typename T>
double Analytics::statStd (
  const DataSet<T> & dataSet,
  const typename vector<T>::size_type & beg,
  const typename vector<T>::size_type & end,
  const double & avg
) {
  // the inner product between the two vectors
  double innerProd = 0;
  auto itr = dataSet.beg + beg * dataSet.incr;

  while (itr < dataSet.beg + end * dataSet.incr) {
    innerProd += *itr * *itr;

    // increment iterator
    itr += dataSet.incr;
  }

  // statistically average
  return sqrt(innerProd / (end - beg) - avg * avg);
}
template double Analytics::statStd<double> (
  const DataSet<double> &,
  const vector<double>::size_type &,
  const vector<double>::size_type &,
  const double &
);