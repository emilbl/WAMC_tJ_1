#include "Analytics.h"

using namespace std;

////
////
////
template<typename T>
double Analytics::bootstrap (
  const function<double(const vector<DataSet<T> > &,
                        const typename vector<T>::size_type &,
                        const typename vector<T>::size_type &)> & quantity,
  const vector<DataSet<T> > &                                     dataSets,
  const typename vector<T>::size_type &                           size,
  const double &                                                  meanVal,       // mean value of the quantity
  const string &                                                  quantityName   // the name of quantity being bootstrapped which might be used in dialog outputs
) {
  ////
  //// settings
  ////
  const unsigned minNumSegs = 1000,
                 numBootstrapIters = 1000;

  ////
  //// calculate the (maximum) correlation length for the data
  ////
  vector<unsigned> taus = {};
  for (const auto dataSet : dataSets) {
    const auto tau = Analytics::correlationLength(dataSet,
                                                  size,
                                                  exp(-1));

    taus.push_back(tau);
  }
  const auto tau = *max_element(taus.begin(), taus.end());

  ////
  //// divide data into n segments
  ////
  vector<double> meanVals;
  const bool succDiv = Analytics::divide<T, double>(quantity,
                                                    dataSets,
                                                    size,
                                                    minNumSegs,
                                                    2 * tau,
                                                    meanVals);

  // if the division did go to plan or not
  if ( ! succDiv) {
    cout << "Analytics::bootstrap: problem \"dividing\" " << quantityName << endl;
  }

  ////
  //// seed the random number generator
  ////
  Pseudorandom pseudoRandom{(uint_fast64_t) meanVal, true};

  ////
  //// do the bootstrapping
  ////
  const unsigned meanValsSize = meanVals.size();
  double mean_std = 0;
  for (unsigned i = 0; i < numBootstrapIters; i++) {
    double sum = 0;
    for (unsigned j = 0; j < meanValsSize; j++) {
      unsigned k = pseudoRandom.template Uint<unsigned>(0, meanValsSize - 1);
      sum += meanVals[k];
    }
    mean_std += pow(sum / (double) meanValsSize - meanVal, 2.);
  }

  return sqrt(mean_std / (double) numBootstrapIters);
}
template double Analytics::bootstrap<double> (
  const function<double(const vector<DataSet<double> > &,
                        const typename vector<double>::size_type &,
                        const typename vector<double>::size_type &)> &,
  const vector<DataSet<double> > &,
  const vector<double>::size_type &,
  const double & meanVal,
  const string &
);
template double Analytics::bootstrap<int> (
  const function<double(const vector<DataSet<int> > &,
                        const typename vector<int>::size_type &,
                        const typename vector<int>::size_type &)> &,
  const vector<DataSet<int> > &,
  const vector<int>::size_type &,
  const double & meanVal,
  const string &
);


////
////
////
template<typename T>
double Analytics::bootstrap (
  const function<double(const vector<DataSet<T> > &,
                        const typename vector<T>::size_type &,
                        const typename vector<T>::size_type &)> & quantity,      // quantity of interest, a function of the data
  const vector<DataSet<T> > &                                     dataSets,      // data
  const typename vector<T>::size_type &                           size,          // size of data
  const double &                                                  meanVal,       // mean value of the quantity
  const typename vector<T>::size_type &                           tau,           // correlation length
  const string &                                                  quantityName   // the name of quantity being bootstrapped which might be used in dialog outputs
) {
  ////
  //// settings
  ////
  const unsigned minNumSegs = 1000,
                 numBootstrapIters = 1000;

  ////
  //// divide the data into n segments and evaluate the quantity on these segments
  ////
  vector<double> meanVals;
  const bool succDiv = Analytics::divide<T, double>(quantity,
                                                    dataSets,
                                                    size,
                                                    minNumSegs,
                                                    2 * tau,
                                                    meanVals);

  // if the division did go to plan or not
  if ( ! succDiv) {
    cout << "Analytics::bootstrap: problem \"dividing\" " << quantityName << endl;
  }

  ////
  //// seed the random number generator
  ////
  Pseudorandom pseudoRandom{(uint_fast64_t) meanVal, true};

  ////
  //// do the bootstrapping
  ////
  const unsigned meanValsSize = meanVals.size();
  double mean_std = 0;
  for (unsigned i = 0; i < numBootstrapIters; i++) {
    double sum = 0;
    for (unsigned j = 0; j < meanValsSize; j++) {
      unsigned k = pseudoRandom.template Uint<unsigned>(0, meanValsSize - 1);
      sum += meanVals[k];
    }
    mean_std += pow(sum / (double) meanValsSize - meanVal, 2.);
  }

  return sqrt(mean_std / (double) numBootstrapIters);
}
template double Analytics::bootstrap<double> (
  const function<double(const vector<DataSet<double> > &,
                        const typename vector<double>::size_type &,
                        const typename vector<double>::size_type &)> &,
  const vector<DataSet<double> > &,
  const vector<double>::size_type &,
  const double &,
  const typename vector<double>::size_type &,
  const string &
);
template double Analytics::bootstrap<int> (
  const function<double(const vector<DataSet<int> > &,
                        const typename vector<int>::size_type &,
                        const typename vector<int>::size_type &)> &,
  const vector<DataSet<int> > &,
  const vector<int>::size_type &,
  const double &,
  const typename vector<int>::size_type &,
  const string &
);
template double Analytics::bootstrap<unsigned> (
  const function<double(const vector<DataSet<unsigned> > &,
                        const typename vector<unsigned>::size_type &,
                        const typename vector<unsigned>::size_type &)> &,
  const vector<DataSet<unsigned> > &,
  const vector<unsigned>::size_type &,
  const double &,
  const typename vector<unsigned>::size_type &,
  const string &
);