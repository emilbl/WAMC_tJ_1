#ifndef ANALYTICS_H
#define ANALYTICS_H

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>   // max_element
#include <numeric>     // accumulate
#include <string>
#include <utility>
#include <cmath>
#include <functional>   // std::function

#include "../etc/cout.h"
#include "../Pseudorandom/Pseudorandom.h"


class Analytics {
  public:

    template<typename T>
    struct DataSet {
      typename std::vector<T>::const_iterator beg;
      unsigned incr;
    };

    struct Compressed  {
      const double mean,
                   std,
                   mean_std;
      const unsigned corrLen;
    };

  private:

  public:

    // calculates the correlation length given a vector of numeric data
    template<typename T>
    static typename std::vector<T>::size_type correlationLength (const DataSet<T> &,                           // the data
                                                                 const typename std::vector<T>::size_type &,   // number of elements to include (size of list, padding excluded)
                                                                 const double &);                              // threshold

    // returns <x>, asymmetric std(x), symmetric std(<x>) and correlation length
    template<typename T>
    static Compressed compress (const DataSet<T> &,                           // the data
                                const typename std::vector<T>::size_type &,   // number of elements to include (size of list, padding excluded)
                                const std::string & quantityName = "");       // the name of quantity being compressed which might be used in dialog outputs

    // when we have to multiply with configuration weight sign
    template<typename T>
    static Compressed compress_signed (const DataSet<T> &,                           // the data
                                       const DataSet<T> &,                           // the signs
                                       const typename std::vector<T>::size_type &,   // number of elements to include (size of list, padding excluded)
                                       const std::string & quantityName = "");       // the name of quantity being compressed which might be used in dialog outputs

    template<typename T>
    static double bootstrap (const std::function<double(const std::vector<DataSet<T> > &,
                                                        const typename std::vector<T>::size_type &,
                                                        const typename std::vector<T>::size_type &)> &,
                             const std::vector<DataSet<T> > &,
                             const typename std::vector<T>::size_type &,
                             const double & meanVal,
                             const std::string & quantityName = "<unknown>");   // the name of quantity being bootstrapped which might be used in dialog outputs

    template<typename T>
    static double bootstrap (const std::function<double(const std::vector<DataSet<T> > &,
                                                        const typename std::vector<T>::size_type &,
                                                        const typename std::vector<T>::size_type &)> &,
                             const std::vector<DataSet<T> > &,
                             const typename std::vector<T>::size_type &,
                             const double &,                                         // the mean value of the quantity
                             const typename std::vector<T>::size_type &,        // correlation length
                             const std::string & quantityName = "<unknown>");   // the name of quantity being bootstrapped which might be used in dialog outputs

    template<typename T>
    static double statProdAvg (const DataSet<T> &,
                               const DataSet<T> &,
                               const typename std::vector<T>::size_type &,
                               const typename std::vector<T>::size_type &);

    template<typename T>
    static double statAvg (const DataSet<T> &,
                           const typename std::vector<T>::size_type &,
                           const typename std::vector<T>::size_type &);

    template<typename To, typename Ti1, typename Ti2>
    static std::vector<To> elemProd (const DataSet<Ti1> &,
                                    const DataSet<Ti2> &,
                                    const typename std::vector<To>::size_type &);

    template<typename T>
    static double statStd (const DataSet<T> &,
                           const typename std::vector<T>::size_type &,
                           const typename std::vector<T>::size_type &,
                           const double &);

    template<typename T>
    static void histogram (const DataSet<T> &,
                           const typename std::vector<T>::size_type &,
                           const unsigned,
                           std::vector<double> &,
                           std::vector<double> &);

  private:

    template<typename Ti, typename To>
    static bool divide (const std::function<To(const std::vector<DataSet<Ti> > &,
                                               const typename std::vector<Ti>::size_type &,
                                               const typename std::vector<Ti>::size_type &)> &,   // function to evaluate on the divided data sets
                        const std::vector<DataSet<Ti> > &,                                        // data sets
                        const typename std::vector<Ti>::size_type &,                              // size of the data set
                        const typename std::vector<To>::size_type &,                              // minimum number of segments
                        const typename std::vector<To>::size_type &,                              // minimum segment size
                        std::vector<To> &);                                                       // output data


};

#endif