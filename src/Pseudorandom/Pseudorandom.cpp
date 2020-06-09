#include "Pseudorandom.h"

using namespace std;

Pseudorandom::Pseudorandom (
  uint_fast64_t defualtSeed,
  bool uniqueSeed
) {
  using namespace std::chrono;

  if (uniqueSeed) {
    ////
    //// create a new unique seed
    ////

    // nanoseconds since 1970
    high_resolution_clock::time_point now = high_resolution_clock::now();
    duration<unsigned long, nano> duration = duration_cast<nanoseconds>(now.time_since_epoch());

    // concatenate the "semi-unique" additional seed and convert to string
    this->seed = defualtSeed + hash<string> {} (to_string(duration.count()));
  } else {
    ////
    //// use the default seed
    ////

    this->seed = defualtSeed;
  }

  // seed the Mersenne Twister 64b engine
  this->mt.seed(this->seed);
}


uint_fast64_t Pseudorandom::getSeed () const {
  return this->seed;
}


template<typename T>
T Pseudorandom::Uint (T lowerBound, T upperBound) {
  if (settings::mode::debug) {
    if (lowerBound > upperBound) {
      cout << settings::cout::enterRed << "Pseudorandom::Uint: ERROR: lowerBound > upperBound  ->  EXIT" << settings::cout::resetStyle << endl;
      exit(EXIT_FAILURE);
    }
  }

  uniform_int_distribution<T> distribution(lowerBound, upperBound);
  return distribution(this->mt);
}
template unsigned Pseudorandom::Uint<unsigned> (unsigned, unsigned);
template int Pseudorandom::Uint<int> (int, int);
template long Pseudorandom::Uint<long> (long, long);


template<typename T>
T Pseudorandom::U (T lowerBound, T upperBound) {
  if (settings::mode::debug) {
    if (lowerBound > upperBound) {
      cout << settings::cout::enterRed << "Pseudorandom::Uint: ERROR: lowerBound > upperBound  ->  EXIT" << settings::cout::resetStyle << endl;
      exit(EXIT_FAILURE);
    }
  }

  uniform_real_distribution<T> distribution(lowerBound, upperBound);
  return distribution(this->mt);
}
template double Pseudorandom::U<double> (double, double);



template<typename T>
T Pseudorandom::Exp (
  const T & upperBound,
  const double & lambda
) {
  // make sure that the interval is defined properly
  if (settings::mode::debug && 0 > upperBound) {
    cout << settings::cout::enterRed << "Pseudorandom::Uint: ERROR: 0 > upperBound = " << upperBound << "  ->  EXIT" << settings::cout::resetStyle << endl;
    exit(EXIT_FAILURE);
  }

  // uniform drawn number between 0 and 1
  double u = this->U<double>(0, 1);



  // exponential drawn number between 0 and upperBound
  T generated;

  double exponent = -lambda * upperBound;
  if (exponent > 500) {
    // handle overflow from large exponent
    generated = upperBound - log((1 - u) * exp(-exponent) + u) / lambda;
  } else if (abs(lambda) < pow(10., -10.)) {
    // handle overflow from dividing by zero
    // expansion to first order in lambda
    generated = u * upperBound; //  - 0.5 * (1 - u) * u * upperBound * upperBound * lambda;
  } else {
    // no risk ov owerflow
    generated = - log(1 - u * (1 - exp(exponent))) / lambda;
  }


  if ( ! (generated >= 0 && generated <= upperBound)) {

    cout << "=========== (Pseudorandom::Exp)" << endl
         << "upperBound=" << upperBound << endl
         << "lambda=" << lambda << endl
         << (abs(lambda) <= pow(10., -10.)) << endl
         << generated << endl
         << exponent << endl
         << lambda << endl
         << u << endl
         << -lambda * (double) upperBound << endl
         << exp(-lambda * (double) upperBound) << endl
         << (1 - exp(-lambda * (double) upperBound)) << endl
         << 1 - u * (1 - exp(-lambda * (double) upperBound)) << endl
         << log(1 - u * (1 - exp(-lambda * (double) upperBound))) << endl
         << "===========" << endl;
    exit(EXIT_FAILURE);
  }

  // make sure that the generated number is inside the bounds
  if (settings::mode::debug && ! (generated >= 0 && generated <= upperBound)) {
    cout << settings::cout::enterRed << "Pseudorandom::Exp: ERROR: generated number outside of interval  ->  EXIT" << settings::cout::resetStyle << endl;
    cout << generated << endl;
    exit(EXIT_FAILURE);
  }

  return generated;
}
template long Pseudorandom::Exp<long> (const long &, const double &);
template double Pseudorandom::Exp<double> (const double &, const double &);