#ifndef PSEUDORANDOM_H
#define PSEUDORANDOM_H

#include <iostream>
#include <random>       // mt19937_64, etc.
#include <chrono>       // seed mt
#include <string>
#include <sstream>      // ostringstream
#include <functional>   // hash
#include <mpi.h>

#include "../settings/settings.h"

class Pseudorandom {

  public:

  private:

    uint_fast64_t seed;
    std::mt19937_64 mt;

  public:
    Pseudorandom (uint_fast64_t, bool);

    template<typename T>
    T Uint (T lowerBound, T upperBound);

    template<typename T>
    T U (T lowerBound, T upperBound);

    template<typename T>
    T Exp (const T &,
           const double &);

    uint_fast64_t getSeed () const;

  private:
};

#endif