#ifndef SETTINGS_H
#define SETTINGS_H

#include <limits.h>   // ULONG_MAX
#include <string>
#include <cmath>
#include <limits.h>   // ULONG_MAX

// precompiled parameters and other settings
namespace settings {


  // program mode
  namespace mode {

    // tries to go back after each update + checks other details in the update methods
    constexpr bool debug = false;

    // recomputes and compare worm weight difference, worm configuration and quantities
    constexpr bool debug_major = false;

    constexpr unsigned long long debugFrom = 0; // 1e0;

    // perform some checks every Nth update
    constexpr unsigned long long verifyEvery = 1e8;


    constexpr bool shutItDown = false;
    constexpr bool verbose    = false;


  }


  // model settings
  namespace model {

    // choose the model you want to simulate
    enum Type { tJ };

    constexpr unsigned numComps = 2;

    constexpr Type modelType    = tJ;
    const std::string modelName = "tJ";

    constexpr bool has_t    = true;
    constexpr bool has_J    = true;
    constexpr bool has_U    = false;
    constexpr bool has_U_nn = true;
    constexpr bool has_mu   = true;


    constexpr unsigned numInters         = (numComps * numComps + numComps) / 2;   // diagonal + above diagonal
    constexpr unsigned numExternalInters = (numComps * numComps - numComps) / 2;   // above diagonal
  }


  // precompiled worm parameters
  namespace worm {
    // the discretization. absolutely no greater than ULONG_MAX/6
    // must not be unsigned since used in calculations with signed integers
    const long tMax = std::min((long) 1e10, LONG_MAX/100);

    // the reserved size of the segment containers
    constexpr unsigned initialSegmentContainerSize = 1E3;   // linear dependency on beta?

    // sensitivity threshold value for which the ratio in
    // Worm::checkDetailedBalance should not be larger
    constexpr double detailedBalanceMaxDiff = 1.E-10;

    // sensitivity threshold value for the acceptance
    // ratio difference in the anti-update test
    constexpr double diffTresh = 1.E-10;

    // the reserved size of data containers
    // after which they need to be reallocated
    constexpr unsigned long long reservedHistSize = 1048576;

    // whether multicomponent worm should be enabled (true) or not (false)
    constexpr bool allowMultiCompWorm = false;

    // whether worms with negative populations should be enabled (true) or not (false)
    // (this only makes sense when allowMultiCompWorm = false)
    constexpr bool allowAntiWorm = false;

    // whether an exponential distribution (true) or a uniform
    // distribution (false) should be used in the update procedures
    constexpr bool expDistrEnabled = true;
  }


  // styles for the cout (ANSI escape code)
  namespace cout {
    const std::string enterBlue   = "\033[38;5;12m",
                      enterYellow = "\033[38;5;11m",
                      enterGreen  = "\033[38;5;10m",
                      enterRed    = "\033[38;5;196m",
                      enterOrange = "\033[38;5;202m",
                      enterBold   = "\033[1m",
                      resetStyle  = "\033[0m";
  }
}

#endif