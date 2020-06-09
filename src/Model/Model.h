#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <utility>
#include <string>
#include <algorithm>

// #include <mpi.h>   // MPI

#include "../etc/cout.h"
#include "../Pseudorandom/Pseudorandom.h"
#include "../settings/settings.h"


////
//// The original Hamiltonin and the definition of parameters
////
////
//// H = ∑  ∑         + tᵢⱼ c⁺ᵃᵢ cᵃⱼ              // kinetic energy  (might differ a factor 2 or ½ to other definitions)
////     ᵃ ⁽ⁱʲ⁾
////     ∑  ∑         - μᵃᵢ nᵃᵢ                   // chemical potential
////     ᵃ  ⁱ
////     ∑  ∑         + ½ Uᵃ nᵃᵢ (nᵃᵢ - 1)        // intracomponent interaction
////     ᵃ  ⁱ
////     ∑    ∑       + Uᵃᵇ nᵃᵢ nᵇᵢ               // intercomponent interaction
////    ᵃ!⁼ᵇ  ⁱ
////     ∑  ∑         + ηᵃ cᵃᵢ + H.c              // worm operator
////     ᵃ  ⁱ
////

// numComps and numInters
using namespace settings::model;

class Model {
  public:

  private:

    // proper structures
    const std::array<bool, numComps> isBosonic, isCanonical;
    const std::array<unsigned, numComps> Ns;


    const unsigned numSites,
                   numDims;     // number of dimensions

    // whether or not mu and t depends on the position (i, j)
    const bool uniform;
    const bool equivalentComponents;

    const unsigned maxNumWorms;

    Pseudorandom & pseudoRandom;

    ////
    //// should be adjustable, hence non-constant
    ////
    double beta;
    std::vector<double> ts, Js, Us, mus, etas;

    // const double L = std::pow((double) V, 1 / (double) d);   // characteristic number of sites across system in some direction

    double (Model::*tFuncPtr)   (const unsigned, const std::vector<double> &, const std::vector<double> &) const;
    double (Model::*JfuncPtr)   (const unsigned, const unsigned, const std::vector<double> &, const std::vector<double> &) const;
    double (Model::*muFuncPtr)  (const std::vector<double> &, const unsigned);
    double (Model::*UfuncPtr)   (const unsigned, const unsigned) const;
    double (Model::*UnnFuncPtr) (const unsigned, const unsigned) const;
    double (Model::*etaFuncPtr) (const unsigned) const;


  public:
    Model (const std::string &,
           const std::array<bool, numComps> &,       // isBosonic
           const std::array<bool, numComps> &,       // isCanonical
           const std::array<unsigned, numComps> &,   // Ns
           const std::vector<double> &,
           const std::vector<double> &,
           const std::vector<double> &,
           const std::vector<double> &,
           const std::vector<double> &,
           const double &,                  // beta
           const unsigned,                  // numSites
           const unsigned,                  // numDims
           const bool,                      // uniform in space
           const bool,                      // equivalent components
           const unsigned,
           Pseudorandom &);                 // random number generator

    void loadParameters (const double &,
                         const std::vector<double> &,
                         const std::vector<double> &,
                         const std::vector<double> &,
                         const std::vector<double> &,
                         const std::vector<double> &);

    double t (const unsigned,
              const std::vector<double> &,
              const std::vector<double> &) const;

    double J (const unsigned,
              const unsigned,
              const std::vector<double> &,
              const std::vector<double> &) const;

    double mu (const std::vector<double> &,   // position vector
               const unsigned);               // component index

    double U (const unsigned,
              const unsigned) const;

    double U_nn (const unsigned,
                 const unsigned) const;

    double eta (const unsigned) const;

    bool hasEquivalentComponents () const;

    bool isUniform () const;

    std::array<bool, numComps> getIsBosonic () const;
    std::array<bool, numComps> getIsCanonical () const;
    std::array<unsigned, numComps> getNs () const;

    double getBeta () const;

    unsigned getMaxNumWorms () const;

  private:

    ////
    //// t-J model
    ////
    void tJ ();
    double tJ_t    (const unsigned, const std::vector<double> &, const std::vector<double> &) const,
           tJ_J    (const unsigned, const unsigned, const std::vector<double> &, const std::vector<double> &) const,
           tJ_mu   (const std::vector<double> &, const unsigned),
           tJ_U_nn (const unsigned, const unsigned) const,
           tJ_eta  (const unsigned) const;



    static std::vector<double> diff (const std::vector<double> &,
                                       const std::vector<double> &);

    static std::vector<double> sum (const std::vector<double> &,
                                      const std::vector<double> &);

    static double dot (const std::vector<double> &,
                       const std::vector<double> &);

    static double norm (const std::vector<double> &);
};

#endif