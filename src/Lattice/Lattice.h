#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <utility>
#include <algorithm>

#include "../etc/cout.h"
// #include "../Config/Config.h"

#include "../settings/settings.h"

class Lattice {
  public:
    const std::string name;

  private:

    // pointers to functions which are lattice dependent
    std::vector<unsigned> (Lattice::*i2isFuncPtr) (unsigned) const;
    void (Lattice::*is2RfuncPtr) (const std::vector<int> &,
                                  std::vector<double> &) const;

    // the spatial dimensions of the lattice
    unsigned numDimensions,
             numSites;

    // in which directions the lattice has periodic BC
    std::vector<bool> isPeriodic;

    std::vector<unsigned> size;

    // set of primitive vectors
    std::vector<std::array<double, 3> > pVs;

    // the positions of the lattice sites
    std::vector<std::vector<double> > Rs;

    // the nearest neighbors sites of each lattice site
    std::vector<std::vector<unsigned> > NNs;

    // the second nearest neighbors sites of each lattice site
    std::vector<std::vector<unsigned> > N2Ns;

    // the third nearest neighbors sites of each lattice site
    std::vector<std::vector<unsigned> > N3Ns;

    // the distance to the nearest neighbors
    std::vector<std::vector<std::vector<double> > > dists;

    // the boundary crossings in jumps from one site to another
    // indices: i*numSites + j (jumping from site j to i)
    // possible values: 0 (no crossing) or ±1 (crossing) in each dimension
    std::vector<std::vector<int> > BCs;

    // if a boundary is crossed going from j -> i
    std::vector<bool> BC;

  public:
    // Lattice (std::string,
    //          std::pair<std::vector<unsigned>, std::vector<double> >,
    //          Config &);

    Lattice (const std::string &,
             const std::vector<unsigned> &,
             const std::vector<double> &,
             const std::vector<bool> &);

    std::vector<int> boundaryCrossings (const unsigned,
                                        const unsigned) const;

    bool boundaryCrossed (const unsigned,
                          const unsigned) const;

    unsigned getNumDimensions () const;

    void getNNsAndDists (const unsigned,
                         std::vector<unsigned> &,
                         std::vector<std::vector<double> > &) const;

    bool getIsPeriodic (const unsigned) const;

    std::vector<bool> getIsPeriodic () const;

    bool getIsTranslInvariant () const;

    unsigned getNumNNs (const unsigned) const;

    void getNNs (const unsigned,
                 std::vector<unsigned> &) const;

    void getR (const unsigned,
               std::vector<double> &) const;
    std::vector<double> getR (const unsigned) const;

    void getRdiff (const unsigned,
                   const unsigned,
                   std::vector<double> &) const;

    unsigned getNumSites () const;

    std::vector<unsigned> getSize () const;

    std::vector<unsigned> getCenterSiteIs () const;


    std::vector<unsigned> i2is (unsigned) const;

    void getGhostSiteR (const unsigned,
                        const unsigned,
                        std::vector<double> &) const;

    std::vector<std::array<double, 3> > getPvs () const;

  private:

    void rectangular2d (unsigned,
                        unsigned,
                        double,
                        double);

    // used for the periodic boundary conditions
    static unsigned positiveModulo (const int,
                                    const int);

    // wrap each spatial lattice site index to be inside the domain
    void wrapIs (const std::vector<int> &,
                 std::vector<unsigned> &) const;

    std::vector<unsigned> wrapIs (const std::vector<int> &) const;


    ////
    //// 2D rectangular lattice
    ////
    unsigned rectangular2dIs2i (const std::vector<unsigned> &) const;

    template<typename T>
    std::vector<T> rectangular2dI2is (const unsigned) const;

    void rectangular2dIs2R (const std::vector<int> &,
                            std::vector<double> &) const;

    void rectangular2dNNsAndDists (const unsigned,
                                   std::vector<unsigned> &,
                                   std::vector<std::vector<double> > &) const;

    std::vector<int> rectangular2dBoundaryCrossings (const unsigned,
                                                     const unsigned) const;

    std::vector<unsigned> rectangular2dN2N (const unsigned) const;

    std::vector<unsigned> rectangular2dN3N (const unsigned) const;



    ////
    //// friend classes
    ////
    friend class Worm;

};

#endif