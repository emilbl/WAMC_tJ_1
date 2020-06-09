#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <map>
#include <cmath>

// #include <mpi.h>   // MPI

#include "../etc/cout.h"
 #include "../settings/settings.h"
#include "../Lattice/Lattice.h"
#include "../Model/Model.h"
#include "../Worm/Segment/Segment.h"
#include "../Worm/Node/Node.h"



using namespace settings::model;

class Hamiltonian {

  public:
    const Lattice & lattice;
    Model & model;

  private:

    // if the chemical potential is uniform in space
    const bool equivalentComponents = model.hasEquivalentComponents();
    const bool uniform = model.isUniform();
    const std::array<bool, numComps> isCanonical = model.getIsCanonical();

    const unsigned numSites = lattice.getNumSites();

    // const double beta = model.getBeta();

    // const double etaFactor = 1 / sqrt(numSites * beta);

    // number of vector elements per jump
    const unsigned numVecElemsPerJump     = numComps;
    const unsigned numVecElemsPerExchange = std::pow(numComps, 2);

    std::vector<double> t_ija  = std::vector<double>(numVecElemsPerJump     * (uniform ? 1 : numSites * numSites), 0),
                        J_ijab = std::vector<double>(numVecElemsPerExchange * (uniform ? 1 : numSites * numSites), 0);


    // number of vector elements per site
    const unsigned numVecElemsPerSite = numComps;

    std::vector<double> mu_ia = std::vector<double>(numVecElemsPerSite * (uniform ? 1 : numSites), 0),   // original H quantity
                        V_ia  = std::vector<double>(numVecElemsPerSite * (uniform ? 1 : numSites), 0);   // transformed H quantity


    std::array<std::array<double, numComps>, numComps> U_inter_ab = {},   // original H quantity
                                                       U_ab       = {},   // transformed H quantity
                                                       U_nn_ab    = {};   // nearest neighbor potential interaction

    std::array<double, numComps> eta_a = {},
                                 U_a   = {};                              // original H quantity

  public:

    Hamiltonian (const Lattice &,   // the lattice
                 Model &);          // the model

    void loadParameters ();

    bool isUniform () const;

    bool hasEquivalentComponents () const;

    void setMu (const unsigned, const double &);

    std::vector<double> getMu () const;

    void correctChemicalPotential (const double &,                                          // beta
                                   const std::array<unsigned long long, 2 * numComps> &);   // Gs counts

    // used in displayLattice
    static unsigned i_ab (const unsigned,
                          const unsigned);

    // used in displayLattice
    std::vector<double>::size_type i_ia (const unsigned,
                                         const unsigned) const;

    std::vector<double>::size_type i_ija (const unsigned,
                                          const unsigned,
                                          const unsigned) const;

    std::vector<double>::size_type i_ijab (const unsigned,
                                           const unsigned,
                                           const unsigned,
                                           const unsigned) const;




    std::array<double, numComps> getChemicalPotential () const;

    std::vector<double> getMu_ia () const;
    std::vector<double> getV_ia () const;
    std::vector<double> gett_ija () const;
    std::vector<double> getJ_ijab () const;
    std::array<double, numComps> getU_a () const;
    std::array<std::array<double, numComps>, numComps> getU_inter_ab () const;
    std::array<std::array<double, numComps>, numComps> getU_ab () const;
    std::array<std::array<double, numComps>, numComps> getU_nn_ab () const;
    std::array<double, numComps> getEta_a () const;
    double getEta (unsigned) const;

    double getDicsoSquredFact (const unsigned a, const double &) const;

    void calculatePotenEnergies_inhom (const std::vector<unsigned long long> &,
                                       std::vector<double>::iterator,
                                       const bool = false) const;

    void calculatePotenEnergies_hom (const std::array<unsigned long long, numComps> &,
                                     std::vector<double>::iterator,
                                     const bool = false) const;

    void calculateInterEnergies (const std::array<unsigned long long, numInters> &,
                                 const std::array<unsigned long long, numComps> &,
                                 std::vector<double>::iterator,
                                 const bool = false) const;


    double disco (std::array<unsigned, numComps> &,
                  std::array<unsigned, numComps> &,
                  const double &) const;

    void kinet (unsigned,
                unsigned,
                std::array<unsigned, numComps> &,
                std::array<unsigned, numComps> &,
                std::array<unsigned, numComps> &,
                std::array<unsigned, numComps> &,
                double &) const;

    double jump (const unsigned,
                 const std::array<unsigned, numComps> &,
                 const std::array<unsigned, numComps> &,
                 const unsigned,
                 const std::array<unsigned, numComps> &,
                 const std::array<unsigned, numComps> &) const;


    double poten (const unsigned, const std::array<unsigned, numComps> &) const;
    double inter (const std::array<unsigned, numComps> &) const;

    double inter_nn (const std::array<unsigned, numComps> &,
                     const std::array<unsigned, numComps> &) const;


    ////
    //// local differences
    ////
    void potenDiff (const unsigned,                      // site index
                    const std::array<int, numComps> &,   // worm population
                    const int,                           // add (+1) or subtract (-1) worm pop
                    double &,                            // base
                    double &) const;                     // exponent
    void potenDiff (const unsigned,    // site index
                    const unsigned,    // active component
                    const int,         // worm population of active component
                    const int,         // add (+1) or subtract (-1) worm pop
                    double &,          // base
                    double &) const;   // exponent
    void potenDiff (const unsigned,                      // proposed site index
                    const unsigned,                      // current site index
                    const std::array<int, numComps> &,   // worm population
                    double &,                            // base
                    double &) const;                     // exponent
    void potenDiff (const unsigned,    // proposed site index
                    const unsigned,    // current site index
                    const unsigned,    // active component
                    const int,         // worm population of active component
                    double &,          // base
                    double &) const;   // exponent

    double interDiff_nn (const std::array<unsigned, numComps> &,
                         const std::array<unsigned, numComps> &,
                         const std::array<unsigned, numComps> &) const;


    double interDiff (const std::array<unsigned, numComps> &,   // site population before update
                      const std::array<int, numComps> &,        // worm population
                      const int) const;                         // add (+1) or subtract (-1) worm pop
    double interDiff (const std::array<unsigned, numComps> &,   // site population of active component before update
                      const unsigned,                           // active component
                      const int,                                // worm population of active component
                      const int) const;                         // add (+1) or subtract (-1) worm pop
    double interDiff (const std::array<unsigned, numComps> &,     // proposed site population before update
                      const std::array<unsigned, numComps> &,     // current site population before update
                      const std::array<int, numComps> &) const;   // worm population
    double interDiff (const std::array<unsigned, numComps> &,   // proposed site population of active component before update
                      const std::array<unsigned, numComps> &,   // current site population of active component before update
                      const unsigned,                           // active component
                      const int) const;                         // worm population of active component


    void discoDiff (const std::array<unsigned, numComps> &,   // site population before update
                    const std::array<int, numComps> &,        // worm population
                    const int,                                // add (+1) or subtract (-1) worm pop
                    const double &,                           // the current temperature
                    double &,                                 // base
                    double &) const;                          // exponent
    void discoDiff (const unsigned,    // site population of active component before update
                    const unsigned,    // active component
                    const int,         // worm population of active component
                    const int,         // add (+1) or subtract (-1) worm pop
                    const double &,    // the current temperature
                    double &,          // base
                    double &) const;   // exponent
    void discoDiffHead (const std::array<unsigned, numComps> &,   // the new population of the head segment
                        const std::array<unsigned, numComps> &,   // the current population of the head segment
                        const std::array<int, numComps> &,        // worm population
                        const int,                                // add (+1) or subtract (-1) worm pop
                        double &,                                 // base
                        double &) const;                          // exponent
    void discoDiffHead (const unsigned,    // the new population of the head segment
                        const unsigned,    // the current population of the head segment
                        const int,         // worm population of active component
                        const int,         // add (+1) or subtract (-1) worm pop
                        double &,          // base
                        double &) const;   // exponent
    void discoDiffTail (const std::array<unsigned, numComps> &,   // the new population of the segment before the tail segment
                        const std::array<unsigned, numComps> &,   // the current population of the segment before the tail segment
                        const std::array<int, numComps> &,        // worm population
                        const int,                                // add or subtract the worm population of the "new population" (0 do nothing)
                        double &,                                 // base
                        double &) const;                          // exponent
    void discoDiffTail (const unsigned,    // the new population of the segment before the tail segment
                        const unsigned,    // the current population of the segment before the tail segment
                        const int,         // worm population of active component
                        const int,         // add or subtract the worm population of the "new population" (0 do nothing)
                        double &,          // base
                        double &) const;   // exponent
    void discoDiffHeadOrTail (const unsigned,
                              const unsigned,
                              const int,
                              const int,
                              double &) const;


    void kinetDiff (const unsigned,                           // proposed site index
                    const unsigned,                           // current site index
                    const std::array<unsigned, numComps> &,   // proposed site population before update
                    const std::array<unsigned, numComps> &,   // current site population before update
                    const std::array<int, numComps> &,        // worm population
                    double &) const;                          // base
    void kinetDiff (const unsigned,    // proposed site index
                    const unsigned,    // current site index
                    const unsigned,    // proposed site population before update
                    const unsigned,    // current site population before update
                    const unsigned,    // active component
                    const int,         // worm population of active component
                    double &) const;   // base
    void kinetDiff (const unsigned,                              // lattice site index 1
                    const unsigned,                              // lattice site index 2
                    const std::array<unsigned, numComps> &,      // population at lattice site indexed 1 before jump
                    const std::array<unsigned, numComps> &,      // population at lattice site indexed 1 after jump
                    const std::array<unsigned, numComps> &,      // population at lattice site indexed 2 before jump
                    const std::array<unsigned, numComps> &,      // population at lattice site indexed 2 after jump
                    const std::array<int, numComps> &,           // worm population
                    const int,                                   // add or subtract worm to/from site indexed 1 before jump
                    const int,                                   // add or subtract worm to/from site indexed 1 before jump
                    const int,                                   // add or subtract worm to/from site indexed 2 before jump
                    const int,                                   // add or subtract worm to/from site indexed 2 before jump
                    double &) const;                             // base
    void kinetDiff (const unsigned,    // lattice site index 1
                    const unsigned,    // lattice site index 2
                    const unsigned,    // population at lattice site indexed 1 before jump
                    const unsigned,    // population at lattice site indexed 1 after jump
                    const unsigned,    // population at lattice site indexed 2 before jump
                    const unsigned,    // population at lattice site indexed 2 after jump
                    const unsigned,    // active component
                    const int,         // worm population of active component
                    const int,         // add or subtract worm to/from site indexed 1 before jump
                    const int,         // add or subtract worm to/from site indexed 1 before jump
                    const int,         // add or subtract worm to/from site indexed 2 before jump
                    const int,         // add or subtract worm to/from site indexed 2 before jump
                    double &) const;   // base

    void jumpDiff (const unsigned,
                   const std::array<unsigned, numComps> &,
                   const std::array<unsigned, numComps> &,
                   const unsigned,
                   const std::array<unsigned, numComps> &,
                   const std::array<unsigned, numComps> &,
                   const unsigned,
                   const int,
                   double &) const;

  private:

};

#endif