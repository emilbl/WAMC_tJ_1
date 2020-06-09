#ifndef WORM_H
#define WORM_H

#include <iostream>
#include <iomanip>     // setw
#include <sstream>
#include <vector>
#include <memory>
#include <chrono>      // seed mt
#include <algorithm>   // any_of
#include <sstream>
#include <string>
#include <utility>
#include <map>

// #include <mpi.h>   // MPI

#ifdef __clang__
#pragma GCC diagnostic push                      // disables warning temporarily
#pragma GCC diagnostic ignored "-Wunsequenced"   //
#endif
#include <nlohmann/json.hpp>                     // https://github.com/nlohmann/json
#ifdef __clang__
#pragma GCC diagnostic pop                       // enables warning once again
#endif


#include "../settings/settings.h"
#include "../etc/cout.h"
#include "../etc/sgn.h"
#include "../Lattice/Lattice.h"
#include "Segment/Segment.h"
#include "Node/Node.h"
#include "../Model/Model.h"
#include "../Hamiltonian/Hamiltonian.h"
#include "../Pseudorandom/Pseudorandom.h"
#include "../Write-to-file/Write-to-file.h"
#include "../Analytics/Analytics.h"
#include "../Tic/Tic.h"

extern std::string dateAndTime;

using namespace settings::mode;
using namespace settings::worm;
using namespace settings::model;

template<typename Ti, typename To>
extern To positiveModulo (
  const Ti i,
  const Ti n
) {
  return (i % n + n) % n;
}

class Worm {
  public:

  private:
    ////
    //// settings
    ////

    static constexpr long long unsigned numUpdatesPerTimeCheck = debug_major ? 1000 : 1e6;
    static constexpr double printInterval = 10 * 60;



    ////
    //// constructor arguments and constants dependent upon these arguments
    ////
    const bool           devMode;
    Hamiltonian &        H;
    const Lattice &      lattice;
    Pseudorandom &       pseudoRandom;
    const WriteToFile &  writeToFile;
    const double         beta_target;
    double               beta;

    const std::array<bool, numComps>       isBosonic;
    const std::array<bool, numComps>       isCanonical;
    const std::array<unsigned, numComps>   Ns;
    /*const*/ std::array<double, numComps> maxNsDiff;
    const std::array<bool, numComps>       muOptimization;
    const std::array<double, numComps>     targetAvgNsDiff;
    // converted to discretized time
    std::array<long, numComps> _maxNsDiff;

    // whether or not we need to check for fermionic exchange
    const bool hasFermionicExchange;

    const unsigned maxNumWorms;

    const unsigned              numSites = lattice.getNumSites();
    const unsigned              numDims  = lattice.getNumDimensions();
    const double                L = pow( (double) numSites, 1 / (double) numDims);
    const std::vector<unsigned> latticeSize = lattice.getSize();
    const std::vector<bool>     isPeriodic = lattice.getIsPeriodic();

    // if the lattice is translationally invariant in any direction
    const bool isTranslInvariant = lattice.getIsTranslInvariant();


    const unsigned largeDataSamplingPeriod;
    const bool sample_sign;
    const bool sample_numParticles;
    const bool sample_numWinds;
    const bool sample_kinetEnergy;
    const bool sample_exchaEnergy;
    const bool sample_potenEnergy;
    const bool sample_interEnergy;
    const bool sample_totEnergy;
    const bool sample_nnInterEnergy;
    const bool sample_numParticlesAtSite;
    const bool sample_flow;
    const bool sample_instantParticleNum;
    const bool sample_instantParticleNumConv;
    const bool sample_GreensFunction;
    const bool sample_fourPointCorrelator;
    const bool sample_densityDensityCorrelator;

    // all energy components must be sampled in order for the total energy to be sampled as well
    const bool sample_totalEnergy = sample_totEnergy &&
                                    sample_kinetEnergy &&
                                    sample_potenEnergy &&
                                    ( ! has_J    || sample_exchaEnergy) &&
                                    ( ! has_U    || sample_interEnergy) &&
                                    ( ! has_U_nn || sample_nnInterEnergy);

    // if we are sampling any quantity other than the Green's function
    const bool sample_hist =    sample_sign
                             || sample_numParticles
                             || sample_numWinds
                             || sample_kinetEnergy
                             || sample_exchaEnergy
                             || sample_potenEnergy
                             || sample_interEnergy
                             || sample_totEnergy
                             || sample_nnInterEnergy
                             || sample_numParticlesAtSite
                             || sample_flow
                             || sample_instantParticleNum
                             || sample_instantParticleNumConv;


    ////
    //// computable quantities independent of sampling quantities
    ////
    const bool compute_numParticles;
    const bool compute_numWinds;
    const bool compute_kinetEnergy;
    const bool compute_exchaEnergy;
    const bool compute_potenEnergy;
    const bool compute_interEnergy;
    const bool compute_nnInterEnergy;
    const bool compute_osParticleProd;
    const bool compute_nnParticleProd;
    const bool compute_nnnParticleProd;
    const bool compute_C1;
    const bool compute_C2;
    const bool compute_C3;
    const bool compute_particleCorr;
    const bool compute_particleConv;
    const bool compute_spatialCorrelators =  compute_C1
                                          || compute_C2
                                          || compute_C3
                                          || compute_particleCorr
                                          || compute_particleConv;

    // weight using all of imaginary time
    const bool compute_average;


    ////
    ////
    ////
    const bool saveWarmUpData;
    const bool savePresimulationData;


    unsigned numTransl;


    ////
    //// debug quantities
    ////

    // the following variables are used for the anti-update test
    long t_prev,
         t1_prev,
         t2_prev,
         tMini_prev,
         tMinj_prev,
         tMin_prev,
         tMaxi_prev,
         tMaxj_prev,
         tMax_prev,
         intLength_prev,
         modifiedLength_prev,
         modifiedLength1_prev,
         modifiedLength2_prev,
         tMin_old_prev,
         tMax_old_prev;
    unsigned i_prev,
             j_prev,
             actComp_prev;
    int actPop_prev;
    double lambda_prev,
           R_prev;
    std::array<int, numComps> wormPop_prev;
    bool subtract_prev,
         goingForward_prev;

    // used in Worm::preCheckDetailedBalance
    int preWormWeightSign = 1;
    double preWormWeightBase = 1,
           preWormWeightExponent = 0;

    ////
    //// the current system configuration
    ////

    // one vector for each site holding the segments in a chronological order
    std::vector<std::vector<std::shared_ptr<Segment> > > sites;

    // the number of worms present in the system
    unsigned numWorms = 0;

    // holds the head and tail segments of the worm
    std::vector<std::shared_ptr<Segment> > headSegments, tailSegments;

    // corresponding segment index in the site vector of the head and tail
    std::vector<unsigned> headSegmentIndices, tailSegmentIndices;

    // in the case of single component worm (allowMultiCompWorm = false), this stores the active component index
    // to avoid having to loop through all components in various calculations
    std::vector<unsigned> actComps;

    // the worm population of the active component
    std::vector<int> actPops;

    /* TO BE UPDATED */   // holds the population of the active worm
    /* TO BE UPDATED */   std::array<int, numComps> wormPop = {};


    ////
    //// data extracted from the current system configuration
    ////

    bool recomputeFermionicExchangeSign = true;
    int prevFermionicExchangeSign = 1;

    // sign of the weight
    int sign = 1;

    // current number of particles in the system per component and site (times tMax)
    std::vector<unsigned long long> numParticlesAtSite = std::vector<unsigned long long>(numSites * numComps);

    // current number of particles in the system per component
    std::array<unsigned long long, numComps> numParticles = {};

    // current number of particles squared (inter and intra) in the system (times tMax)
    // (n_a * n_b)[a, b]
    std::array<unsigned long long, numInters> numParticlesSquared = {};

    // // current number of nearest neighbor particles squared (inter and intra) in the system (times tMax)
    // // (n_a,i * n_b,j)[a, b]
    // std::array<unsigned long long, numInters> numParticlesSquared_nn = {};

    // number of kinks
    unsigned numKinks = 0;

    // current number of jumps in the system
    std::array<unsigned, numComps> numJumps = {};

    // current number of exchanges in the system
    std::array<unsigned, numExternalInters> numExchanges = {};

    // the current winding number in the system
    std::array<std::vector<int>, numComps> numWinds;

    // the current particle "flow" at each site in the system
    // [site][component][spatial dimension]
    // particle -> increase in particular direction
    // hole -> decrease in particular direction
    std::vector<double> flow = std::vector<double>(numSites * numComps * numDims);

    // nearest neighbor interaction energy
    double U_nn = 0;

    array<double, 2 * numComps> muIntervals;

    ////
    //// history of data extracted from system configurations up until now
    ////

    unsigned long posSignCount = 0,
                  negSignCount = 0;

    ////
    //// small data history vectors
    ////
    std::vector<int>     signHist_small;
    std::vector<int>     numParticlesHist;
    std::vector<int>     numWindsHist;
    std::vector<double>  kinetEnergyHist,
                         exchaEnergyHist,
                         potenEnergyHist,
                         interEnergyHist,
                         totEnergyHist,
                         nnInterEnergyHist;

    ////
    //// large data history vectors
    ////
    std::vector<int>      signHist_large;
    std::vector<double>   numParticlesAtSiteHist;
    std::vector<double>   flowHist;
    std::vector<unsigned> instantParticleNumHist;
    std::vector<unsigned> instantParticleNumConvHist;


    ////
    ////
    ////
    long long Z_bin;
    unsigned long long Z_bin_count;

    ////
    //// Green's function histogram
    ////
    // number of bind per translation of the Green's function
    static constexpr unsigned G_hist_N = 200;  // make even!
    std::array<std::vector<long>, numComps> G_hist;
    std::array<std::vector<unsigned long>, numComps> G_hist_count;


    ////
    //// four point correlator histogram
    ////
    static constexpr unsigned fpc_hist_N = 20; // in each dimension: T, t1, t2
    std::array<std::vector<long>, numComps> fpc_hist;
    std::array<std::vector<unsigned long>, numComps> fpc_hist_count;


    ////
    //// density-density correlator
    ////
    static constexpr unsigned ddc_hist_N_o = 20;   // length of overlap
    static constexpr unsigned ddc_hist_N_d = 20;   // length of worm not contributing to overlap
    std::array<std::vector<long>, numComps> ddc_hist;
    std::array<std::vector<unsigned long>, numComps> ddc_hist_count;


    ////
    //// compute quantities
    ////
    std::vector<double> avgNumParticlesOne     = std::vector<double>(4 * numComps);
    std::vector<double> avgNumParticlesTwo     = std::vector<double>(4 * numInters);
    std::vector<double> avgNumParticlesAll     = std::vector<double>(4 * 1);
    std::vector<double> avgNumWinds            = std::vector<double>(4 * numInters);
    std::vector<double> avgNumWindsCyclic      = std::vector<double>(4 * numComps);
    std::vector<double> avgCounterflow         = std::vector<double>(4 * 1);
    std::vector<double> avgPairwiseCounterflow = std::vector<double>(4 * 1);
    std::vector<double> avgCoCounterflow       = std::vector<double>(4 * 1);
    std::vector<double> avgCoflow              = std::vector<double>(4 * 1);
    std::vector<double> avgKinetEnergy         = std::vector<double>(4 * numComps);
    std::vector<double> avgExchaEnergy         = std::vector<double>(4 * numExternalInters);
    std::vector<double> avgPotenEnergy         = std::vector<double>(4 * numComps);
    std::vector<double> avgInterEnergy         = std::vector<double>(4 * numInters);
    std::vector<double> avgNnInterEnergy       = std::vector<double>(4 * 1);

    std::vector<double> osParticleProd         = std::vector<double>(4 * numInters);
    std::vector<double> nnParticleProd         = std::vector<double>(4 * numInters);
    std::vector<double> nnnParticleProd        = std::vector<double>(4 * numInters);


    ////
    //// average NNs correlations
    ////
    const unsigned halfNumNNs = 0.5 * lattice.getNumNNs(0);
    const std::vector<unsigned> Ls = lattice.getSize();
    const std::vector<unsigned> Cs = lattice.getCenterSiteIs();
    std::vector<std::shared_ptr<Segment> > segments = std::vector<std::shared_ptr<Segment> >(numSites);
    std::vector<int> instC1  = std::vector<int>(numSites * halfNumNNs),
                     _instC1 = std::vector<int>(numSites * halfNumNNs);
    std::vector<double> avgC1; // = std::vector<double>(4 * numSites * halfNumNNs);   // save signs separate

    std::vector<int> instC2  = std::vector<int>(numSites * halfNumNNs),
                     _instC2 = std::vector<int>(numSites * halfNumNNs);
    std::vector<double> avgC2; // = std::vector<double>(4 * numSites * halfNumNNs);   // save signs separate

    std::vector<int> instC3  = std::vector<int>(numSites * halfNumNNs),
                     _instC3 = std::vector<int>(numSites * halfNumNNs);
    std::vector<double> avgC3; // = std::vector<double>(4 * numSites * halfNumNNs);   // save signs separate

    std::vector<unsigned> instantParticleNums  = std::vector<unsigned>(numComps * numSites),
                          _instantParticleNums = std::vector<unsigned>(numComps * numSites);
    std::vector<double> avgParticleCorr = std::vector<double>(4 * numComps * numSites);   // save signs separate

    std::vector<unsigned long> avgParticleConv = std::vector<unsigned long>(4 * (numComps + 1) * numSites);






    // in the case of a canonical component, this keeps track of how many Gs configuration in that component are
    // below as well as above the specified particle number in order to automatically optimize the chemical potential
    // std::array<unsigned long long, 2 * numComps> numGsCount = {};

    // used to optimize the chemical components for the canonical components
    std::array<double, numComps> avgParticleDiff          = {};
    std::array<long long, numComps> avgParticleDiff_count = {};

    // simulated annealing statistics
    nlohmann::json warmupStatistics;

    ////
    //// other member variables
    ////
    std::array<unsigned long long, numComps> wormCount;

    // update and bin counter
    unsigned long long currUpdateCount = 0;
    unsigned long long currGsCount = 0;
    unsigned long long currZsCount = 0;
    unsigned long currBinCount0 = 0;
    unsigned long currBinCount1 = 0;
    unsigned long currBinCount2 = 0;

    // print timer
    Tic printTimer{};

    // how frequently one should bin. Will be modified by the warm up.
    unsigned long samplingPeriod = 1;
    // unsigned long long minNumUpdatesPerZ = 9602;
    // unsigned long long updateNumAtPrevBin = 1;

    // how many times one has tried to bin but being denied due to the bin period time
    unsigned currTriedBinCount = 0;


    // weights for choosing the corresponding update procedure (to be read from file)
    unsigned Wremove,
             Wjump,
             WantiJump,
             Wreconnect,
             WantiReconnect,
             WtimeShift,
             Wovertake,
             Wfollow,
             Wtot;

    // what is the current ongoing update
    unsigned currentUpdate;

    // select, progress and accept frequency for each update
    std::array<unsigned long long, 9 * 3> updateStatistics;

    // to be filled with update procedures corresponding to their weight
    std::vector<void (Worm::*)(const double &, const unsigned, const double &, const bool)> weightedUpdateProcedures = {};

  public:


    Worm (Hamiltonian &,
          const Lattice &,
          const Model &,
          Pseudorandom &,
          const WriteToFile &,
          const nlohmann::json &,
          const nlohmann::json &,
          const nlohmann::json &,
          const std::array<double, numComps> &,
          const std::array<bool, numComps> &,
          const std::array<double, numComps> &,
          const bool);

    void loadConfiguration (const nlohmann::json &);

    ////
    //// getters
    ////
    const Lattice & getLattice () const;
    const std::vector<std::shared_ptr<Segment> > & getSegments (unsigned) const;
    double getBeta () const;
    double getBeta_target () const;
    void setBeta (const double &);
    void setMaxNsDiff (const std::array<double, numComps> &);
    long gettMax () const;
    unsigned long long getCurrUpdateCount () const;
    unsigned getNumWorms () const;
    std::array<int, numComps> getWormPop (const unsigned) const;

     /* OLD ? */ bool isHeadSegment (std::shared_ptr<Segment>) const;
     /* OLD ? */ bool isTailSegment (std::shared_ptr<Segment>) const;
     /* OLD ? */ bool isActiveComponent (unsigned) const;

    std::shared_ptr<Segment> getHeadSegment (const unsigned) const;
    std::shared_ptr<Segment> getTailSegment (const unsigned) const;

    void shutDown () const;



    void searchOptimalMu (const double &,
                          const unsigned,
                          double &,
                          double &,
                          const double &,
                          const unsigned long &,
                          const unsigned,
                          unsigned long &,
                          unsigned long &,
                          const long unsigned &,
                          const unsigned);
    void searchOptimalMu (const double &,
                          const unsigned,
                          double &,
                          double &,
                          const double &,
                          const unsigned long &,
                          const unsigned);

    double testMu (const double &,
                   const unsigned,
                   const double &,
                   const unsigned long &,
                   const double &,
                   unsigned long &,
                   unsigned long &,
                   const long unsigned &,
                   const unsigned);

    void warmUp (const unsigned,
                 const double &,
                 const unsigned);

    void singleUpdate ();
    void multipleUpdates (const unsigned,
                          const unsigned,
                          const unsigned long &);

    void optimizeMu (const unsigned long &,
                     const double &,
                     const array<unsigned long long, 2 * numComps> &);

  private:

    static unsigned i_ab (const unsigned,
                          const unsigned);
    static unsigned i_ab_external (const unsigned,
                                   const unsigned);

    double int2time (const double &) const;

    // returns the later segment
    std::shared_ptr<Segment> splitSegment (const unsigned,         // lattice site index
                                           const unsigned,         // segment index
                                           const unsigned long);   // splitting time

    void removeSegment (const unsigned,    // lattice site index
                        const unsigned);   // segment index

    void update ();

    void insert        (const double &, const double &, const bool),
         remove        (const double &, const unsigned, const double &, const bool),
         jump          (const double &, const unsigned, const double &, const bool),
         antiJump      (const double &, const unsigned, const double &, const bool),
         reconnect     (const double &, const unsigned, const double &, const bool),
         antiReconnect (const double &, const unsigned, const double &, const bool),
         timeShift     (const double &, const unsigned, const double &, const bool),
         overtake      (const double &, const unsigned, const double &, const bool),
         follow        (const double &, const unsigned, const double &, const bool);

    double getInsertProposalProbability (const std::vector<unsigned> &) const;

    // axillary function for insert and remove
    void proposeWormPopulationAndInsertionType (const std::vector<unsigned> &,
                                                const std::array<unsigned, numComps> &,
                                                std::array<int, numComps> &,
                                                unsigned &,
                                                int &,
                                                double &,
                                                bool &) const;
    void proposeWormPopulationAndInsertionType (const std::vector<unsigned> &,
                                                const std::array<unsigned, numComps> &,
                                                const unsigned,
                                                double &) const;


    // computes the nearest neighbor interaction energy for the who;e system
    double compute_U_nn () const;

    double computeNNinteraction (const unsigned,
                                 const std::array<unsigned, numComps> &,
                                 const long &,
                                 const long &) const;


    double computeNNinteractionDiff (const unsigned,
                                     const std::array<unsigned, numComps> &,
                                     const unsigned,
                                     const int,
                                     const long &,
                                     const long &);
    double computeNNinteractionDiff (const unsigned,
                                     const std::array<unsigned, numComps> &,
                                     const std::array<unsigned, numComps> &,
                                     const long &,
                                     const long &);
    double computeNNinteractionDiff (const unsigned,
                                     const std::array<unsigned, numComps> &,
                                     const unsigned,
                                     const int,
                                     const long &,
                                     const long &,
                                     const unsigned);
    double computeNNinteractionDiff (const unsigned,
                                     const std::array<unsigned, numComps> &,
                                     const std::array<unsigned, numComps> &,
                                     const long &,
                                     const long &,
                                     const unsigned);
    double computeNNinteractionDiff (const unsigned,
                                     const std::array<unsigned, numComps> &,
                                     const std::array<unsigned, numComps> &,
                                     const unsigned,
                                     const std::array<unsigned, numComps> &,
                                     const std::array<unsigned, numComps> &,
                                     const long &,
                                     const long &);
    double computeNNinteractionDiff (const unsigned,
                                     const std::array<unsigned, numComps> &,
                                     const unsigned,
                                     const int,
                                     const unsigned,
                                     const std::array<unsigned, numComps> &,
                                     const unsigned,
                                     const int,
                                     const long &,
                                     const long &);


    void computeWormQuantities (unsigned &,
                                double &,
                                std::array<unsigned, numComps> &,
                                std::array<unsigned, numExternalInters> &,
                                std::array<unsigned long long, numComps> &,
                                std::array<unsigned long long, numInters> &,
                                std::array<std::vector<int>, numComps> &,
                                std::vector<double> &,
                                std::vector<unsigned long long> &) const;


    void findSegmentIndex (const unsigned,        // lattice site index
                           const long,            // time
                           unsigned &) const;     // the segment index


    void findSegment (const unsigned,                          // lattice site index
                      const long,                              // time
                      std::shared_ptr<Segment> &,   // the segment
                      unsigned &) const;                       // the segment index

    void timeModulo (long &) const;

    void checkWormConfiguration () const;
    void checkWormQuantities () const;

    void preCheckWormWeight (const double &);
    void checkWormWeight (const double &,
                          const double &,
                          const double &) const;

    void evaluateWorm (const double &,
                       int &,
                       double &,
                       double &) const;

    void trySampleData (const double &);
    void sampleZ ();
    void sampleG ();
    void sampleFPC ();
    bool sampleDDC ();
    void sampleSpatialCorrelators ();
    void sampleSpatialCorrelators_slice ();
    void sampleSpatialCorrelators_average ();
    void sampleSpatialCorrelators_computeAndStore (const unsigned,
                                                   const std::vector<unsigned> &,
                                                   const unsigned,
                                                   const double,
                                                   const bool);
    void sampleSpatialCorrelators_store (const unsigned,
                                         const unsigned,
                                         const double,
                                         const bool);

    void resetDataContainers ();

    std::tuple<std::shared_ptr<Segment>, std::shared_ptr<Segment>, std::shared_ptr<Segment>, std::shared_ptr<Segment> >
      findHeadsAndTails (const unsigned) const;

    unsigned spatialIndexDifference (const unsigned,
                                     const unsigned) const;

    int computeFermionicExchangeSign () const;


    template<typename To, typename Ti1, typename Ti2>
    static std::array<To, numComps> diff (const std::array<Ti1, numComps> &,
                                          const std::array<Ti2, numComps> &);

    template<typename To, typename Ti1, typename Ti2>
    static std::array<To, numComps> sum (const std::array<Ti1, numComps> &,
                                         const std::array<Ti2, numComps> &);

    static void add (std::array<unsigned, numComps> &,
                     const std::array<int, numComps> &),
                add (std::vector<int> &,
                     const std::vector<int> &),
                subtract (std::array<unsigned, numComps> &,
                          const std::array<int, numComps> &),
                subtract (std::vector<int> &,
                          const std::vector<int> &);

    template<typename T1, typename T2>
    static void addOrSubtract (std::array<T1, numComps> &,
                               const std::array<T2, numComps> &,
                               const int);

    template<typename T1, typename T2>
    static void addOrSubtract (std::vector<T1> &,
                               const std::vector<T2> &,
                               const int);

    static bool mismatchingPopDiff (const std::array<unsigned, numComps> &,
                                   const std::array<unsigned, numComps> &,
                                   const unsigned,
                                   const int);
    static bool mismatchingPopDiff (const std::array<unsigned, numComps> &,
                                    const std::array<unsigned, numComps> &,
                                    const std::array<int, numComps> &);


    static bool hasInvalidPop (const std::array<unsigned, numComps> &);

    static bool wouldHaveInvalidPop (const std::array<unsigned, numComps> &,
                                     const std::array<int, numComps> &,
                                     const int);
    static bool wouldHaveInvalidPop (const std::array<unsigned, numComps> &,
                                     const unsigned,
                                     const int,
                                     const int);


    static bool hasInvalidPopDiff (const std::array<unsigned, numComps> &,
                                   const std::array<unsigned, numComps> &);

    // static bool wouldHaveInvalidPopDiff (const std::array<unsigned, numComps> &,
    //                                      const std::array<unsigned, numComps> &
    //                                      const std::array<int, numComps> &,
    //                                      const int);
    // static bool wouldHaveInvalidPopDiff (const std::array<unsigned, numComps> &,
    //                                      const unsigned,
    //                                      const int,
    //                                      const int);

    bool wouldBreakFermiPrinc (std::array<unsigned, numComps> &,
                               std::array<int, numComps> &);

    void calculateAverageSegmentLengths (const double &,
                                         std::vector<std::vector<double> > &) const;

    void saveData (const unsigned long long &,
                   const unsigned long long &,
                   const double &,
                   const std::string &) const;
    void compressData (nlohmann::json &) const;
    void saveBinaryFiles (const unsigned long long &,
                          const unsigned long long &,
                          const std::string &) const;


    void saveConfiguration (const double &) const;


  friend class ofAppWorm;
  friend class DisplayWorm;
};

#endif