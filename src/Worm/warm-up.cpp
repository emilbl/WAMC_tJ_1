#include "Worm.h"

using namespace std;
using json = nlohmann::json;

void Worm::warmUp (
  const unsigned simulationTime_annealing,
  const double & beta_beg,
  const unsigned simulationTime_optimizeMu
) {

  // timers
  Tic tStart{};


  ////
  //// create initial mu intervals
  ////
  const auto mu_ia = this->H.getMu();
  for (unsigned a = 0; a < numComps; a++) {
    if (this->muOptimization[a]) {

      double mu;
      if (mu_ia.size() == numComps) {
        mu = mu_ia[a];
      } else if (mu_ia.size() == 1) {
        mu = mu_ia[0];
      } else {
        cout << "Worm::warmUp: implement me" << endl;
        exit(EXIT_SUCCESS);
      }

      this->muIntervals[2 * a + 1] = (mu > 0 ? 1.1 : 0.9) * mu + 1e-2;  // upper
      this->muIntervals[2 * a + 0] = (mu > 0 ? 0.9 : 1.1) * mu - 1e-2;  // lower
    }
  }


  if (simulationTime_annealing) {

    cout << settings::cout::enterBlue
         << "############################################" << endl
         << "#                                          #" << endl
         << "#   About to perform simulated annealing   #" << endl
         << "#                                          #" << endl
         << "############################################" << endl
         << settings::cout::resetStyle;

    // number of intermediate betas
    const unsigned numIntermediateBetas = 9;

    ////
    //// maximum number of Z configurations to save
    ////
    const unsigned long maxNumBins = 10000;
    const unsigned samplingPeriod = 50;

    ////
    //// how often one should save to file
    ////
    const unsigned long maxNumBinsPerBeta = max(maxNumBins / this->beta_target, 1.);

    ////
    //// containers for statistics regarding the simulated annealing equilibration
    ////
    vector<double>             betas(numIntermediateBetas + 1, 0);
    vector<double>             numSecs(numIntermediateBetas + 1, 0);
    vector<double>             musHist((numIntermediateBetas + 1) * numComps, 0);
    vector<double>             muIntervalsHist((numIntermediateBetas + 1) * numComps * 2, 0);
    vector<unsigned long>      numZcounts(numIntermediateBetas + 1, 0);
    vector<unsigned long long> numGsCounts((numIntermediateBetas + 1) * 2 * numComps, 0);

    const unsigned numIntervalReductions = 1;

    for (unsigned betaIndex = 0; betaIndex < numIntermediateBetas + 1; betaIndex++) {
      // compute beta and the simulation time
      const double deltaBeta = (this->beta_target - beta_beg) / (numIntermediateBetas + 1.);
      const double _beta = beta_beg + deltaBeta * (betaIndex + 1.);
      const double simulationTime = simulationTime_annealing / (numIntermediateBetas + 1.)
                                  * (beta_beg + deltaBeta * (betaIndex + 1.))
                                  / (beta_beg + 0.5 * deltaBeta * (numIntermediateBetas + 2.));

      // set current beta
      this->setBeta(_beta);


      // start times
      Tic t_single_beta{};

      // keeps track of how many Z configurations have been visited
      unsigned long numZcounter = 0,
                    numBinnedZcounter = 0;


      // perform binary interval reduction to find optimum mu
      if (this->H.model.hasEquivalentComponents()) {
        if (this->muOptimization[0]) {
          this->searchOptimalMu(_beta,
                                0,
                                this->muIntervals[0],
                                this->muIntervals[1],
                                simulationTime,
                                numUpdatesPerTimeCheck,
                                numIntervalReductions,
                                numZcounter,
                                numBinnedZcounter,
                                maxNumBinsPerBeta,
                                samplingPeriod);
        } else {
          // simulate the specified amount of time
          Tic tic{};
          do {
            for (long unsigned i = 0; i < numUpdatesPerTimeCheck; i++) {
              this->update();

              // try to bin every "samplingPeriod"-th Z configurations but no more than "maxNumBinsPerBeta" times
              if (
                ! this->numWorms &&
                ! (numZcounter % samplingPeriod) &&
                numBinnedZcounter < maxNumBinsPerBeta
              ) {
                this->trySampleData(beta);
                numBinnedZcounter++;
              }

              // increase the counter if the current configuration is a Z
              if ( ! this->numWorms) numZcounter++;
            }
          } while (tic.toc() < simulationTime);
        }
      } else {
        // loop over all components which are in the canonical ensemble
        for (unsigned a = 0; a < numComps; a++) {
          if (this->muOptimization[a]) {
            this->searchOptimalMu(_beta,
                                  a,
                                  this->muIntervals[2 * a + 0],
                                  this->muIntervals[2 * a + 1],
                                  simulationTime,
                                  numUpdatesPerTimeCheck,
                                  numIntervalReductions,
                                  numZcounter,
                                  numBinnedZcounter,
                                  maxNumBinsPerBeta,
                                  samplingPeriod);
          }
        }
      }



      ////
      //// store statistics
      ////
      const unsigned i = betaIndex;
      betas[i] = _beta;
      numZcounts[i] = numZcounter;
      numSecs[i] = t_single_beta.toc();
      const auto mus = this->H.getChemicalPotential();
      copy(mus.begin(),
           mus.end(),
           musHist.begin() + i * numComps);
      copy(this->muIntervals.begin(),
           this->muIntervals.end(),
           muIntervalsHist.begin() + i * numComps * 2);


      ////
      //// output to terminal
      ////
      printf("[Simulated annealing: β=%.2f, #Z-configs=%lu completed in %.2fs (%.2f)]\n",
             _beta,
             numZcounter,
             tStart.toc(),
             t_single_beta.toc());
      fflush(stdout);

      // reset
      t_single_beta.reset();
    }

    // store in warm up statistics json object
    this->warmupStatistics["betas"]           = betas;
    this->warmupStatistics["numZcounts"]      = numZcounts;
    this->warmupStatistics["numSecs"]         = numSecs;
    this->warmupStatistics["musHist"]         = musHist;
    this->warmupStatistics["muIntervalsHist"] = muIntervalsHist;

    ////
    //// write warm up statistics to file
    ////
    json statistics;
    statistics["warmup"] = this->warmupStatistics;

    // prettified json string
    stringstream ss;
    ss << setw(2) << statistics << endl;

    // write the json object to file
    this->writeToFile.aString(ss.str(), "statistics.json");



    ////
    //// write simulated annealing data to file
    ////
    if (this->saveWarmUpData) {
      this->saveBinaryFiles(0, this->currBinCount0, "_sa");
    }


    ////
    //// reset bins and bin counter
    ////
    this->resetDataContainers();
  }



  if (simulationTime_optimizeMu) {
    tStart.reset();

    cout << settings::cout::enterBlue
         << "############################################" << endl
         << "#                                          #" << endl
         << "#   About to optimize chemical potential   #" << endl
         << "#                                          #" << endl
         << "############################################" << endl
         << settings::cout::resetStyle;

    const unsigned numIntervalReductions_final = 10;
    // perform binary interval reduction to find optimum mu
    if (this->H.model.hasEquivalentComponents()) {
      if (this->muOptimization[0]) {
        this->searchOptimalMu(this->beta_target,
                              0,
                              this->muIntervals[0],
                              this->muIntervals[1],
                              simulationTime_optimizeMu,
                              numUpdatesPerTimeCheck,
                              numIntervalReductions_final);
      }
    } else {
      // loop over all components which are in the canonical ensemble
      for (unsigned a = 0; a < numComps; a++) {
        if (this->muOptimization[a]) {
          this->searchOptimalMu(this->beta_target,
                                a,
                                this->muIntervals[2 * a + 0],
                                this->muIntervals[2 * a + 1],
                                simulationTime_optimizeMu,
                                numUpdatesPerTimeCheck,
                                numIntervalReductions_final);
        }
      }
    }


    ////
    //// output to terminal
    ////
    printf("[search optimal mu: completed in %.2fs ]\n", tStart.toc());
    fflush(stdout);


    ////
    //// reset bins and bin counter
    ////
    this->resetDataContainers();
  }
}