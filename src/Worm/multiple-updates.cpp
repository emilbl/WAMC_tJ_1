#include "Worm.h"

using namespace std;

void Worm::multipleUpdates (
  const unsigned        simulationTime_preSimulation,
  const unsigned        simulationTime,
  const unsigned long & numDataToAquire
) {

  // timers
  Tic tStart{};


  ////
  //// pre-simulation determining how frequently one should bin
  ////
  cout << settings::cout::enterBlue
       << "######################################" << endl
       << "#                                    #" << endl
       << "#   About to perform equilibration   #" << endl
       << "#                                    #" << endl
       << "######################################" << endl
       << settings::cout::resetStyle;

  // maximum number of Z configurations to save
  const unsigned long maxNumBins = 1000;

  // keeps track of how many Z configurations have been visited
  unsigned long numZcounter = 0;

  // initiate timer
  Tic t0{};

  // set beta
  this->setBeta(this->beta_target);

  do {
    for (unsigned long i = 0; i < this->numUpdatesPerTimeCheck; i++) {
      this->update();

      // try to bin every numZsPerSave-th time
      if (numZcounter < maxNumBins) this->trySampleData(this->beta_target);

      // increase the counter if the current configuration is a Z
      if ( ! this->numWorms) numZcounter++;
    }

  } while (t0.toc() < simulationTime_preSimulation);

  // if we have not managed to sample any partition configurations
  if ( ! numZcounter) {
    cout << settings::cout::enterRed
         << "######################################################" << endl
         << "#                                                    #" << endl
         << "#   No sampled Z-configurations during calibration   #" << endl
         << "#                                                    #" << endl
         << "######################################################" << endl
         << settings::cout::resetStyle;
  }

  ////
  //// estimate how often one should bin
  ////
  if (numDataToAquire == 0) {
    // save all data
    this->samplingPeriod = 1;

    printf("[%lu data points gathered in %.2fs. However the bin period is forced to %lu]\n",
           numZcounter,
           t0.toc(),
           this->samplingPeriod);
    fflush(stdout);
  } else {
    this->samplingPeriod = max(1., numZcounter * simulationTime / (numDataToAquire * t0.toc()));

    printf("[%lu data points gathered in %.2fs making the bin period %lu]\n",
           numZcounter,
           t0.toc(),
           this->samplingPeriod);
    fflush(stdout);
  }


  ////
  //// write pre-simulation data to file
  ////
  if (this->savePresimulationData) {
    this->saveData(0, this->currBinCount0, this->beta, "_ps");
  }

  ////
  //// reset bins and bin counter
  ////
  this->resetDataContainers();

  // also reset update counter
  this->currUpdateCount = 0;





  ////
  //// gather the actual data
  ////
  cout << settings::cout::enterBlue
       << "######################################" << endl
       << "#                                    #" << endl
       << "#   About to start main simulation   #" << endl
       << "#                                    #" << endl
       << "######################################" << endl
       << settings::cout::resetStyle;

  // reset timer
  tStart.reset();


  if (numDataToAquire) {
    ////
    //// if number of data points is specified
    ////
    unsigned numPrints = 100,
             numSaves  = 10;

    // more than one day
    if (simulationTime >= 1 * 24 * 60 * 60) {
      numPrints = 1000,
      numSaves  = 100;
    }

    unsigned long printPeriod = numDataToAquire / numPrints,
                  savePeriod = numDataToAquire / numSaves;


    unsigned long nextPrint = printPeriod,
                  nextSave = savePeriod,
                  prevSave = 0;

    Tic tPrint{};

    // set beta
    this->setBeta(this->beta_target);

    while (this->currBinCount0 < numDataToAquire) {
      // suggest updating the worm
      this->update();

      // try to bin
      this->trySampleData(this->beta_target);

      ////
      //// trySample increments currBinCount0
      ////

      // output status
      if (this->currBinCount0 == nextPrint) {
        // increment trigger
        nextPrint = min(nextPrint + printPeriod, numDataToAquire);

        printf("[%lu of %lu (%.2f%%) data points acquired in %.2fs (%.2f)]\n",
               this->currBinCount0,
               numDataToAquire,
               100 * this->currBinCount0 / (double) numDataToAquire,
               tStart.toc(),
               tPrint.toc());
        fflush(stdout);

        // reset clock
        tPrint.reset();
      }

      // write data to file
      if (this->currBinCount0 == nextSave) {

        // timer
        Tic tSave{};

        // write to file
        this->saveData(prevSave, this->currBinCount0, this->beta_target, "");

        // update
        prevSave = nextSave;

        // increment trigger
        nextSave = min(nextSave + savePeriod, numDataToAquire);

        printf("[%lu / %lu / %lu of %lu / . / . (%.2f%%) data points from %llu updates saved in %.2fs (%.2f)]\n",
               this->currBinCount0,
               this->currBinCount1,
               this->currBinCount2,
               numDataToAquire,
               100 * this->currBinCount0 / (double) numDataToAquire,
               this->currUpdateCount,
               tStart.toc(),
               tSave.toc());
        fflush(stdout);

        // reset print clock so as not to influence gathering time
        tPrint.reset();
      }
    }

  } else {
    ////
    //// if number of data points is not specified
    ////

    // save and print every 10 minutes
    const unsigned printAndSaveInterval = min((unsigned) (this->devMode ? 10 : 10 * 60), simulationTime);

    // set beta
    this->setBeta(this->beta_target);

    // to keep track how long we should simulate
    Tic mainSimulationTimer{};

    // to keep track of how often we should print and save
    Tic printAndSaveTimer{};

    // to keep track how many data points have already been saved
    unsigned long prevNumSavedDataPoints = 0;

    // to ensure that the data has been saved before a short simulation exits
    bool hasSaved = false;

    do {
      ////
      //// perform an amount of updates
      ////
      for (unsigned long i = 0; i < numUpdatesPerTimeCheck; i++) {
        // try update worm
        this->update();

        // try to bin
        this->trySampleData(this->beta_target);
      }


      ////
      //// try print and save
      ////
      if (printAndSaveTimer.toc() > printAndSaveInterval) {
        ////
        //// print
        ////
        printf("[%lu / %lu / %lu data points from %llu updates saved in %.2fs (%.2f)]\n",
               this->currBinCount0,
               this->currBinCount1,
               this->currBinCount2,
               this->currUpdateCount,
               mainSimulationTimer.toc(),
               printAndSaveTimer.toc());

        ////
        //// save
        ////
        this->saveData(prevNumSavedDataPoints, this->currBinCount0, this->beta_target, "");
        prevNumSavedDataPoints = this->currBinCount0;
        hasSaved = true;

        // reset timer
        printAndSaveTimer.reset();


        ////
        //// temp
        ////
        for (unsigned a = 0; a < numComps; a++) {
          cout << "wormCount[" << a << "] = " << this->wormCount[a] << endl;
        }
        cout << "avgParticleDiff[1] = " << this->avgParticleDiff[1] / this->avgParticleDiff_count[1] * this->beta / tMax << endl;

      }

    } while (mainSimulationTimer.toc() < (double) simulationTime || ! hasSaved);


  }


}