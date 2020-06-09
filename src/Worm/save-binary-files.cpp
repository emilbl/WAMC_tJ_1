#include "Worm.h"

using namespace std;


void Worm::saveBinaryFiles (
  const unsigned long long & prevNumData,
  const unsigned long long & newNumData,
  const string &             postfix
) const {

  const bool append = !! prevNumData;

  ////
  //// average sign
  ////
  const double avgSign = ((double) this->posSignCount - (double) this->negSignCount)
                       / (this->posSignCount + this->negSignCount);

  ////
  //// small data
  ////
  if (this->sample_sign) {
    this->writeToFile.template aVector<int>(this->signHist_small.cbegin() + prevNumData,
                                            vector<int>::size_type((newNumData - prevNumData)),
                                            "sign_small" + postfix + ".bin",
                                            append);
  }


  if (this->sample_numParticles) {
    this->writeToFile.template aVector<int>(this->numParticlesHist.cbegin() + prevNumData * numComps,
                                            vector<int>::size_type((newNumData - prevNumData) * numComps),
                                            "numParticles" + postfix + ".bin",
                                            append);
  }

  if (this->sample_numWinds) {
    this->writeToFile.template aVector<int>(this->numWindsHist.cbegin() + prevNumData * numComps * this->numDims,
                                            vector<int>::size_type((newNumData - prevNumData) * numComps * this->numDims),
                                            "numWinds" + postfix + ".bin",
                                            append);
  }

  if (this->sample_kinetEnergy) {
    this->writeToFile.template aVector<double>(this->kinetEnergyHist.cbegin() + prevNumData * numComps,
                                               vector<double>::size_type((newNumData - prevNumData) * numComps),
                                               "kinetEnergy" + postfix + ".bin",
                                               append);
  }

  if (this->sample_potenEnergy) {
    this->writeToFile.template aVector<double>(this->potenEnergyHist.cbegin() + prevNumData * numComps,
                                               vector<double>::size_type((newNumData - prevNumData) * numComps),
                                               "potenEnergy" + postfix + ".bin",
                                               append);
  }

  if (this->sample_totEnergy) {
    this->writeToFile.template aVector<double>(this->totEnergyHist.cbegin() + prevNumData,
                                               vector<double>::size_type((newNumData - prevNumData)),
                                               "totalEnergy" + postfix + ".bin",
                                               append);
  }

  if (has_J) {
    if (this->sample_exchaEnergy) {
      this->writeToFile.template aVector<double>(this->exchaEnergyHist.cbegin() + prevNumData * numExternalInters,
                                                 vector<double>::size_type((newNumData - prevNumData) * numExternalInters),
                                                 "exchaEnergy" + postfix + ".bin",
                                                 append);
    }
  }

  if (has_U) {
    if (this->sample_interEnergy) {
      this->writeToFile.template aVector<double>(this->interEnergyHist.cbegin() + prevNumData * numInters,
                                                 vector<double>::size_type((newNumData - prevNumData) * numInters),
                                                 "interEnergy" + postfix + ".bin",
                                                 append);
    }
  }

  if (has_U_nn) {
    if (this->sample_nnInterEnergy) {
      this->writeToFile.template aVector<double>(this->nnInterEnergyHist.cbegin() + prevNumData,
                                                 vector<double>::size_type((newNumData - prevNumData)),
                                                 "nnInterEnergy" + postfix + ".bin",
                                                 append);
    }
  }


  ////
  //// large data
  ////
  const unsigned long prevNumLargeData = prevNumData / this->largeDataSamplingPeriod,
                      newNumLargeData  = newNumData  / this->largeDataSamplingPeriod;

  if (this->sample_sign) {
    this->writeToFile.template aVector<int>(this->signHist_large.cbegin() + prevNumLargeData,
                                            vector<int>::size_type(newNumLargeData - prevNumLargeData),
                                            "sign_large" + postfix + ".bin",
                                            append);
  }

  if(this->sample_numParticlesAtSite) {
    this->writeToFile.template aVector<double>(this->numParticlesAtSiteHist.cbegin() + prevNumLargeData * this->numSites * numComps,
                                               vector<double>::size_type((newNumLargeData - prevNumLargeData) * this->numSites * numComps),
                                               "numParticlesAtSite" + postfix + ".bin",
                                               append);
  }

  if(this->sample_flow) {
    this->writeToFile.template aVector<double>(this->flowHist.cbegin() + prevNumLargeData * this->numSites * numComps * this->numDims,
                                               vector<double>::size_type((newNumLargeData - prevNumLargeData) * this->numSites * numComps * this->numDims),
                                               "flow" + postfix + ".bin",
                                               append);
  }

  if(this->sample_instantParticleNum) {
    this->writeToFile.template aVector<unsigned>(this->instantParticleNumHist.cbegin() + prevNumLargeData * this->numSites * numComps,
                                                 vector<unsigned>::size_type((newNumLargeData - prevNumLargeData) * this->numSites * numComps),
                                                 "instantParticleNum" + postfix + ".bin",
                                                 append);
  }

  if(this->sample_instantParticleNumConv) {
    this->writeToFile.template aVector<unsigned>(this->instantParticleNumConvHist.cbegin() + prevNumLargeData * this->numSites * numInters,
                                                 vector<unsigned>::size_type((newNumLargeData - prevNumLargeData) * this->numSites * numInters),
                                                 "instantParticleNumConv" + postfix + ".bin",
                                                 append);
  }


  ////
  //// C1, C2 and C3
  ////
  const auto numLargeData = newNumLargeData;
  if (modelType == tJ && this->compute_C1 && this->Ns[1] == 2) {
    vector<double> avgC1_normalized(this->avgC1.size());
    transform(this->avgC1.begin(),
              this->avgC1.end(),
              avgC1_normalized.begin(),
              [&] (const double & x) { return x / numLargeData / avgSign; });
    this->writeToFile.aVector<double>(avgC1_normalized.begin(),
                                      avgC1_normalized.size(),
                                      "C1" + postfix + ".bin",
                                      false);
  }
  if (modelType == tJ && this->compute_C2 && this->Ns[1] == 2) {
    vector<double> avgC2_normalized(this->avgC2.size());
    transform(this->avgC2.begin(),
              this->avgC2.end(),
              avgC2_normalized.begin(),
              [&] (const double & x) { return x / numLargeData / avgSign; });
    this->writeToFile.aVector<double>(avgC2_normalized.begin(),
                                      avgC2_normalized.size(),
                                      "C2" + postfix + ".bin",
                                      false);
  }
  if (modelType == tJ && this->compute_C3 && this->Ns[1] == 2) {
    vector<double> avgC3_normalized(this->avgC3.size());
    transform(this->avgC3.begin(),
              this->avgC3.end(),
              avgC3_normalized.begin(),
              [&] (const double & x) { return x / numLargeData / avgSign; });
    this->writeToFile.aVector<double>(avgC3_normalized.begin(),
                                      avgC3_normalized.size(),
                                      "C3" + postfix + ".bin",
                                      false);
  }


  ////
  //// Green's function histogram
  ////
  if (this->sample_GreensFunction && this->isTranslInvariant) {

    // first compute the partition function
    const double Z      = this->Z_bin; // / (double) this->Z_bin_count;
    const double Z_bose = this->Z_bin_count; // / (double) this->Z_bin_count;

    unsigned G_size = this->G_hist[0].size();
    vector<double> G(numComps * G_size);
    vector<double> G_bose(numComps * G_size);

    for (unsigned a = 0; a < numComps; a++) {

      const double dt = (this->isCanonical[a] ? 2 * this->maxNsDiff[a] : this->beta)
                      / (double) G_hist_N;


      const double etaFactor = this->H.getDicsoSquredFact(a, this->beta);


      ////
      //// Since we are averaging out a few degrees of freedom we must "normalize" accordingly, i.e. <x> = 1/N \sum_{i=1}^N x_i
      ////
      ////   G(t_2 - t_1 = \Delta) = 1/\beta \int_0^\beta dt G(t_2 = \Delta + t, t_1 = t)
      ////
      const double dofFactor = this->numSites   // translation invariance
                             * this->beta;      // invariant under overall shift in tau


      const double normalization = Z * dt * etaFactor * dofFactor;
      transform(this->G_hist[a].begin(),
                this->G_hist[a].end(),
                G.begin() + a * G_size,
                [&] (const long & num) { return num / normalization; });

      const double normalization_bose = Z_bose * dt * etaFactor * dofFactor;
      transform(this->G_hist_count[a].begin(),
                this->G_hist_count[a].end(),
                G_bose.begin() + a * G_size,
                [&] (const long & num) { return num / normalization_bose; });

      // cout << "count of G[" << a << "]: " << accumulate(this->G_hist_count[a].begin(), this->G_hist_count[a].end(), 0l) << endl;
    }

    this->writeToFile.template aVector<double>(G.cbegin(), G.size(), "G" + postfix + ".bin", false);
    this->writeToFile.template aVector<double>(G_bose.cbegin(), G_bose.size(), "G_bose" + postfix + ".bin", false);
  }



  ////
  //// four point correlator
  ////
  if (this->sample_fourPointCorrelator) {

    // first compute the partition function
    const double Z      = this->Z_bin;
    const double Z_bose = this->Z_bin_count;

    vector<double> fpc;
    vector<double> fpc_bose;

    for (unsigned a = 0; a < numComps; a++) {
      // if we are sampling for the particular component
      if (this->fpc_hist[a].size()) {

        const double dT = 0.5 * this->beta / (double) fpc_hist_N;
        const double dt = this->maxNsDiff[a] / (double) fpc_hist_N;

        const double etaFactor = pow(this->H.getDicsoSquredFact(a, this->beta), 2.0);

        const double dofFactor = this->numSites   // translation invariance
                               * this->numSites   // such that at fpc(T=0, t1=0, t2=0) = 1
                               * this->beta;      // invariant under overall shift in tau



        // reserve so you don't have grow the vector geometrically.
        fpc.reserve(fpc.size() + this->fpc_hist[a].size());
        fpc_bose.reserve(fpc_bose.size() + this->fpc_hist[a].size());

        const double normalization = Z * dT * dt * dt * etaFactor * dofFactor;
        transform(this->fpc_hist[a].begin(),
                  this->fpc_hist[a].end(),
                  back_inserter(fpc),
                  [&] (const long & num) {
                    // obtain index of element
                    size_t I = &num - &fpc_hist[a][0];

                    // compute i0, i1, i2
                    const unsigned i0 = I / fpc_hist_N / fpc_hist_N;
                    const unsigned i1 = (I / fpc_hist_N - i0 * fpc_hist_N);
                    const unsigned i2 = (I - i0 * fpc_hist_N * fpc_hist_N - i1 * fpc_hist_N);

                    // additional normalization
                    const double addNorm = (i1 != i2 ? 2. : 1.)                      // <- t1 ≥ t2 ordering -> twice the probability mass
                                         * (i1 + i2 == fpc_hist_N - 1 ? 0.5 : 1.);   // <- elements on the antidiagonal are half the size
                                                                                     //    (no contribution from t1 + t2 > maxNsDiff)

                    return num / (normalization * addNorm);
                  });


        const double normalization_bose = Z_bose * dT * dt * dt * etaFactor * dofFactor;
        transform(this->fpc_hist_count[a].begin(),
                  this->fpc_hist_count[a].end(),
                  back_inserter(fpc_bose),
                  [&] (const unsigned long & num) {
                    // obtain index of element
                    size_t I = &num - &fpc_hist_count[a][0];

                    // compute i0, i1, i2
                    const unsigned i0 = I / fpc_hist_N / fpc_hist_N;
                    const unsigned i1 = (I / fpc_hist_N - i0 * fpc_hist_N);
                    const unsigned i2 = (I - i0 * fpc_hist_N * fpc_hist_N - i1 * fpc_hist_N);

                    // additional normalization
                    const double addNorm = (i1 != i2 ? 2. : 1.)                      // <- t1 ≥ t2 ordering -> twice the probability mass
                                         * (i1 + i2 == fpc_hist_N - 1 ? 0.5 : 1.);   // <- elements on the antidiagonal are half the size
                                                                                     //    (no contribution from t1 + t2 > maxNsDiff)

                    return num / (normalization_bose * addNorm);
                  });
      }

      this->writeToFile.template aVector<double>(fpc.cbegin(), fpc.size(), "fpc" + postfix + ".bin", false);
      this->writeToFile.template aVector<double>(fpc_bose.cbegin(), fpc_bose.size(), "fpc_bose" + postfix + ".bin", false);
    }


  }


  ////
  //// density-density correlator
  ////
  if (this->sample_densityDensityCorrelator) {

    // first compute the partition function
    const double Z      = this->Z_bin;
    const double Z_bose = this->Z_bin_count;

    vector<double> ddc;
    vector<double> ddc_bose;

    for (unsigned a = 0; a < numComps; a++) {
      // if we are sampling for the particular component
      if (this->ddc_hist[a].size()) {

        const double dt_o = 0.5 * this->beta   / (double) ddc_hist_N_o;
        const double dt_d = this->maxNsDiff[a] / (double) ddc_hist_N_d;

        const double etaFactor = pow(this->H.getDicsoSquredFact(a, this->beta), 2.0);

        const double dofFactor = this->numSites   // translation invariance
                               * this->beta;      // invariant under overall shift in tau

        // reserve so you don't have grow the vector geometrically.
        ddc.reserve(     ddc.size()      + this->ddc_hist[a].size());
        ddc_bose.reserve(ddc_bose.size() + this->ddc_hist[a].size());

        const double normalization = Z * dt_o * dt_d * etaFactor * dofFactor;
        transform(this->ddc_hist[a].begin(),
                  this->ddc_hist[a].end(),
                  back_inserter(ddc),
                  [&] (const long & num) { return num; /* / normalization; */ });

        const double normalization_bose = Z_bose * dt_o * dt_d * etaFactor * dofFactor;
        transform(this->ddc_hist_count[a].begin(),
                  this->ddc_hist_count[a].end(),
                  back_inserter(ddc_bose),
                  [&] (const unsigned long & num) { return num; /* / normalization_bose; */ });
      }

      this->writeToFile.template aVector<double>(ddc.cbegin(), ddc.size(), "ddc" + postfix + ".bin", false);
      this->writeToFile.template aVector<double>(ddc_bose.cbegin(), ddc_bose.size(), "ddc_bose" + postfix + ".bin", false);
    }


  }



}