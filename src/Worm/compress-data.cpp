#include "Worm.h"

using namespace std;
using json = nlohmann::json;


////
//// allow Analytics::Compressed to be converted into a json
////
void to_json (
  nlohmann::json & j,
  const Analytics::Compressed & c
) {
  j = nlohmann::json{{"mean",     c.mean},
                     {"std",      c.std},
                     {"mean_std", c.mean_std},
                     {"corrLen",  c.corrLen}};
}


void Worm::compressData (
  json & statistics
) const {

  ////
  //// compressed data
  ////
  clock_t tStart = clock();
  bool printElapsedTime = false;


  ////
  //// convert sign vector to other types to be used in calculations
  ////
  vector<double> signs_small_double(this->signHist_small.size(), 1);
  if (this->sample_sign) {
    for (unsigned i = 0; i < signs_small_double.size(); i++) {
      signs_small_double[i] = this->signHist_small[i];
    }
  }


  ////
  //// sign
  ////
  if (this->sample_sign) {
    printElapsedTime = true;
    const auto sign = Analytics::compress<int>({this->signHist_small.begin(), 1},
                                               this->currBinCount0,
                                               "sign");
    statistics["sign"].push_back(sign);
  }


  ////
  //// kinetic energy
  ////
  if (this->sample_kinetEnergy) {
    printElapsedTime = true;
    statistics["E_Ks"] = {};
    for (unsigned a = 0; a < numComps; a++) {
      const auto E_K = Analytics::compress_signed<double>({this->kinetEnergyHist.begin() + a, numComps},
                                                          {signs_small_double.begin(), 1},
                                                          this->currBinCount0,
                                                          "kinetic energy (" + to_string(a) + ")");
      statistics["E_Ks"].push_back(E_K);
    }
  }


  ////
  //// exchange energy
  ////
  if (has_J && this->sample_exchaEnergy) {
    printElapsedTime = true;
    statistics["E_Js"] = {};
    for (unsigned ab = 0; ab < numExternalInters; ab++) {
      const auto E_J = Analytics::compress_signed<double>({this->exchaEnergyHist.begin() + ab, numExternalInters},
                                                          {signs_small_double.begin(), 1},
                                                          this->currBinCount0,
                                                          "exchange energy (" + to_string(ab) + ")");
      statistics["E_Js"].push_back(E_J);
    }
  }


  ////
  //// potential energy
  ////
  if (this->sample_potenEnergy) {
    printElapsedTime = true;
    statistics["E_mus"] = {};
    for (unsigned a = 0; a < numComps; a++) {
      const auto E_mu = Analytics::compress_signed<double>({this->potenEnergyHist.begin() + a, numComps},
                                                           {signs_small_double.begin(), 1},
                                                           this->currBinCount0,
                                                           "potential energy (" + to_string(a) + ")");
      statistics["E_mus"].push_back(E_mu);
    }
  }


  ////
  //// on-site interaction
  ////
  if (has_U && this->sample_interEnergy) {
    printElapsedTime = true;
    statistics["E_Us"] = {};
    for (unsigned a = 0; a < numComps; a++) {
      for (unsigned b = a; b < numComps; b++) {
        if (has_U) {
          const auto i_ab = Worm::i_ab(a, b);
          const auto E_U = Analytics::compress_signed<double>({this->interEnergyHist.begin() + i_ab, numInters},
                                                              {signs_small_double.begin(), 1},
                                                              this->currBinCount0,
                                                              "on-site interaction energy (" + to_string(a) + ", " + to_string(b) + ")");
          statistics["E_Us"].push_back(E_U);
        }
      }
    }
  }


  ////
  //// NN interaction energy
  ////
  if (has_U_nn && this->sample_interEnergy) {
    printElapsedTime = true;
    statistics["E_U_nn"] = {};
    const auto E_U_nn = Analytics::compress_signed<double>({this->nnInterEnergyHist.begin(), 1},
                                                           {signs_small_double.begin(), 1},
                                                           this->currBinCount0,
                                                           "nearest neighbor interaction energy");
    statistics["E_U_nn"].push_back(E_U_nn);
  }


  ////
  //// total energy
  ////
  if (this->sample_totalEnergy) {
    printElapsedTime = true;
    statistics["E"] = {};
    const auto E = Analytics::compress_signed<double>({this->totEnergyHist.begin(), 1},
                                                      {signs_small_double.begin(), 1},
                                                      this->currBinCount0,
                                                      "total energy");
    statistics["E"].push_back(E);
  }


  ////
  //// number of particles
  ////
  if (this->sample_numParticles) {
    printElapsedTime = true;
    statistics["Ns"] = {};
    for (unsigned a = 0; a < numComps; a++) {
      const auto N = Analytics::compress_signed<int>({this->numParticlesHist.begin() + a, numComps},
                                                     {this->signHist_small.begin(), 1},
                                                     this->currBinCount0,
                                                     "particle number (" + to_string(a) + ")");
      statistics["Ns"].push_back(N);
    }
  }


  ////
  //// winding number
  ////
  if (this->sample_numWinds) {
    printElapsedTime = true;
    statistics["Ws"] = {};
    for (unsigned a = 0; a < numComps; a++) {
      for (unsigned b = a; b < numComps; b++) {
        for (unsigned d = 0; d < this->numDims; d++) {
          const auto WW = Analytics::elemProd<int, int, int>({this->numWindsHist.begin() + a * this->numDims + d, numComps * this->numDims},
                                                             {this->numWindsHist.begin() + b * this->numDims + d, numComps * this->numDims},
                                                             this->currBinCount0);
          const auto W = Analytics::compress_signed<int>({WW.begin(), 1},
                                                         {this->signHist_small.begin(), 1},
                                                         this->currBinCount0,
                                                         "winding number (" + to_string(a) + ", " + to_string(b) +  + ", " + to_string(d) + ")");
          statistics["Ws"].push_back(W);
        }
      }
    }
  }



  ////
  //// calculate the drag interaction varrho_ab and L times rho_{tot, SCF}
  ////
  if (numComps > 1 && this->sample_numWinds) {
    printElapsedTime = true;
    ////
    //// varrho description
    ////
    auto varrho = [&](
      const vector<Analytics::DataSet<int> > & dataSets,
      const vector<int>::size_type & from,
      const vector<int>::size_type & to
    ) {
      //the data is assumed to be in the order Wax_1, wax_2, ..., Wbx_1, Wbx_2, ...

      double nominator = 0,
             denominator_a = 0,
             denominator_b = 0;

      for (unsigned d = 0; d < this->numDims; d++) {
        nominator += Analytics::statProdAvg<int>(dataSets[d],
                                                 dataSets[this->numDims + d],
                                                 from,
                                                 to);

        denominator_a += Analytics::statProdAvg<int>(dataSets[d],
                                                     dataSets[d],
                                                     from,
                                                     to);

        denominator_b += Analytics::statProdAvg<int>(dataSets[this->numDims + d],
                                                     dataSets[this->numDims + d],
                                                     from,
                                                     to);
      }

      return nominator / sqrt(denominator_a * denominator_b);
    };

    ////
    //// rho_tot and rho_SCF description
    ////
    auto Lrho_tot = [&](
      const vector<Analytics::DataSet<int> > & dataSets,
      const vector<int>::size_type & from,
      const vector<int>::size_type & to
    ) {
      //the data is assumed to be in the order Wax_1, wax_2, ..., Wbx_1, Wbx_2, ...

      double Waa = 0,
             Wbb = 0,
             Wab = 0;

      for (unsigned d = 0; d < this->numDims; d++) {
        Wab += Analytics::statProdAvg<int>(dataSets[d],
                                           dataSets[this->numDims + d],
                                           from,
                                           to);

        Waa += Analytics::statProdAvg<int>(dataSets[d],
                                           dataSets[d],
                                           from,
                                           to);

        Wbb += Analytics::statProdAvg<int>(dataSets[this->numDims + d],
                                           dataSets[this->numDims + d],
                                           from,
                                           to);
      }

      return this->L * (Waa + Wbb + 2 * Wab) / (2 * beta);
    };
    auto Lrho_SCF = [&](
      const vector<Analytics::DataSet<int> > & dataSets,
      const vector<int>::size_type & from,
      const vector<int>::size_type & to
    ) {
      //the data is assumed to be in the order Wax_1, wax_2, ..., Wbx_1, Wbx_2, ...

      double Waa = 0,
             Wbb = 0,
             Wab = 0;

      for (unsigned d = 0; d < this->numDims; d++) {
        Wab += Analytics::statProdAvg<int>(dataSets[d],
                                           dataSets[this->numDims + d],
                                           from,
                                           to);

        Waa += Analytics::statProdAvg<int>(dataSets[d],
                                           dataSets[d],
                                           from,
                                           to);

        Wbb += Analytics::statProdAvg<int>(dataSets[this->numDims + d],
                                           dataSets[this->numDims + d],
                                           from,
                                           to);
      }

      return this->L * (Waa + Wbb - 2 * Wab) / (2 * beta);
    };




    ////
    //// loads data container: Wa_x, Wa_y, Wa_z, ..., Wb_x, Wb_y, Wb_z, ..., Wc_x, Wc_y, Wc_z, ...
    ////
    vector<Analytics::DataSet<int> > dataSets = {};
    for (unsigned a = 0; a < numComps; a++) {
      for (unsigned d = 0; d < this->numDims; d++) {
        dataSets.push_back({this->numWindsHist.cbegin() + a * this->numDims + d, numComps * this->numDims});
      }
    }

    statistics["varrhos"]   = vector<json>{};
    statistics["Lrho_tots"] = vector<json>{};
    statistics["Lrho_SCFs"] = vector<json>{};

    for (unsigned a = 0; a < numComps; a++) {
      for (unsigned b = a + 1; b < numComps; b++) {

        // create a subset of data: Wa_x, Wa_y, Wa_z, ..., Wb_x, Wb_y, Wb_z, ...
        vector<Analytics::DataSet<int> > subDataSets = {};
        for (unsigned d = 0; d < this->numDims; d++) {
          subDataSets.push_back(dataSets[a * numDims + d]);
        }
        for (unsigned d = 0; d < this->numDims; d++) {
          subDataSets.push_back(dataSets[b * numDims + d]);
        }


        ////
        //// calculate mean value
        ////
        const auto varrho_mean = varrho(subDataSets, 0, this->currBinCount0);
        const auto Lrho_tot_mean = Lrho_tot(subDataSets, 0, this->currBinCount0);
        const auto Lrho_SCF_mean = Lrho_SCF(subDataSets, 0, this->currBinCount0);

        ////
        //// estimate errors
        ////
        const auto varrho_mean_std = Analytics::bootstrap<int>(varrho,
                                                               subDataSets,
                                                               this->currBinCount0,
                                                               varrho_mean,
                                                               "varrho");
        const auto Lrho_tot_mean_std = Analytics::bootstrap<int>(Lrho_tot,
                                                                 subDataSets,
                                                                 this->currBinCount0,
                                                                 Lrho_tot_mean,
                                                                 "Lrho_tot");
        const auto Lrho_SCF_mean_std = Analytics::bootstrap<int>(Lrho_SCF,
                                                                 subDataSets,
                                                                 this->currBinCount0,
                                                                 Lrho_SCF_mean,
                                                                 "Lrho_SCF");

        ////
        //// store data in json object
        ////
        statistics["varrhos"].push_back({{"mean", varrho_mean},
                                         {"mean_std", varrho_mean_std}});
        statistics["Lrho_tots"].push_back({{"mean", Lrho_tot_mean},
                                           {"mean_std", Lrho_tot_mean_std}});
        statistics["Lrho_SCFs"].push_back({{"mean", Lrho_SCF_mean},
                                           {"mean_std", Lrho_SCF_mean_std}});
      }
    }
  }

  if (printElapsedTime) {
    printf("[Data compressed in %.2fs]\n",
           (double)(clock() - tStart)/CLOCKS_PER_SEC);
    fflush(stdout);
  }
}