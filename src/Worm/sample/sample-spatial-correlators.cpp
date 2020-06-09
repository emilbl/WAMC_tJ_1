#include "../Worm.h"

using namespace std;


static inline
void innerProd (
  const vector<unsigned> & V1,
  const vector<unsigned> & V2,
  unsigned long *          res
) {
  /* SLOWER */ // unsigned val[numComps] = {};
  /* SLOWER */ // for (unsigned i = 0; i < V1.size(); i++) val[i % numComps] += V1[i] + V2[i];
  /* SLOWER */ // for (unsigned a = 0; a < numComps; a++) *(res + a) += factor * val[a];

  for (unsigned i = 0; i < V1.size(); i++) *(res + i % numComps) += V1[i] * V2[i];
}

static inline
void innerProd_interComponent (
  const vector<unsigned> & V1,
  const vector<unsigned> & V2,
  unsigned long *          res
) {

  for (unsigned i = 0; i < V1.size() / 2; i++) *res += V1[2*i] * V2[2*i + 1];
}

static inline
void convolute (
  const vector<unsigned> & M,
  const vector<unsigned> & L,
  unsigned long *          container
) {
  // PROFILING
  // profiler("convolute");

  if (L.size() > 2) {
    cout << "Worm::convolute: implement dimensions > 2  ->  EXIT" << endl;
    exit(EXIT_SUCCESS);
  }


  ////
  //// the other layer to be "slided"
  ////
  auto _M = M;


  ////
  //// roll to starting position by rolling (L_i + 1) / 2 back in each direction
  ////
  // major
  rotate(_M.begin(),
         _M.begin() + (L[0] + 1) / 2 * (L.size() == 2 ? L[1] : 1) * numComps,
         _M.end());
  // minor
  if (L.size() == 2) {
    for (unsigned i = 0; i < L[0]; i++) {
      rotate(_M.begin() + numComps * ( i      * L[1]),
             _M.begin() + numComps * ( i      * L[1] + (L[1] + 1) / 2),
             _M.begin() + numComps * ((i + 1) * L[1]));
    }
  }



  ////
  //// "slide back"
  ////
  if (L.size() == 2) {
    for (unsigned i = 0; i < L[0]; i++) {
      for (unsigned j = 0; j < L[1]; j++) {
        // perform inner product
        innerProd(M, _M, container + (i * L[1]+ j) * (numComps + 1));

        // inter component inner product
        innerProd_interComponent(M, _M, container + (i * L[1]+ j) * (numComps + 1) + 2);

        // rotate minor axis
        for (unsigned i = 0; i < L[0]; i++) {
          rotate(_M.begin() + numComps * ( i      * L[1]),
                 _M.begin() + numComps * ((i + 1) * L[1] - 1),
                 _M.begin() + numComps * ((i + 1) * L[1]));
        }
      }

      // rotate major axis
      rotate(_M.begin(), _M.end() - L[1] * numComps, _M.end());
    }
  } else {
    for (unsigned i = 0; i < L[0]; i++) {
      // perform inner product
      innerProd(M, _M, container + i * (numComps + 1));

      // inter component inner product
      innerProd_interComponent(M, _M, container + i * (numComps + 1) + 2);

      // rotate major axis
      rotate(_M.begin(), _M.end() -numComps, _M.end());
    }
  }

}


void Worm::sampleSpatialCorrelators_store (
  const unsigned shift,
  const unsigned holeIndexDiff,
  const double   dt,
  const bool     performConvolution
) {
  ////
  //// C1, C2, C3: append to outer container
  ////
  if (this->Ns[1] == 1) {
    if (this->compute_C1) {
      transform(this->_instC1.begin(),
                this->_instC1.end(),
                this->avgC1.begin() + shift * this->_instC1.size(),
                this->avgC1.begin() + shift * this->_instC1.size(),
                [&] (int correlation, double val) { return val + 2 * (correlation - 0.5) * dt; });   // "- 0.5" in order for the noise to vanish
                                                                                                     // "* 2" in order for the correlation to lie in between [-1, 1]
    }
    if (this->compute_C2) {
      transform(this->_instC2.begin(),
                this->_instC2.end(),
                this->avgC2.begin() + shift * this->_instC2.size(),
                this->avgC2.begin() + shift * this->_instC2.size(),
                [&] (int correlation, double val) { return val + 2 * (correlation - 0.5) * dt; });   // "- 0.5" in order for the noise to vanish
                                                                                                     // "* 2" in order for the correlation to lie in between [-1, 1]
    }
    if (this->compute_C3) {
      transform(this->_instC3.begin(),
                this->_instC3.end(),
                this->avgC3.begin() + shift * this->_instC3.size(),
                this->avgC3.begin() + shift * this->_instC3.size(),
                [&] (int correlation, double val) { return val + 2 * (correlation - 0.5) * dt; });   // "- 0.5" in order for the noise to vanish
                                                                                                     // "* 2" in order for the correlation to lie in between [-1, 1]
    }
  } else if (this->Ns[1] == 2) {
    const int sign = (shift == 0 || shift == 3) ? 1 : -1;
    if (this->compute_C1) {
      transform(this->_instC1.begin(),
                this->_instC1.end(),
                this->avgC1.begin() + holeIndexDiff * this->_instC1.size(),
                this->avgC1.begin() + holeIndexDiff * this->_instC1.size(),
                [&] (int correlation, double val) { return val + 2 * (correlation - 0.5) * dt * sign; });   // "- 0.5" in order for the noise to vanish
                                                                                                            // "* 2" in order for the correlation to lie in between [-1, 1]
    }
    if (this->compute_C2) {
      transform(this->_instC2.begin(),
                this->_instC2.end(),
                this->avgC2.begin() + holeIndexDiff * this->_instC2.size(),
                this->avgC2.begin() + holeIndexDiff * this->_instC2.size(),
                [&] (int correlation, double val) { return val + 2 * (correlation - 0.5) * dt * sign; });   // "- 0.5" in order for the noise to vanish
                                                                                                            // "* 2" in order for the correlation to lie in between [-1, 1]
    }
    if (this->compute_C3) {
      transform(this->_instC3.begin(),
                this->_instC3.end(),
                this->avgC3.begin() + holeIndexDiff * this->_instC3.size(),
                this->avgC3.begin() + holeIndexDiff * this->_instC3.size(),
                [&] (int correlation, double val) { return val + 2 * (correlation - 0.5) * dt * sign; });   // "- 0.5" in order for the noise to vanish
                                                                                                            // "* 2" in order for the correlation to lie in between [-1, 1]
    }
  }


  ////
  //// instantParticleNums: append to outer container
  ////
  if (this->compute_particleCorr) {
    transform(this->_instantParticleNums.begin(),
              this->_instantParticleNums.end(),
              this->avgParticleCorr.begin() + shift * this->_instantParticleNums.size(),
              this->avgParticleCorr.begin() + shift * this->_instantParticleNums.size(),
              [&] (unsigned particleNum, double val) { return val + particleNum * dt; });
  }


  ////
  //// perform convolution
  ////
  if (this->compute_particleConv && performConvolution) {
    convolute(this->instantParticleNums,
              this->Ls,
              &this->avgParticleConv[shift * this->avgParticleConv.size() / 4]);
  }
}



void Worm::sampleSpatialCorrelators_computeAndStore (
  const unsigned           shift,
  const vector<unsigned> & holeIs,
  const unsigned           holeIndexDiff,
  const double             dt,
  const bool               performConvolution
) {
  ////
  //// the data containers:
  //// instC1, instC2, instC3, instantParticleNums
  //// must already be prepared
  ////


  ////
  ////              0           0
  ////        ┌─────┴─────┬─────┴─────┬───
  ////        │ 0,0       │ 0,1       │ 0,2
  ////        │           │           │
  ////      1 ┤         1 ┤         1 ┤
  ////        │           │           │
  ////        │     0     │     0     │
  ////        ├─────┴─────┼─────┴─────┼─         ->       [(0,0:0), (0,0:1), (0,1:0), (0,1:1), ... ]
  ////        │ 1,0       │ 1,1
  ////        │           │
  ////      1 ┤         1 ┤
  ////        │           │
  ////        │     0     │
  ////        ├─────┴─────┼─
  ////        │ 2,0
  ////
  ////


  ////
  //// C1, C2, C3: center hole
  ////
  if (this->compute_C1 || this->compute_C2 || this->compute_C3) {
    if (this->Cs.size() >= 1) {
      // center major
      const int diff = ((int) this->Cs[0] - (int) holeIs[0]) * (this->Cs.size() == 2 ? this->Ls[1] : 1) * this->halfNumNNs;

      if (this->compute_C1) {
        rotate_copy(this->instC1.begin(),
                    this->instC1.begin() - diff + (diff > 0 ? this->numSites * this->halfNumNNs : 0),
                    this->instC1.end(),
                    this->_instC1.begin());
      }
      if (this->compute_C2) {
        rotate_copy(this->instC2.begin(),
                    this->instC2.begin() - diff + (diff > 0 ? this->numSites * this->halfNumNNs : 0),
                    this->instC2.end(),
                    this->_instC2.begin());
      }
      if (this->compute_C3) {
        rotate_copy(this->instC3.begin(),
                    this->instC3.begin() - diff + (diff > 0 ? this->numSites * this->halfNumNNs : 0),
                    this->instC3.end(),
                    this->_instC3.begin());
      }
    }
    if (this->Cs.size() == 2 && holeIs[1] != this->Cs[1]) {
      // center minor
      const int diff = ((int) this->Cs[1] - (int) holeIs[1]) * this->halfNumNNs;
      for (unsigned i = 0; i < this->Ls[0]; i++) {
        if (this->compute_C1) {
          rotate(this->_instC1.begin() +  i      * this->Ls[1] * this->halfNumNNs,
                 this->_instC1.begin() +  i      * this->Ls[1] * this->halfNumNNs - diff + (diff > 0 ? this->Ls[1] * this->halfNumNNs : 0),
                 this->_instC1.begin() + (i + 1) * this->Ls[1] * this->halfNumNNs);
        }
        if (this->compute_C2) {
          rotate(this->_instC2.begin() +  i      * this->Ls[1] * this->halfNumNNs,
                 this->_instC2.begin() +  i      * this->Ls[1] * this->halfNumNNs - diff + (diff > 0 ? this->Ls[1] * this->halfNumNNs : 0),
                 this->_instC2.begin() + (i + 1) * this->Ls[1] * this->halfNumNNs);
        }
        if (this->compute_C3) {
          rotate(this->_instC3.begin() +  i      * this->Ls[1] * this->halfNumNNs,
                 this->_instC3.begin() +  i      * this->Ls[1] * this->halfNumNNs - diff + (diff > 0 ? this->Ls[1] * this->halfNumNNs : 0),
                 this->_instC3.begin() + (i + 1) * this->Ls[1] * this->halfNumNNs);
        }
      }
    }
  }


  ////
  //// instantParticleNums: center hole
  ////
  if (this->compute_particleCorr) {
    if (this->Cs.size() >= 1) {
      // center major
      const int diff = ((int) this->Cs[0] - (int) holeIs[0]) * (this->Cs.size() == 2 ? this->Ls[1] : 1);

      rotate_copy(this->instantParticleNums.begin(),
                  this->instantParticleNums.begin() + ((diff > 0 ? this->numSites : 0) - diff) * numComps,
                  this->instantParticleNums.end(),
                  this->_instantParticleNums.begin());
    }
    if (this->Cs.size() == 2 && holeIs[1] != this->Cs[1]) {
      // center minor
      const int diff = ((int) this->Cs[1] - (int) holeIs[1]);
      for (unsigned i = 0; i < this->Ls[0]; i++) {
        rotate(this->_instantParticleNums.begin() + numComps * ( i      * this->Ls[1]),
               this->_instantParticleNums.begin() + numComps * ( i      * this->Ls[1] - diff + (diff > 0 ? this->Ls[1] : 0)),
               this->_instantParticleNums.begin() + numComps * ((i + 1) * this->Ls[1]));
      }
    }
  }


  ////
  //// C1, C2, C3: append to outer container
  ////
  this->sampleSpatialCorrelators_store(shift, holeIndexDiff, dt, performConvolution);
}



void Worm::sampleSpatialCorrelators_average () {
  // PROFILING
  // profiler("Worm::sampleSpatialCorrelators_average");



  if (this->Cs.size() != 2) {
    cout << "Worm::sampleSpatialCorrelators: implement/verify for dimensions other than 2!" << endl;
    exit(EXIT_SUCCESS);
  }


  unsigned numConvolutions = 0;
  const long tConv = this->compute_particleConv ?
                     this->pseudoRandom.template Uint<long>(0, tMax - 1) :
                     0;

  const unsigned shift = this->prevFermionicExchangeSign == 1 ?
                         (this->sign == 1 ? 0 : 1) :
                         (this->sign == 1 ? 2 : 3);

  ////
  //// tJ
  ////
  if (modelType == tJ) {

    // load t=0 segments
    for (unsigned i = 0; i < this->numSites; i++) this->segments[i] = this->sites[i][0];

    // load t=0 nearest neighbor correlations
    if (this->compute_C1) {
      vector<unsigned> NNs;
      for (unsigned i = 0; i < this->numSites; i++) {
        this->lattice.getNNs(i, NNs);
        for (unsigned k = 0; k < this->halfNumNNs; k++) {
          const unsigned j = NNs[k];
          this->instC1[i * this->halfNumNNs + k] = (this->segments[i]->pop[0] == this->segments[j]->pop[0]) ? 1 : 0;
          // this->instC1[i * this->halfNumNNs + k] = (this->segments[i]->pop[1] || this->segments[j]->pop[1]) ? -1 : i * this->halfNumNNs + k;
        }
      }
    }

    // load t=0 second nearest neighbor correlations
    if (this->compute_C2) {
      for (unsigned i = 0; i < this->numSites; i++) {
        for (unsigned k = 0; k < this->halfNumNNs; k++) {
          const unsigned j = this->lattice.N2Ns[i][k];
          this->instC2[i * this->halfNumNNs + k] = (this->segments[i]->pop[0] == this->segments[j]->pop[0]) ? 1 : 0;
        }
      }
    }

    // load t=0 third nearest neighbor correlations
    if (this->compute_C3) {
      for (unsigned i = 0; i < this->numSites; i++) {
        for (unsigned k = 0; k < this->halfNumNNs; k++) {
          const unsigned j = this->lattice.N3Ns[i][k];
          this->instC3[i * this->halfNumNNs + k] = (this->segments[i]->pop[0] == this->segments[j]->pop[0]) ? 1 : 0;
        }
      }
    }

    // load t=0 instantParticleNums
    if (this->compute_particleCorr || this->compute_particleConv) {
      for (unsigned i = 0; i < this->numSites; i++) {
        for (unsigned a = 0; a < numComps; a++) {
          this->instantParticleNums[i * numComps + a] = this->segments[i]->pop[a];
        }
      }
    }


    ////
    //// iterate through time
    ////
    long t_beg = 0;
    int nums = -1;
    while (true) {
      nums++;
      long t_end = tMax;
      unsigned i = 0;

      // initiate with negative value so that we can figure out if they have updated
      int holeSiteIndex_1 = -1;
      int holeSiteIndex_2 = -1;


      for (const auto & segment : this->segments) {
        // locate the hole
        if (segment->pop[1] == 1) {
          if (holeSiteIndex_1 == -1) holeSiteIndex_1 = segment->siteIndex;
          else                       holeSiteIndex_2 = segment->siteIndex;
        }

        // find next site change
        if (segment->end && segment->end->t < t_end) {
          t_end = segment->end->t;
          i = segment->siteIndex;
        }
      }


      // hole position
      vector<unsigned> holeIs;
      unsigned holeIndexDiff = 0;   // used only for CX in the case of two carriers
      if (this->Ns[1] == 1) {
        // DEBUG
        if (debug) {
          if (holeSiteIndex_1 == -1) {
            cout << "Worm::sampleSpatialCorrelators_average: ERROR: holeSiteIndex_1 < 0." << endl;
            exit(EXIT_SUCCESS);
          }
          if (holeSiteIndex_2 != -1) {
            cout << "Worm::sampleSpatialCorrelators_average: ERROR: holeSiteIndex_2 != -1." << endl;
            exit(EXIT_SUCCESS);
          }
        }

        holeIs = this->lattice.i2is(holeSiteIndex_1);
      } else if (this->Ns[1] == 2) {
        // DEBUG
        if (debug) {
          if (holeSiteIndex_1 == -1) {
            cout << "Worm::sampleSpatialCorrelators_average: ERROR: holeSiteIndex_1 < 0." << endl;
            exit(EXIT_SUCCESS);
          }
          if (holeSiteIndex_2 == -1) {
            cout << "Worm::sampleSpatialCorrelators_average: ERROR: holeSiteIndex_2 < 0." << endl;
            exit(EXIT_SUCCESS);
          }
        }

        const auto diff = this->spatialIndexDifference(holeSiteIndex_1, holeSiteIndex_2);
        const auto CI = this->Cs[0] * this->latticeSize[1] + this->Cs[1];
        if (diff <= CI) {
          holeIs        = this->lattice.i2is(holeSiteIndex_2);
          holeIndexDiff = diff;
        } else {
          holeIs        = this->lattice.i2is(holeSiteIndex_1);
          holeIndexDiff = this->spatialIndexDifference(holeSiteIndex_2, holeSiteIndex_1);;
        }
      } else {
        cout << "Worm::sampleSpatialCorrelators: ERROR: can only handle 1 or 2 hole(s)." << endl;
        exit(EXIT_SUCCESS);
      }

      // length of interval
      const double dt = (t_end - t_beg) / (double) tMax;

      bool performConvolution = t_beg <= tConv && tConv < t_end;
      if (performConvolution) numConvolutions++;


      ////
      //// compute correlators and convolutions
      ////
      this->sampleSpatialCorrelators_computeAndStore(shift,
                                                     holeIs,
                                                     holeIndexDiff,
                                                     dt,
                                                     performConvolution);

      ////
      //// prepare for next iteration
      ////
      if (t_end != tMax) {
        // find also the connected site
        const auto j = this->segments[i]->end->conn->inco->siteIndex;

        // update segments list
        this->segments[i] = this->segments[i]->end->outg;
        this->segments[j] = this->segments[j]->end->outg;

        // update quantities
        vector<unsigned> NNs;
        for (const unsigned _i : {i, j}) {
          ////
          //// update instC1
          ////
          if (this->compute_C1) {
            this->lattice.getNNs(_i, NNs);
            for (unsigned _k = 0; _k < NNs.size(); _k++) {
              const auto _j = NNs[_k];
              if (_k < this->halfNumNNs) {
                // belongs to _i
                this->instC1[_i * this->halfNumNNs + _k] = (this->segments[_i]->pop[0] == this->segments[_j]->pop[0]) ? 1 : 0;
                // this->instC1[_i * this->halfNumNNs + _k] = (this->segments[_i]->pop[1] || this->segments[_j]->pop[1]) ? -1 : _i * this->halfNumNNs + _k;
              } else {
                // belongs to _jconst
                const auto _l = _k - this->halfNumNNs;   // this is valid for a 1D and 2D square lattice
                this->instC1[_j * this->halfNumNNs + _l] = (this->segments[_i]->pop[0] == this->segments[_j]->pop[0]) ? 1 : 0;
                // this->instC1[_j * this->halfNumNNs + _l] = (this->segments[_i]->pop[1] || this->segments[_j]->pop[1]) ? -1 : _j * this->halfNumNNs + _l;
              }
            }
          }

          ////
          //// update instC2
          ////
          if (this->compute_C2) {
            for (unsigned _k = 0; _k < NNs.size(); _k++) {
              const auto _j = this->lattice.N2Ns[_i][_k];
              if (_k < this->halfNumNNs) {
                // belongs to _i
                this->instC2[_i * this->halfNumNNs + _k] = (this->segments[_i]->pop[0] == this->segments[_j]->pop[0]) ? 1 : 0;
                // this->instC2[_i * this->halfNumNNs + _k] = (this->segments[_i]->pop[1] || this->segments[_j]->pop[1]) ? -1 : _i * this->halfNumNNs + _k;
              } else {
                // belongs to _jconst
                const auto _l = _k - this->halfNumNNs;   // this is valid for a 1D and 2D square lattice
                this->instC2[_j * this->halfNumNNs + _l] = (this->segments[_i]->pop[0] == this->segments[_j]->pop[0]) ? 1 : 0;
                // this->instC2[_j * this->halfNumNNs + _l] = (this->segments[_i]->pop[1] || this->segments[_j]->pop[1]) ? -1 : _j * this->halfNumNNs + _l;
              }
            }
          }

          ////
          //// update instC3
          ////
          if (this->compute_C3) {
            for (unsigned _k = 0; _k < NNs.size(); _k++) {
              const auto _j = this->lattice.N3Ns[_i][_k];
              if (_k < this->halfNumNNs) {
                // belongs to _i
                this->instC3[_i * this->halfNumNNs + _k] = (this->segments[_i]->pop[0] == this->segments[_j]->pop[0]) ? 1 : 0;
                // this->instC3[_i * this->halfNumNNs + _k] = (this->segments[_i]->pop[1] || this->segments[_j]->pop[1]) ? -1 : _i * this->halfNumNNs + _k;
              } else {
                // belongs to _jconst
                const auto _l = _k - this->halfNumNNs;   // this is valid for a 1D and 2D square lattice
                this->instC3[_j * this->halfNumNNs + _l] = (this->segments[_i]->pop[0] == this->segments[_j]->pop[0]) ? 1 : 0;
                // this->instC3[_j * this->halfNumNNs + _l] = (this->segments[_i]->pop[1] || this->segments[_j]->pop[1]) ? -1 : _j * this->halfNumNNs + _l;
              }
            }
          }

          ////
          //// update instantParticleNums
          ////
          if (this->compute_particleCorr) {
            for (unsigned a = 0; a < numComps; a++) {
              this->instantParticleNums[i * numComps + a] = this->segments[i]->pop[a];
              this->instantParticleNums[j * numComps + a] = this->segments[j]->pop[a];
            }
          }
        }

        ////
        //// DEBUG
        ////
        if (debug_major) {
          ////
          //// C1
          ////
          if (this->compute_C1) {
            vector<int> currentNNsCorrelations(this->numSites * this->halfNumNNs);
            for (unsigned i = 0; i < this->numSites; i++) {
              this->lattice.getNNs(i, NNs);
              for (unsigned k = 0; k < this->halfNumNNs; k++) {
                const unsigned j = NNs[k];
                currentNNsCorrelations[i * this->halfNumNNs + k] = (this->segments[i]->pop[0] == this->segments[j]->pop[0]) ? 1 : 0;
                // this->instC1[i * this->halfNumNNs + k] = (this->segments[i]->pop[1] || this->segments[j]->pop[1]) ? -1 : i * this->halfNumNNs + k;
              }
            }

            if (this->instC1 != currentNNsCorrelations) {
              cout << settings::cout::enterRed
                   << "Worm::sampleSpatialCorrelators: this->instC1 != currentNNsCorrelations" << endl
                   << settings::cout::resetStyle;
              if (shutItDown) this->shutDown();
            }
          }


          ////
          //// C2
          ////
          if (this->compute_C2) {
            vector<int> currentC2(this->numSites * this->halfNumNNs);
            for (unsigned i = 0; i < this->numSites; i++) {
              for (unsigned k = 0; k < this->halfNumNNs; k++) {
                const unsigned j = this->lattice.N2Ns[i][k];
                currentC2[i * this->halfNumNNs + k] = (this->segments[i]->pop[0] == this->segments[j]->pop[0]) ? 1 : 0;
              }
            }

            if (this->instC2 != currentC2) {
              cout << settings::cout::enterRed
                   << "Worm::sampleSpatialCorrelators: this->instC2 != currentC2" << endl
                   << settings::cout::resetStyle;
              if (shutItDown) this->shutDown();
            }
          }


          ////
          //// C3
          ////
          if (this->compute_C3) {
            vector<int> currentC3(this->numSites * this->halfNumNNs);
            for (unsigned i = 0; i < this->numSites; i++) {
              for (unsigned k = 0; k < this->halfNumNNs; k++) {
                const unsigned j = this->lattice.N3Ns[i][k];
                currentC3[i * this->halfNumNNs + k] = (this->segments[i]->pop[0] == this->segments[j]->pop[0]) ? 1 : 0;
              }
            }

            if (this->instC3 != currentC3) {
              cout << settings::cout::enterRed
                   << "Worm::sampleSpatialCorrelators: this->instC3 != currentC3" << endl
                   << settings::cout::resetStyle;
              if (shutItDown) this->shutDown();
            }
          }


          ////
          //// instantParticleNums
          ////
          if (this->compute_particleCorr) {
            vector<unsigned> currentParticleNums(numComps * this->numSites);
            for (unsigned i = 0; i < this->numSites; i++) {
              for (unsigned a = 0; a < numComps; a++) {
                currentParticleNums[i * numComps + a] = this->segments[i]->pop[a];
              }
            }

            if (this->instantParticleNums != currentParticleNums) {
              cout << settings::cout::enterRed
                   << "Worm::sampleSpatialCorrelators: this->instantParticleNums != currentParticleNums" << endl
                   << settings::cout::resetStyle;
              if (shutItDown) this->shutDown();
            }
          }

        }

        // update beginning time
        t_beg = t_end;
      } else {
        // were done
        break;
      }

    }

  }



  if (this->compute_particleConv && numConvolutions != 1) {
    cout << "Worm::sampleSpatialCorrelators: numConvolutions != 1  ->  EXIT" << endl;
    exit(EXIT_SUCCESS);
  }
}








void Worm::sampleSpatialCorrelators_slice () {
  // PROFILING
  // profiler("Worm::sampleSpatialCorrelators_slice");

  ////
  //// tJ
  ////
  if (modelType != tJ) {
    cout << "Worm::sampleSpatialCorrelators_slice: implement models other than tJ  ->  EXIT" << endl;
    exit(EXIT_SUCCESS);
  }

  // the fermionic sign shift in data container
  const unsigned shift = this->prevFermionicExchangeSign == 1 ?
                         (this->sign == 1 ? 0 : 1) :
                         (this->sign == 1 ? 2 : 3);

  // select a time slice
  const long tSlice = this->pseudoRandom.template Uint<long>(0, tMax - 1);

  // load segments at chosen time slice
  for (unsigned i = 0; i < this->numSites; i++) {
    unsigned segmentIndex;
    this->findSegmentIndex(i, tSlice, segmentIndex);
    this->segments[i] = this->sites[i][segmentIndex];
  }

  // compute nearest neighbor correlations
  if (this->compute_C1) {
    vector<unsigned> NNs;
    for (unsigned i = 0; i < this->numSites; i++) {
      this->lattice.getNNs(i, NNs);
      for (unsigned k = 0; k < this->halfNumNNs; k++) {
        const unsigned j = NNs[k];
        this->instC1[i * this->halfNumNNs + k] = (this->segments[i]->pop[0] == this->segments[j]->pop[0]) ? 1 : 0;
      }
    }
  }

  // compute second nearest neighbor correlations
  if (this->compute_C2) {
    for (unsigned i = 0; i < this->numSites; i++) {
      for (unsigned k = 0; k < this->halfNumNNs; k++) {
        const unsigned j = this->lattice.N2Ns[i][k];
        this->instC2[i * this->halfNumNNs + k] = (this->segments[i]->pop[0] == this->segments[j]->pop[0]) ? 1 : 0;
      }
    }
  }

  // compute third nearest neighbor correlations
  if (this->compute_C3) {
    for (unsigned i = 0; i < this->numSites; i++) {
      for (unsigned k = 0; k < this->halfNumNNs; k++) {
        const unsigned j = this->lattice.N3Ns[i][k];
        this->instC3[i * this->halfNumNNs + k] = (this->segments[i]->pop[0] == this->segments[j]->pop[0]) ? 1 : 0;
      }
    }
  }

  // compute instantParticleNums
  if (this->compute_particleCorr || this->compute_particleConv) {
    for (unsigned i = 0; i < this->numSites; i++) {
      for (unsigned a = 0; a < numComps; a++) {
        this->instantParticleNums[i * numComps + a] = this->segments[i]->pop[a];
      }
    }
  }


  // initiate with negative value so that we can figure out if they have updated
  int holeSiteIndex_1 = -1;
  int holeSiteIndex_2 = -1;
  for (const auto & segment : this->segments) {
    // locate the hole
    if (segment->pop[1] == 1) {
      if (holeSiteIndex_1 == -1) {
        holeSiteIndex_1 = segment->siteIndex;
      } else {
        holeSiteIndex_2 = segment->siteIndex;
        break;
      }
    }
  }

  // hole position
  vector<unsigned> holeIs;
  unsigned holeIndexDiff = 0;   // used only for CX in the case of two carriers
  if (this->Ns[1] == 1) {
    // DEBUG
    if (debug) {
      if (holeSiteIndex_1 == -1) {
        cout << "Worm::sampleSpatialCorrelators_average: ERROR: holeSiteIndex_1 < 0." << endl;
        exit(EXIT_SUCCESS);
      }
      if (holeSiteIndex_2 != -1) {
        cout << "Worm::sampleSpatialCorrelators_average: ERROR: holeSiteIndex_2 != -1." << endl;
        exit(EXIT_SUCCESS);
      }
    }

    holeIs = this->lattice.i2is(holeSiteIndex_1);
  } else if (this->Ns[1] == 2) {
    // DEBUG
    if (debug) {
      if (holeSiteIndex_1 == -1) {
        cout << "Worm::sampleSpatialCorrelators_average: ERROR: holeSiteIndex_1 < 0." << endl;
        exit(EXIT_SUCCESS);
      }
      if (holeSiteIndex_2 == -1) {
        cout << "Worm::sampleSpatialCorrelators_average: ERROR: holeSiteIndex_2 < 0." << endl;
        exit(EXIT_SUCCESS);
      }
    }

    const auto diff = this->spatialIndexDifference(holeSiteIndex_1, holeSiteIndex_2);
    const auto CI = this->Cs[0] * this->latticeSize[1] + this->Cs[1];
    if (diff <= CI) {
      holeIs        = this->lattice.i2is(holeSiteIndex_2);
      holeIndexDiff = diff;
    } else {
      holeIs        = this->lattice.i2is(holeSiteIndex_1);
      holeIndexDiff = this->spatialIndexDifference(holeSiteIndex_2, holeSiteIndex_1);;
    }
  } else {
    cout << "Worm::sampleSpatialCorrelators: ERROR: can only handle 1 or 2 hole(s)." << endl;
    exit(EXIT_SUCCESS);
  }




  ////
  //// compute correlators and convolutions
  ////
  this->sampleSpatialCorrelators_computeAndStore(shift,
                                                 holeIs,
                                                 holeIndexDiff,
                                                 1,
                                                 true);
}










void Worm::sampleSpatialCorrelators () {
  if (this->compute_average) {
    ////
    //// average over imaginary time
    ////
    this->sampleSpatialCorrelators_average();
  } else {
    ////
    //// choose a time slice at random
    ////
    this->sampleSpatialCorrelators_slice();
  }
}