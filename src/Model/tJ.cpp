#include "Model.h"

using namespace std;

void Model::tJ (
) {
  this->tFuncPtr   = &Model::tJ_t;
  this->JfuncPtr   = &Model::tJ_J;
  this->muFuncPtr  = &Model::tJ_mu;
  this->UnnFuncPtr = &Model::tJ_U_nn;
  this->etaFuncPtr = &Model::tJ_eta;


  if ( ! (    has_t    == true
           && has_J    == true
           && has_U    == false
           && has_U_nn == true
           && has_mu   == true )
  ) {
    cout << settings::cout::enterRed
         << "Model::tJ: ERROR: incorrect parameters provided  ->  EXIT" << endl
         << settings::cout::resetStyle;
    exit(EXIT_FAILURE);
  }


  ////
  //// [↑↓, hole] -> [boson, fermion], [grand canonical, canonical]
  ////

  // the must be two components
  if (numComps != 2) {
    cout << settings::cout::enterRed
         << "Model::tJ: ERROR: the t-J requires 2 components  ->  EXIT" << endl
         << settings::cout::resetStyle;
    exit(EXIT_FAILURE);
  }

  // the components are not equivalent
  if (this->equivalentComponents) {
    cout << settings::cout::enterRed
         << "Model::tJ: ERROR: the t-J components should not be equivalent  ->  EXIT" << endl
         << settings::cout::resetStyle;
    exit(EXIT_FAILURE);
  }

  // the first component must be bosonic and the second fermionic
  if ( ! this->isBosonic[0] || this->isBosonic[1]) {
    cout << settings::cout::enterRed
         << "Model::tJ: ERROR: incorrect statistics isBosonic=" << this->isBosonic << "!=[1, 0]  ->  EXIT" << endl
         << settings::cout::resetStyle;
    exit(EXIT_FAILURE);
  }

  // the first component must be in the grand canonical ensemble whilst the second must not
  if (this->isCanonical[0] || ! this->isCanonical[1]) {
    cout << settings::cout::enterRed
         << "Model::tJ: ERROR: incorrect ensembles isCanonical=" << this->isCanonical << "!=[0, 1]  ->  EXIT" << endl
         << settings::cout::resetStyle;
    exit(EXIT_FAILURE);
  }

  // no more than two simultaneous worms
  if (this->maxNumWorms > 2) {
    cout << settings::cout::enterRed
         << "Model::tJ: ERROR: this->maxNumWorms=" << this->maxNumWorms << "<= 2  ->  EXIT" << endl
         << settings::cout::resetStyle;
    exit(EXIT_FAILURE);
  }

  // // must allow two simultaneous worms
  // if (this->maxNumWorms != 2) {
  //   cout << settings::cout::enterRed
  //        << "Model::tJ: ERROR: this->maxNumWorms=" << this->maxNumWorms << "!= 2  ->  EXIT" << endl
  //        << settings::cout::resetStyle;
  //   exit(EXIT_FAILURE);
  // }
}

double Model::tJ_t (
  const unsigned a,
  const vector<double> & Ri,
  const vector<double> & Rj
) const {
  if (a == 0) return 0.5 * this->Js[0];
  else        return this->ts[0];
}

double Model::tJ_J (
  const unsigned a,
  const unsigned b,
  const vector<double> & Ri,
  const vector<double> & Rj
) const {
  ////
  //// the direction is determine by the spin
  //// i.e. J_ij if the spin jumps from j -> i
  ////
  if (this->equivalentComponents) {
    return this->ts[0];
  } else {
    // diagonal elements are not allowed since it should
    // not be possible to exchange identical particles
    if (a == b) return 0;

    return this->ts[0];
  }
}

double Model::tJ_mu (
  const vector<double> & R,
  const unsigned a
) {
  if (a == 0) return 0;
  else        return this->mus[a];
}

double Model::tJ_U_nn (
  const unsigned a,
  const unsigned b
) const {
  if (a == b) return 0;
  else        return -0.5 * this->Js[0];
}

double Model::tJ_eta (
  const unsigned a
) const {
  return this->etas[a];
}