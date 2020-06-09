#include "Hamiltonian.h"
using namespace std;


Hamiltonian::Hamiltonian (
  const Lattice & lattice,
  Model & model
) :
  lattice{lattice},
  model{model}
{
  this->loadParameters();
}

////
//// original Hamiltonian quantities
////
vector<double> Hamiltonian::gett_ija() const {
  return this->t_ija;
}

vector<double> Hamiltonian::getJ_ijab() const {
  return this->J_ijab;
}

vector<double> Hamiltonian::getMu_ia() const {
  return this->mu_ia;
}

array<double, numComps> Hamiltonian::getU_a () const {
  return this->U_a;
}

array<array<double, numComps>, numComps> Hamiltonian::getU_inter_ab () const {
  return this->U_inter_ab;
}

array<double, numComps> Hamiltonian::getEta_a () const {
  return this->eta_a;
}

double Hamiltonian::getEta (unsigned a) const {
  return this->eta_a[a];
}

double Hamiltonian::getDicsoSquredFact (
  const unsigned a,
  const double & beta
) const {
  return pow(this->eta_a[a], 2.) / (beta * this->numSites);
}


////
//// Transformed Hamiltonian quantities
////
vector<double> Hamiltonian::getV_ia() const {
  return this->V_ia;
}


array<array<double, numComps>, numComps> Hamiltonian::getU_ab () const {
  return this->U_ab;
}


array<array<double, numComps>, numComps> Hamiltonian::getU_nn_ab () const {
  return this->U_nn_ab;
}




void Hamiltonian::calculatePotenEnergies_inhom (
  const vector<unsigned long long> & numParticlesAtSite,
  std::vector<double>::iterator      itr,
  const bool                         incr
) const {
  ////
  //// OBS: in original Hamiltonian quantities
  ////
  for (unsigned a = 0; a < numComps; a++) {
    double val = 0;
    for (unsigned i = 0; i < this->numSites; i++) {
      val -= this->mu_ia[this->i_ia(i, a)] * numParticlesAtSite[i * numComps + a];
    }
    *itr = (incr ? *itr : 0) + val / (double) settings::worm::tMax;
    itr++;
  }
}


void Hamiltonian::calculatePotenEnergies_hom (
  const array<unsigned long long, numComps> & numParticles,
  vector<double>::iterator                    itr,
  const bool                                  incr
) const {
  ////
  //// OBS: in original Hamiltonian quantities
  ////
  for (unsigned a = 0; a < numComps; a++) {
    *itr = (incr ? *itr : 0) - this->mu_ia[this->i_ia(0, a)] * numParticles[a] / (double) settings::worm::tMax;
    itr++;
  }
}


void Hamiltonian::calculateInterEnergies (
  const array<unsigned long long, numInters> & numParticlesSquared,
  const array<unsigned long long, numComps> &  numParticles,
  std::vector<double>::iterator                itr,
  const bool                                   incr
) const {
  //
  //        b →
  //          0    1    2    3
  //  a     ┌──────────────────
  //  ↓   0 │ 0    1    2    3
  //      1 │      4    5    6
  //      2 │           7    8
  //      3 │                9
  //

  ////
  //// OBS: in original Hamiltonian quantities
  ////
  for (unsigned a = 0; a < numComps; a++) {
    // intracomponent
    *itr = (incr ? *itr : 0) + 0.5 * this->U_a[a]
                             * (numParticlesSquared[Hamiltonian::i_ab(a, a)] - numParticles[a])
                             / (double) settings::worm::tMax;
    itr++;

    // intercomponent
    for (unsigned b = a + 1; b < numComps; b++) {
      *itr = (incr ? *itr : 0) + (this->U_inter_ab[a][b] * numParticlesSquared[Hamiltonian::i_ab(a, b)])
                               / (double) settings::worm::tMax;
      itr++;
    }
  }
}


unsigned Hamiltonian::i_ab (
  const unsigned a,
  const unsigned b
) {
  return 0.5 * a * (2 * numComps - 1 - a) + b;
}

bool Hamiltonian::isUniform () const {
  return this->uniform;
}

bool Hamiltonian::hasEquivalentComponents () const {
  return this->equivalentComponents;
}


array<double, numComps> Hamiltonian::getChemicalPotential () const {
  if (this->uniform) {
    // copy to array
    array<double, numComps> arr;
    // if (this->equivalentComponents) {
    //   arr.fill(this->mu_ia[0]);
    // } else {
      copy(this->mu_ia.begin(), this->mu_ia.end(), arr.begin());
    // }

    return arr;
  } else {
    cout << "Hamiltonian::getChemicalPotential: implementation needed for uniform = false!" << endl;
    return {};
  }
}

vector<double>::size_type Hamiltonian::i_ia (
  const unsigned i,
  const unsigned a
) const {
  if (this->uniform) {
    // return this->equivalentComponents ? 0 : a;
    return a;
  } else {
    // if (this->equivalentComponents) {
    //   return i;
    // } else {
    return i * numComps + a;
    // }
  }
}

vector<double>::size_type Hamiltonian::i_ija (
  const unsigned i,
  const unsigned j,
  const unsigned a
) const {
  if (this->uniform) {
    // return this->equivalentComponents ? 0 : a;
    return a;
  } else {
    // if (this->equivalentComponents) {
    //   return i * this->numSites + j;
    // } else {
    //   return (i * this->numSites + j) * numComps + a;
    // }
    return (i * this->numSites + j) * numComps + a;
  }
}


vector<double>::size_type Hamiltonian::i_ijab (
  const unsigned i,
  const unsigned j,
  const unsigned a,
  const unsigned b
) const {
  if (this->uniform) {
    // if (this->equivalentComponents) {
    //   return 0;
    // } else {
    return a * numComps + b;
    // }
  } else {
    // if (this->equivalentComponents) {
    //   return i * this->numSites + j;
    // } else {
      return (i * this->numSites + j) * pow(numComps, 2) + a * numComps + b;
    // }
  }
}

vector<double> Hamiltonian::getMu () const {
  return this->mu_ia;
}