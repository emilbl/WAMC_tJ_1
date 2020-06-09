#include "Lattice.h"

using namespace std;

Lattice::Lattice (
  const string &           name,
  const vector<unsigned> & size,
  const vector<double> &   parameters,
  const vector<bool> &     isPeriodic
) :
  name{name},
  isPeriodic{isPeriodic}
{

  if (name == "2d-rectangular") {
    this->rectangular2d(size[0],
                        size[1],
                        parameters[0],
                        parameters[1]);
  }

}

unsigned Lattice::positiveModulo (
  const int i,
  const int n
) {
  return (i % n + n) % n;
}

void Lattice::wrapIs (
  const vector<int> & signedIs,
  vector<unsigned> &  unsignedIs
) const {
  // copy the vector in order to preserv indicies which label lattice points within the Bravais cell
  unsignedIs = vector<unsigned>(signedIs.begin(), signedIs.end());

  // loop only over the indices corresponding to dimension
  for (unsigned d = 0; d < this->numDimensions; d++) {
    unsignedIs[d] = Lattice::positiveModulo(signedIs[d], this->size[d]);
  }
}


vector<unsigned> Lattice::wrapIs (
  const vector<int> & signedIs
) const {
  // copy the vector in order to preserv indicies which label lattice points within the Bravais cell
  vector<unsigned> unsignedIs(signedIs.size());

  // loop only over the indices corresponding to dimension
  for (unsigned d = 0; d < this->numDimensions; d++) {
    unsignedIs[d] = Lattice::positiveModulo(signedIs[d], this->size[d]);
  }

  return unsignedIs;
}



unsigned Lattice::getNumDimensions () const {
  return this->numDimensions;
}

unsigned Lattice::getNumNNs (
  const unsigned i
) const {
  return this->NNs[i].size();
}

void Lattice::getNNs (
  const unsigned i,
  vector<unsigned> & NNs
) const {
  NNs = this->NNs[i];
}

vector<unsigned> Lattice::getCenterSiteIs () const {
  auto center = this->getSize();
  for (auto & val : center) val /= 2;
  return center;
}

void Lattice::getNNsAndDists (
  const unsigned i,
  vector<unsigned> & NNs,
  vector<vector<double> > & dists
) const {
  NNs = this->NNs[i];
  dists = this->dists[i];
}

unsigned Lattice::getNumSites () const {
  return this->numSites;
}


vector<unsigned> Lattice::getSize () const {
  return this->size;
}

vector<int> Lattice::boundaryCrossings (
  const unsigned i,   // jumping to
  const unsigned j    // jumping from
) const {
  return this->BCs[i * this->numSites + j];
}


bool Lattice::boundaryCrossed (
  const unsigned i,   // jumping to
  const unsigned j    // jumping from
) const {
  return this->BC[i * this->numSites + j];
}



void Lattice::getR (
  const unsigned i,
  vector<double> & R
) const {
  // convert single index to spatial indices
  vector<unsigned> is = (this->*this->i2isFuncPtr)(i);

  // allow for negative indices
  vector<int> signedIs(is.begin(), is.end());

  (this->*this->is2RfuncPtr)(signedIs, R);
}

vector<double> Lattice::getR (
  const unsigned i
) const {
  // convert single index to spatial indices
  vector<unsigned> is = (this->*this->i2isFuncPtr)(i);

  // allow for negative indices
  vector<int> signedIs(is.begin(), is.end());

  vector<double> R;
  (this->*this->is2RfuncPtr)(signedIs, R);

  return R;
}


void Lattice::getRdiff (
  const unsigned i,
  const unsigned j,
  vector<double> & Rdiff
) const {

}

void Lattice::getGhostSiteR (
  const unsigned depSiteIndex,   // the ghost site unit cell is touching the unit cell containing this site
  const unsigned arrSiteIndex,   // the ghost site index positive modulo the system site is equal to this one
  vector<double> & R             // the position of the ghost site
) const {

  auto depIs = (this->*this->i2isFuncPtr)(depSiteIndex),
       arrIs = (this->*this->i2isFuncPtr)(arrSiteIndex);

  // the ghost site may be negative
  vector<int> ghostIs(arrIs.begin(), arrIs.end());

  for (unsigned d = 0; d < this->numDimensions; d++) {
    int diff = (int) depIs[d] - (int) arrIs[d];

    if (abs(diff) == (int) (this->size[d] - 1)) {
      // the two sites are on opposite sides of the system in the direction of this dimension

      if (diff > 0) {
        // crossing upper boundary
        ghostIs[d] += this->size[d];
      } else {
        // crossing lower boundary
        ghostIs[d] -= this->size[d];
      }


    } else if (diff != 0) {
      cout << settings::cout::enterRed << "Lattice::getGhostSiteR: ERROR: the two unit cells do not touch one another  ->  EXIT" << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) exit (EXIT_FAILURE);
    }
  }

  (this->*this->is2RfuncPtr)(ghostIs, R);
}


vector<unsigned> Lattice::i2is (unsigned i) const {
  return (this->*this->i2isFuncPtr)(i);
}


bool Lattice::getIsPeriodic (const unsigned d) const {
  return this->isPeriodic[d];
}


vector<bool> Lattice::getIsPeriodic () const {
  return this->isPeriodic;
}


bool Lattice::getIsTranslInvariant () const {
  return find(this->isPeriodic.begin(), this->isPeriodic.end(), true) != this->isPeriodic.end();
}


std::vector<std::array<double, 3> > Lattice::getPvs () const {
  return this->pVs;
}