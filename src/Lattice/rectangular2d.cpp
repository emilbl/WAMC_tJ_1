#include "Lattice.h"

using namespace std;

void Lattice::rectangular2d (
  unsigned Nx,
  unsigned Ny,
  double   a,
  double   b
) {
  this->i2isFuncPtr = &Lattice::rectangular2dI2is;
  this->is2RfuncPtr = &Lattice::rectangular2dIs2R;

  this->numDimensions = 2;

  // the primitive vectors
  this->pVs = {
    {{a, 0, 0}},
    {{0, b, 0}}
  };

  this->size = {Nx, Ny};

  this->numSites = Nx * Ny;

  this->Rs.resize(this->numSites);
  this->NNs.resize(this->numSites);
  this->dists.resize(this->numSites);

  for (unsigned i = 0; i < this->numSites; i++) {
    this->getR(i, this->Rs[i]);
    this->rectangular2dNNsAndDists(i, this->NNs[i], this->dists[i]);
  }

  // make correct size
  this->BCs.resize(pow(this->numSites, 2));
  this->BC.resize(pow(this->numSites, 2));

  // locate which jumps go across a border
  vector<int> zeroVector = vector<int>(this->numDimensions, 0);
  for (unsigned j = 0; j < this->numSites; j++) {
    for (unsigned i : this->NNs[j]) {
      const auto crossings = this->rectangular2dBoundaryCrossings(i, j);
      this->BCs[i * this->numSites + j] = crossings;
      this->BC[i * this->numSites + j] = (crossings != zeroVector);
    }
  }

  // load N2Ns
  this->N2Ns.resize(this->numSites);
  for (unsigned i = 0; i < this->numSites; i++) {
    this->N2Ns[i] = this->rectangular2dN2N(i);
  }

  // load N3Ns
  this->N3Ns.resize(this->numSites);
  for (unsigned i = 0; i < this->numSites; i++) {
    this->N3Ns[i] = this->rectangular2dN3N(i);
  }
}


unsigned Lattice::rectangular2dIs2i (
  const vector<unsigned> & is
) const {
  // return is[0] + is[1] * (this->size[0]);
  return is[0] * (this->size[1]) + is[1];
}


template<typename T>
vector<T> Lattice::rectangular2dI2is (
  const unsigned i
) const {
  vector<T> is = {0, 0};

  // is[1] = i / this->size[0];
  // is[0] = i - is[1] * this->size[0];
  is[0] = i / this->size[1];
  is[1] = i - is[0] * this->size[1];

  return is;
}


void Lattice::rectangular2dIs2R (
  const vector<int> & is,
  vector<double> & R
) const {
  R = vector<double>(this->numDimensions);

  for (unsigned d = 0; d < this->numDimensions; d++) {
    // which dimension
    for (unsigned j = 0; j < this->numDimensions; j++) {
      // which primitive vector
      R[d] += is[j] * this->pVs[j][d];
    }
  }
}


void Lattice::rectangular2dNNsAndDists (
  const unsigned i,
  vector<unsigned> & NNs,
  vector<vector<double> > & dists
) const {

  // convert total lattice index to spatial lattice indices
  const auto is = this->rectangular2dI2is<int>(i);

  // the neighboring spatial lattice site indices
  // (should be anti clockwise starting from the top)
  vector<vector<int> > NNsIs = { { is[0] - 1,   is[1]     },
                                 { is[0],       is[1] - 1 },
                                 { is[0] + 1,   is[1]     },
                                 { is[0],       is[1] + 1 } };

  ////
  //// filter away periodic boundary conditions
  ////
  unsigned j = 0;
  while (j < NNsIs.size()) {
    if ( ( ! this->isPeriodic[0]) && (NNsIs[j][0] < 0 || NNsIs[j][0] > (int) this->size[0] - 1)) {
      // crossing boundary -> remove
      NNsIs.erase(NNsIs.begin() + j, NNsIs.begin() + j + 1);
      continue;
    }
    if ( ( ! this->isPeriodic[1]) && (NNsIs[j][1] < 0 || NNsIs[j][1] > (int) this->size[1] - 1)) {
      // crossing boundary -> remove
      NNsIs.erase(NNsIs.begin() + j, NNsIs.begin() + j + 1);
      continue;
    }

    // increment
    j++;
  }

  // set appropriate size
  NNs.resize(NNsIs.size());
  dists.resize(NNsIs.size());

  // find the position
  vector<double> Ri;
  this->rectangular2dIs2R(is, Ri);

  // calculate the distance to the neighboring sites
  for (unsigned j = 0; j < NNsIs.size(); j++) {
    // fetch neighboring lattice site position
    vector<double> Rj;
    this->rectangular2dIs2R(NNsIs[j], Rj);

    // calculate the distance
    vector<double> dist(this->numDimensions);
    transform(Rj.begin(), Rj.end(), Ri.begin(), dist.begin(), minus<double>());

    // translate so that the spatial lattice dimension are within the system
    vector<unsigned> wrappedIs;
    this->wrapIs(NNsIs[j], wrappedIs);

    // append to output
    NNs[j] = this->rectangular2dIs2i(wrappedIs);
    dists[j] = dist;
  }
}


vector<int> Lattice::rectangular2dBoundaryCrossings (
  const unsigned i,   // jumping to
  const unsigned j    // jumping from
) const {
  // the spatial index representation
  const auto is = this->rectangular2dI2is<unsigned>(i),
             js = this->rectangular2dI2is<unsigned>(j);

  // the change in each index
  int d0 = (int) js[0] - (int) is[0],
      d1 = (int) js[1] - (int) is[1];

  // the change need to be larger then one in order for there to be a jump across a border
  int w0 = abs(d0) > 1 ? (d0 > 0) - (d0 < 0) : 0,
      w1 = abs(d1) > 1 ? (d1 > 0) - (d1 < 0) : 0;

  return {w0, w1};
}


vector<unsigned> Lattice::rectangular2dN2N (
  const unsigned I
) const {
  // convert total lattice index to spatial lattice indices
  const auto is = this->rectangular2dI2is<int>(I);

  vector<vector<int> > N2NsIs = { { is[0] - 1,   is[1] - 1 },
                                  { is[0] + 1,   is[1] - 1 },
                                  { is[0] + 1,   is[1] + 1 },
                                  { is[0] - 1,   is[1] + 1 } };


  ////
  //// try remove those crossing the boundary
  ////
  unsigned j = 0;
  while (j < N2NsIs.size()) {
    if ( ( ! this->isPeriodic[0]) && (N2NsIs[j][0] < 0 || N2NsIs[j][0] > (int) this->size[0] - 1)) {
      // crossing boundary -> remove
      N2NsIs.erase(N2NsIs.begin() + j, N2NsIs.begin() + j + 1);
      continue;
    }
    if ( ( ! this->isPeriodic[1]) && (N2NsIs[j][1] < 0 || N2NsIs[j][1] > (int) this->size[1] - 1)) {
      // crossing boundary -> remove
      N2NsIs.erase(N2NsIs.begin() + j, N2NsIs.begin() + j + 1);
      continue;
    }

    // increment
    j++;
  }


  ////
  //// convert to global index
  ////
  vector<unsigned> N2Ns;
  for (const auto js : N2NsIs) {
    N2Ns.push_back(this->rectangular2dIs2i(this->wrapIs(js)));
  }

  return N2Ns;
}


vector<unsigned> Lattice::rectangular2dN3N (
  const unsigned I
) const {
  // convert total lattice index to spatial lattice indices
  const auto is = this->rectangular2dI2is<int>(I);

  vector<vector<int> > N3NsIs = { { is[0] - 2,   is[1]     },
                                  { is[0],       is[1] - 2 },
                                  { is[0] + 2,   is[1]     },
                                  { is[0],       is[1] + 2 } };


  ////
  //// try remove those crossing the boundary
  ////
  unsigned j = 0;
  while (j < N3NsIs.size()) {
    if ( ( ! this->isPeriodic[0]) && (N3NsIs[j][0] < 0 || N3NsIs[j][0] > (int) this->size[0] - 1)) {
      // crossing boundary -> remove
      N3NsIs.erase(N3NsIs.begin() + j, N3NsIs.begin() + j + 1);
      continue;
    }
    if ( ( ! this->isPeriodic[1]) && (N3NsIs[j][1] < 0 || N3NsIs[j][1] > (int) this->size[1] - 1)) {
      // crossing boundary -> remove
      N3NsIs.erase(N3NsIs.begin() + j, N3NsIs.begin() + j + 1);
      continue;
    }

    // increment
    j++;
  }


  ////
  //// convert to global index
  ////
  vector<unsigned> N3Ns;
  for (const auto js : N3NsIs) {
    N3Ns.push_back(this->rectangular2dIs2i(this->wrapIs(js)));
  }

  return N3Ns;
}