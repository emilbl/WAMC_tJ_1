#include "Model.h"

using namespace std;

Model::Model (
  const string &                    modelType,
  const array<bool, numComps>     & isBosonic,
  const array<bool, numComps>     & isCanonical,
  const array<unsigned, numComps> & Ns,
  const vector<double> &            ts,
  const vector<double> &            Js,
  const vector<double> &            Us,
  const vector<double> &            mus,
  const vector<double> &            etas,
  const double &                    beta,
  const unsigned                    numSites,
  const unsigned                    numDims,
  const bool                        uniform,
  const bool                        equivalentComponents,
  const unsigned                    maxNumWorms,
  Pseudorandom &                    pseudoRandom
) :
  isBosonic{isBosonic},
  isCanonical{isCanonical},
  Ns{Ns},
  numSites{numSites},
  numDims{numDims},
  uniform{uniform},
  equivalentComponents{equivalentComponents},
  maxNumWorms{maxNumWorms},
  pseudoRandom{pseudoRandom}
{
  ////
  //// load parameters
  ////
  this->loadParameters(beta,
                       ts,
                       Js,
                       Us,
                       mus,
                       etas);

  ////
  //// set up the corresponding model
  ////
  if (modelType == "tJ") {
    this->tJ();
  } else {
    cout << settings::cout::enterRed
         << "Model::Model: model type \"" << modelType << "\" not found."
         << settings::cout::resetStyle << endl;
    exit(EXIT_SUCCESS);
  }
}


double Model::t (
  unsigned a,
  const vector<double> & Rarr,
  const vector<double> & Rdep
) const {
  return (this->*tFuncPtr)(a, Rarr, Rdep);
}

double Model::J (
  unsigned a,
  unsigned b,
  const vector<double> & Rarr,
  const vector<double> & Rdep
) const {
  return (this->*JfuncPtr)(a, b, Rarr, Rdep);
}


double Model::mu (
  const vector<double> & R,
  const unsigned comp
) {
    return (this->*muFuncPtr)(R, comp);
}


double Model::U (
  const unsigned comp1,
  const unsigned comp2
) const {
  // need to be symmetric!
  // however only the upper triangular part will be used   <-- THE HECK DID I MEAN BY THIS!?
  return (this->*UfuncPtr)(comp1, comp2);
}


double Model::U_nn (
  const unsigned comp1,
  const unsigned comp2
) const {
  return (this->*UnnFuncPtr)(comp1, comp2);
}


double Model::eta (const unsigned comp) const {
  return (this->*etaFuncPtr)(comp);
}



vector<double> Model::diff (
  const vector<double> & V1,
  const vector<double> & V2
) {
  vector<double> V = vector<double>(V1.size());

  for (unsigned i = 0; i < V.size(); i++) {
    V[i] = V1[i] - V2[i];
  }

  return V;
}

vector<double> Model::sum (
  const vector<double> & V1,
  const vector<double> & V2
) {
  vector<double> V = vector<double>(V1.size());

    for (unsigned i = 0; i < V.size(); i++) {
    V[i] = V1[i] + V2[i];
  }

  return V;
}

double Model::dot (
  const vector<double> & V1,
  const vector<double> & V2
) {
  double dot = 0;

  for (unsigned i = 0; i < V1.size(); i++) {
    dot += V1[i] * V2[i];
  }

  return dot;
}

double Model::norm (
  const vector<double> & V
) {
  double norm = 0;

  for (unsigned i = 0; i < V.size(); i++) {
    norm += V[i] * V[i];
  }

  return norm;
}


bool Model::isUniform () const {
  return this->uniform;
}

bool Model::hasEquivalentComponents () const {
  return this->equivalentComponents;
}


std::array<bool, numComps> Model::getIsBosonic () const {
  return this->isBosonic;
}

std::array<bool, numComps> Model::getIsCanonical () const {
  return this->isCanonical;
}

std::array<unsigned, numComps> Model::getNs () const {
  return this->Ns;
}

double Model::getBeta () const {
  return this->beta;
}

unsigned Model::getMaxNumWorms () const {
  return this->maxNumWorms;
}
