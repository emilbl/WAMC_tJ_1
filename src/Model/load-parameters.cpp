#include "Model.h"

using namespace std;

void Model::loadParameters (
  const double &         beta,
  const vector<double> & ts,
  const vector<double> & Js,
  const vector<double> & Us,
  const vector<double> & mus,
  const vector<double> & etas
) {
  this->beta = beta;
  this->etas = etas;

  if (has_t)  this->ts  = ts;
  if (has_J)  this->Js  = Js;
  if (has_U)  this->Us  = Us;
  if (has_mu) this->mus = mus;
}