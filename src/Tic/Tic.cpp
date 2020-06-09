#include "Tic.h"

using namespace std;


Tic::Tic () {
  this->reset();
}


double Tic::toc () const {
  const auto now = chrono::high_resolution_clock::now();

  const auto time_span = chrono::duration_cast<chrono::duration<double> >(now - start);

  return time_span.count();
}


void Tic::reset () {
  this->start = chrono::high_resolution_clock::now();
}