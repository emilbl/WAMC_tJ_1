#include "Hamiltonian.h"

using namespace std;


void Hamiltonian::setMu (
  const unsigned comp,
  const double & mu
) {
  if (this->uniform) {
    if (this->model.hasEquivalentComponents()) {
      for (unsigned a = 0; a < numComps; a++) {
        // get difference
        const auto dMu = mu - this->mu_ia[a];

        // add difference to both
        this->mu_ia[a] += dMu;
        this->V_ia[a]  += dMu;
      }
    } else {
      // get difference
      const auto dMu = mu - this->mu_ia[comp];

      // add difference to both
      this->mu_ia[comp] += dMu;
      this->V_ia[comp]  += dMu;
    }

  } else {
    cout << "Hamiltonian::setMu: implementation needed when uniform = false!" << endl;
    exit(EXIT_SUCCESS);
  }

}