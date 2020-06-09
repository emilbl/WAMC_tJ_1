#include "Hamiltonian.h"

using namespace std;

void Hamiltonian::correctChemicalPotential (
  const double & beta,
  const array<unsigned long long, 2 * numComps> & numGsCount
) {
  if (this->uniform) {

    double avgRelDiff = 0;

    for (unsigned a = 0; a < numComps; a++) {

      const double numGsAbove = numGsCount[a * 2 + 1],
                   numGsBelow = numGsCount[a * 2 + 0],
                   relDiff = (numGsAbove - numGsBelow) / max(1., numGsAbove + numGsBelow);

      if ( ! this->equivalentComponents && this->isCanonical[a]) {
        // the components are not equivalent and this one is canonical
        const double dMu = - relDiff / beta;

        this->mu_ia[a] += dMu;
        this->V_ia[a]  += dMu;
      } else if (this->equivalentComponents && this->isCanonical[0]) {
        // prepare the average relative difference in population
        avgRelDiff += relDiff / numComps;
      }
    }

    if (this->equivalentComponents && this->isCanonical[0]) {
      // the components are equivalent
      const double dMu = - avgRelDiff / beta;

      cout << "mu: " << this->mu_ia << " -> ";
      for_each(this->mu_ia.begin(), this->mu_ia.end(), [&](double & d) { d += dMu; });
      for_each(this->V_ia.begin(),  this->V_ia.end(),  [&](double & d) { d += dMu; });
      cout << this->mu_ia << endl;
    }
  } else {
    cout << "Hamiltonian::correctChemicalPotential: implementation needed when uniform = false!" << endl;
    exit(EXIT_SUCCESS);
  }
}