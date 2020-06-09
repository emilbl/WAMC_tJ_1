#include "Hamiltonian.h"

using namespace std;


double Hamiltonian::jump (
  const unsigned                    j,
  const array<unsigned, numComps> & pop_j_inco,
  const array<unsigned, numComps> & pop_j_outg,
  const unsigned                    i,
  const array<unsigned, numComps> & pop_i_inco,
  const array<unsigned, numComps> & pop_i_outg
) const {
  ////
  //// for each model determine what type of jump: kinetic, exchange, etc
  ////


  ////
  //// tJ
  ////
  if (modelType == tJ) {
    // what is the population going from j to i
    const int dPop0 = pop_j_inco[0] - pop_j_outg[0],
              dPop1 = pop_j_inco[1] - pop_j_outg[1];


    if (dPop0 != 0 && dPop1 == 0) {
      // spin jump
      return (dPop0 > 0 ? this->t_ija[this->i_ija(i, j, 0)] : this->t_ija[this->i_ija(j, i, 0)])
             * -1;   // the minus sign (-1)^n
    } else if (dPop0 == 0 && dPop1 != 0) {
      // hole jump
      return (dPop1 > 0 ? this->t_ija[this->i_ija(i, j, 1)] : this->t_ija[this->i_ija(j, i, 1)])
             * -1;   // the minus sign (-1)^n
    } else {
      // spin hole exchange
      return (dPop0 > 0 ? this->J_ijab[this->i_ijab(i, j, 0, 1)] : this->J_ijab[this->i_ijab(j, i, 0, 1)])
             * -1;   // the minus sign (-1)^n
    }
  }


  // DEBUG
  cout << settings::cout::enterRed
       << "Hamiltonian::jump: ERROR: unknown model \"" << modelName << "\"  ->  EXIT"<< endl
       << settings::cout::resetStyle;
  exit(EXIT_SUCCESS);
}



void Hamiltonian::jumpDiff (
  const unsigned                    j,
  const array<unsigned, numComps> & pop_j_inco,
  const array<unsigned, numComps> & pop_j_outg,
  const unsigned                    i,
  const array<unsigned, numComps> & pop_i_inco,
  const array<unsigned, numComps> & pop_i_outg,
  const unsigned                    actComp,
  const int                         actPop,
  double &                          base
) const {

  ////
  //// for each model determine what type of jump: kinetic, exchange, etc
  ////


  ////
  //// tJ
  ////
  if (modelType == tJ) {

    // current state
    base /= this->jump(j,
                       pop_j_inco,
                       pop_j_outg,
                       i,
                       pop_i_inco,
                       pop_i_outg);

    // proposed state
    auto _pop_j_inco = pop_j_inco; _pop_j_inco[actComp] += actPop;
    auto _pop_i_outg = pop_i_outg; _pop_i_outg[actComp] += actPop;
    base *= this->jump(j,
                       _pop_j_inco,
                       pop_j_outg,
                       i,
                       pop_i_inco,
                       _pop_i_outg);

    return;
  }


  // DEBUG
  cout << settings::cout::enterRed
       << "Hamiltonian::jumpDiff: ERROR: unknown model \"" << modelName << "\"  ->  EXIT"<< endl
       << settings::cout::resetStyle;
  exit(EXIT_SUCCESS);
}