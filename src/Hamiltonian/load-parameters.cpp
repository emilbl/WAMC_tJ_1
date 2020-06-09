#include "Hamiltonian.h"

using namespace std;


void Hamiltonian::loadParameters () {


  ////
  //// TODO: make sure this is actually working
  ////
  if ( ! this->uniform) {
    cout << settings::cout::enterRed
         << "Hamiltonian::Hamiltonian: TODO: implement \"this->uniform == false\"  ->  EXIT"
         << settings::cout::resetStyle << endl;
    exit(EXIT_FAILURE);
  }


  ////
  //// load t_ija
  ////
  if (has_t) {
    if (this->uniform) {
      for (unsigned a = 0; a < this->numVecElemsPerJump; a++) {
        t_ija[a] = this->model.t(a, {}, {});
      }
    } else {
      for (unsigned j = 0; j < this->numSites; j++) {

        // departing site position
        vector<double> Rj;
        this->lattice.getR(j, Rj);

        // loop over allowed sites to which one may jump to
        vector<unsigned> NNs;
        this->lattice.getNNs(j, NNs);
        for (unsigned i : NNs) {

          // arriving site position
          vector<double> Ri;
          if (this->lattice.boundaryCrossed(i, j)) {
            // boundare crossed, locate ghost cell position
            this->lattice.getGhostSiteR(j, i, Ri);
          } else {
            // boundary not crossed
            this->lattice.getR(i, Ri);
          }

          for (unsigned a = 0; a < this->numVecElemsPerJump; a++) {
            this->t_ija[this->i_ija(i, j, a)] = this->model.t(a, Ri, Rj);
          }
        }
      }
    }
  }


  ////
  //// load J_ijab
  ////
  if (has_J) {
    if (this->uniform) {
      // if (this->equivalentComponents) {
        // J_ijab[0] = this->model.J(0, 0, {}, {});
      // } else {
        for (unsigned a = 0; a < numComps; a++) {
          for (unsigned b = 0; b < numComps; b++) {
            J_ijab[a * numComps + b] = this->model.J(a, b, {}, {});
          }
        }
      // }
    } else {

      if ( ! this->uniform) {
        cout << settings::cout::enterRed
             << "Hamiltonian::Hamiltonian: TODO: implement \"this->uniform = false\" in combination to J  ->  EXIT"
             << settings::cout::resetStyle << endl;
        // MPI_Abort(MPI_COMM_WORLD, 0);
        exit(EXIT_FAILURE);
      }

      // for (unsigned j = 0; j < this->numSites; j++) {

      //   // departing site position
      //   vector<double> Rj;
      //   this->lattice.getR(j, Rj);

      //   // loop over allowed sites to which one may jump to
      //   vector<unsigned> NNs;
      //   this->lattice.getNNs(j, NNs);
      //   for (unsigned i : NNs) {

      //     // arriving site position
      //     vector<double> Ri;
      //     if (this->lattice.boundaryCrossed(i, j)) {
      //       // boundare crossed, locate ghost cell position
      //       this->lattice.getGhostSiteR(j, i, Ri);
      //     } else {
      //       // boundary not crossed
      //       this->lattice.getR(i, Ri);
      //     }

      //     for (unsigned a = 0; a < this->numVecElemsPerJump; a++) {
      //       this->J_ija[this->i_ija(i, j, a)] = this->model.J(a, Ri, Rj);
      //     }
      //   }
      // }
    }
  }


  ////
  //// load mu_ia
  ////
  if (has_mu) {
    if (this->uniform) {
      for (unsigned a = 0; a < this->numVecElemsPerSite; a++) {
        this->mu_ia[a] = this->model.mu({}, a);
      }
    } else {
      for (unsigned i = 0; i < this->numSites; i++) {
        for (unsigned a = 0; a < this->numVecElemsPerSite; a++) {
          vector<double> R;
          this->lattice.getR(i, R);
          this->mu_ia[this->i_ia(i, a)] = this->model.mu(R, a);
        }
      }
    }
  }


  ////
  //// load U_inter and U_a
  //// -> important to fill the whole interaction matrix since depending on
  ////    Worm::allowMultiCompWorm the lower triangular part might be used
  ////
  if (has_U) {
    for (unsigned a = 0; a < numComps; a++) {
      this->U_a[a] = this->model.U(a, a);
      for (unsigned b = 0; b < numComps; b++) {
        this->U_inter_ab[a][b] = a == b ? 0 : this->model.U(a, b);
      }
    }
  }


  ////
  //// load U_nn_ab
  ////
  if (has_U_nn) {
    for (unsigned a = 0; a < numComps; a++) {
      for (unsigned b = 0; b < numComps; b++) {
        this->U_nn_ab[a][b] = this->model.U_nn(a, b);
      }
    }
  }


  ////
  //// load eta_a
  ////
  for (unsigned a = 0; a < numComps; a++) {
    this->eta_a[a] = this->model.eta(a);
  }


  ////
  //// Preform the transformation: μᵃᵢ → μ'ᵃᵢ = μᵃᵢ + ½ Uᵃ
  ////                             Uᵃᵃ = ½ Uᵃ
  //// In order to make the Hamiltonian more compact and the
  //// implementation more efficient. The implemented
  //// Hamiltonian then becomes
  ////
  //// H = ∑  ∑         - tᵢⱼ c⁺ᵃᵢ cᵃⱼ
  ////     ᵃ ⁽ⁱʲ⁾
  ////     ∑  ∑         - μ'ᵃᵢ nᵃᵢ
  ////     ᵃ  ⁱ
  ////     ∑    ∑       + Uᵃᵇ nᵃᵢ nᵇᵢ
  ////    ᵃ>⁼ᵇ  ⁱ
  ////     ∑  ∑         + ηᵃ cᵃᵢ + H.c
  ////     ᵃ  ⁱ
  ////

  ////
  //// load mu'_ia: μ'ᵃᵢ = μᵃᵢ + ½ Uᵃ
  ////
  if (has_U || has_mu) {
    this->V_ia = this->mu_ia;
    if (this->uniform) {
      for (unsigned a = 0; a < this->numVecElemsPerSite; a++) {
        this->V_ia[a] += 0.5 * this->U_a[a];
      }
    } else {
      for (unsigned i = 0; i < this->numSites; i++) {
        for (unsigned a = 0; a < this->numVecElemsPerSite; a++) {
          this->V_ia[this->i_ia(i, a)] += 0.5 * this->U_a[a];
        }
      }
    }
  }

  ////
  //// load U_ab
  ////
  this->U_ab = this->U_inter_ab;
  for (unsigned a = 0; a < numComps; a++) {
    this->U_ab[a][a] = 0.5 * this->U_a[a];
  }



}