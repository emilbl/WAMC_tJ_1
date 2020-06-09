#include "Worm.h"

using namespace std;



static tuple<long, long> transformStartAndEndTime (
  const long & t0,
  const long & t1,
  const long & tMax
) {
  const long t_start = t0                        + (t0 < 0 ? tMax : 0);
  const long t_fin   = t1 + (t0 > t1 ? tMax : 0) + (t0 < 0 ? tMax : 0);

  return make_tuple(t_start, t_fin);
}




double Worm::compute_U_nn () const {
  double exp_U_nn = 0;

  for (unsigned i = 0; i < this->numSites; i++) {
    for (auto seg : this->sites[i]) {

      const long t0 = seg->beg ? seg->beg->t : 0,
                 t1 = seg->end ? seg->end->t : tMax;

      // get nearest neighbors
      vector<unsigned> NNs;
      this->lattice.getNNs(i, NNs);

      for (const auto j : NNs) {

        // no double counting
        if (i > j) continue;

        // fetch the first segment
        unsigned NNsegmentIndex;
        shared_ptr<Segment> NNsegment;
        this->findSegment(j, t0, NNsegment, NNsegmentIndex);

        // loop over NN site segments
        long t_fin = t1 + (t0 > t1 ? tMax : 0);
        long t = t0;
        while (t < t_fin) {

          shared_ptr<Segment> nextNNsegment;
          long t_end;
          if (NNsegment->end) {
            // not crossing t=beta boundary
            nextNNsegment = NNsegment->end->outg;
            t_end         = min(t_fin, nextNNsegment->beg->t);
          } else {
            // crossing t=beta boundary
            if (this->sites[j].size() > 1) {
              // not completely flat site
              nextNNsegment = this->sites[j][1];
              t_end         = min(t_fin, nextNNsegment->beg->t + tMax);
            }  else {
              // complete flat site
              nextNNsegment = NNsegment;
              t_end         = t_fin;
            }
          }

          // update action
          exp_U_nn += this->H.inter_nn(seg->pop, NNsegment->pop)
                    * (t_end - t);

          // prepare for next iteration
          t         = t_end;
          NNsegment = nextNNsegment;
        }
      }
    }
  }

  // normalize
  return exp_U_nn / tMax;
}








double Worm::computeNNinteraction (
  const unsigned                    i,
  const array<unsigned, numComps> & pop,
  const long &                      t0,
  const long &                      t1
) const {
  // get nearest neighbors
  vector<unsigned> NNs;
  this->lattice.getNNs(i, NNs);


  double action = 0;
  for (const auto j : NNs) {
    // to prevent double counting
    if (j >= i) continue;

    // fetch the first segment
    unsigned NNsegmentIndex;
    shared_ptr<Segment> NNsegment;
    this->findSegment(j, t0, NNsegment, NNsegmentIndex);

    // loop over NN site segments
    long t_fin = t1 + (t0 > t1 ? tMax : 0);
    long t = t0;
    while (t < t_fin) {

      shared_ptr<Segment> nextNNsegment;
      long t_end;
      if (NNsegment->end) {
        // not crossing t=beta boundary
        nextNNsegment = NNsegment->end->outg;
        t_end         = min(t_fin, nextNNsegment->beg->t);
      } else {
        // crossing t=beta boundary
        if (this->sites[j].size() > 1) {
          // not completely flat site
          nextNNsegment = this->sites[j][1];
          t_end         = min(t_fin, nextNNsegment->beg->t + tMax);
        }  else {
          // complete flat site
          nextNNsegment = NNsegment;
          t_end         = t_fin;
        }
      }

      // update action
      action += this->H.inter_nn(pop, NNsegment->pop)
              * (t_end - t);

      // prepare for next iteration
      t         = t_end;
      NNsegment = nextNNsegment;
    }
  }

  return action;
}





double Worm::computeNNinteractionDiff (
  const unsigned                    i,
  const array<unsigned, numComps> & pop_aft,
  const array<unsigned, numComps> & pop_bef,
  const long &                      t0,
  const long &                      t1,
  const unsigned                    ignoreSite
) {

  // get nearest neighbors
  vector<unsigned> NNs;
  this->lattice.getNNs(i, NNs);


  // start and end time
  const auto [t_start, t_fin] = transformStartAndEndTime(t0, t1, tMax);


  double action = 0;
  for (const auto j : NNs) {

    // to ignore the site
    if (j == ignoreSite) continue;

    // fetch the first segment
    unsigned NNsegmentIndex;
    shared_ptr<Segment> NNsegment;
    this->findSegment(j, t_start, NNsegment, NNsegmentIndex);

    // loop over NN site segments
    long t = t_start;
    long periodicity = 0;
    while (t < t_fin) {

      shared_ptr<Segment> nextNNsegment;
      long t_end;
      if (NNsegment->end) {
        // not crossing t=beta boundary
        nextNNsegment = NNsegment->end->outg;
        t_end         = min(t_fin, nextNNsegment->beg->t + periodicity);
      } else {
        // crossing t=beta boundary
        periodicity += tMax;

        if (this->sites[j].size() > 1) {
          // not completely flat site
          nextNNsegment = this->sites[j][1];
          t_end         = min(t_fin, nextNNsegment->beg->t + periodicity);
        }  else {
          // complete flat site
          nextNNsegment = NNsegment;
          t_end         = t_fin;
        }
      }

      // update action
      action += this->H.interDiff_nn(pop_aft, pop_bef, NNsegment->pop)
              * (t_end - t);

      // prepare for next iteration
      t         = t_end;
      NNsegment = nextNNsegment;
    }
  }

  return action;
}
double Worm::computeNNinteractionDiff (
  const unsigned                    i,
  const array<unsigned, numComps> & pop_aft,
  const array<unsigned, numComps> & pop_bef,
  const long &                      t0,
  const long &                      t1
) {
  return this->computeNNinteractionDiff(i,
                                        pop_aft,
                                        pop_bef,
                                        t0,
                                        t1,
                                        this->numSites);
}
double Worm::computeNNinteractionDiff (
  const unsigned                    i,
  const array<unsigned, numComps> & pop_bef,
  const unsigned                    actComp,   // difference in population on pop1
  const int                         popDiff,   //
  const long &                      t0,
  const long &                      t1,
  const unsigned                    ignoreSite
) {
  // create population after
  auto pop_aft = pop_bef;
  pop_aft[actComp] += popDiff;

  return this->computeNNinteractionDiff(i, pop_aft, pop_bef, t0, t1, ignoreSite);
}
double Worm::computeNNinteractionDiff (
  const unsigned                    i,
  const array<unsigned, numComps> & pop_bef,
  const unsigned                    actComp,   // difference in population on pop1
  const int                         popDiff,   //
  const long &                      t0,
  const long &                      t1
) {
  return this->computeNNinteractionDiff(i,
                                        pop_bef,
                                        actComp,
                                        popDiff,
                                        t0,
                                        t1,
                                        this->numSites);
}




double Worm::computeNNinteractionDiff (
  const unsigned                    i,
  const array<unsigned, numComps> & pop_aft_i,
  const array<unsigned, numComps> & pop_bef_i,
  const unsigned                    j,
  const array<unsigned, numComps> & pop_aft_j,
  const array<unsigned, numComps> & pop_bef_j,
  const long &                      t0,
  const long &                      t1
) {
  double action = 0;

  action += this->computeNNinteractionDiff(i, pop_aft_i, pop_bef_i, t0, t1, j);
  action += this->computeNNinteractionDiff(j, pop_aft_j, pop_bef_j, t0, t1, i);

  // start and end time
  const auto [t_start, t_fin] = transformStartAndEndTime(t0, t1, tMax);

  // mutual difference
  action += (this->H.inter_nn(pop_aft_i, pop_aft_j) - this->H.inter_nn(pop_bef_i, pop_bef_j)) * (t_fin - t_start);

  return action;
}
double Worm::computeNNinteractionDiff (
  const unsigned                    i,
  const array<unsigned, numComps> & pop_bef_i,
  const unsigned                    actComp_i,
  const int                         popDiff_i,
  const unsigned                    j,
  const array<unsigned, numComps> & pop_bef_j,
  const unsigned                    actComp_j,
  const int                         popDiff_j,
  const long &                      t0,
  const long &                      t1
) {
  // create population after
  auto pop_aft_i = pop_bef_i;
  pop_aft_i[actComp_i] += popDiff_i;
  auto pop_aft_j = pop_bef_j;
  pop_aft_j[actComp_j] += popDiff_j;

  return this->computeNNinteractionDiff(i,
                                        pop_aft_i,
                                        pop_bef_i,
                                        j,
                                        pop_aft_j,
                                        pop_bef_j,
                                        t0,
                                        t1);
}