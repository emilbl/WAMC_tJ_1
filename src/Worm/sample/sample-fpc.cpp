#include "../Worm.h"

using namespace std;


void Worm::sampleFPC () {
  // PROFILING
  // profiler("Worm::sampleFPC");

  if (allowMultiCompWorm) {
    cout << "Worm::sampleFPC: Implement me!" << endl;
    exit(EXIT_SUCCESS);
  }

  ////
  //// both worms must be of the same species
  ////
  if (this->actComps[0] != this->actComps[1]) return;
  const auto actComp = this->actComps[0];


  ////
  //// find true worm heads
  //// (defined in such a way that the worms always have positive population)
  ////
  auto [seg, seg2, temp1, temp2] = this->findHeadsAndTails(actComp);

  const auto t1_end = seg->end->t;
  const auto t2_end = seg2->end->t;

  ////
  //// trace the worm from true head to tail
  ////
  #pragma message("perhaps we need not do this every time, only after having performed certain updates such as insert/remove + (anti)reconnect?")
  do {
    if ( ! seg->beg) {
      // passing tau=0 boundary
      seg = this->sites[seg->siteIndex].back();
    } else if (
         seg->beg->conn
      && seg->beg->inco->pop[actComp] != seg->pop[actComp]   // this component is jumping
      && seg->beg->conn->inco->pop[actComp]                  // possible to jump away
    ) {
      // always jump (must do the opposite of the previous method)
      seg = seg->beg->conn->inco;
    } else {
      if (seg->beg->inco->pop[actComp] >= seg->pop[actComp]) {
        // stays on same site
        seg = seg->beg->inco;
      } else {
        // found a discontinuity (make the worm as short as possible)
        break;
      }
    }
  } while (true);

  ////
  //// the above convention will interpret
  ////
  ////            ∙──────────────────────────∙
  ////                       ∙────∙
  ////
  //// as (both worms located on the same site)
  ////
  ////           ∙───────────────∙
  ////                       ∙───────────────∙
  ////

  // const auto t1_beg = seg->end->t;
  const auto t1_beg = seg->beg->t;
  const auto t1     = t1_end > t1_beg ? t1_end - t1_beg : tMax + t1_end - t1_beg;

  // deduce
  const auto t2 = (long) this->numParticles[actComp] - t1;
  const auto t2_beg = t2_end > t2 ? t2_end - t2 : tMax + t2_end - t2;

  const auto t1_mid = t1_beg + 0.5 * t1;
  const auto t2_mid = t2_beg + 0.5 * t2;


  ////
  //// DEBUG: trace the other worm
  ////
  if (debug) {
    auto _seg = seg2;

    do {
      if ( ! _seg->beg) {
        // passing tau=0 boundary
        _seg = this->sites[_seg->siteIndex].back();
      } else if (_seg->beg->conn && _seg->beg->inco->pop[actComp] != _seg->pop[actComp]) {
        if ( ! _seg->beg->inco->pop[actComp]) {
          // forced jump
          _seg = _seg->beg->conn->inco;
        } else {
          // stay (must do the opposite of the previous method)
          _seg = _seg->beg->inco;
        }
      } else {
        if ( ! _seg->beg->inco->pop[actComp]) {
          // found (can no longer follow)
          break;
        } else {
          // stay on the same site (make the worm as long as possible)
          _seg = _seg->beg->inco;
        }
      }

    } while (true);

    const auto _t2_beg = _seg->beg->t;
    const auto _t2     = t2_end > _t2_beg ? t2_end - _t2_beg : tMax + t2_end - _t2_beg;

    if (t2 != _t2) {
      cerr << "trySample: ERROR: t2=" << t2 << " vs _t2=" << _t2 << endl;
      cout << currUpdateCount << endl;
      exit(EXIT_SUCCESS);
    }
  }



  ////
  //// This wont work due to it mixing the worldlines, even if they are not overlapping
  ////
  ////    i ∙──────────────────────────∙        ∙───────────────∙
  ////                                     ->
  ////    j            ∙────∙                               ∙───────────────∙
  ////
  //// Hence the overlapping computation wont work
  ////
  // {
  //   shared_ptr<Segment> head1 = nullptr;
  //   shared_ptr<Segment> head2 = nullptr;
  //   shared_ptr<Segment> tail1 = nullptr;
  //   shared_ptr<Segment> tail2 = nullptr;
  //   if (this->actPops[0] > 0) {
  //     head1 = this->headSegments[0];

  //     // true head segment of second worm
  //     // (the tails are not necessarily in order with the heads)
  //     if (this->actPops[1] > 0) {
  //       head2 = this->headSegments[1];
  //       tail1 = this->tailSegments[0];
  //       tail2 = this->tailSegments[1];
  //     } else if (this->tailSegments[0]->beg->inco->pop[actComp] > this->tailSegments[0]->pop[actComp]) {
  //       head2 = this->tailSegments[0]->beg->inco;
  //       tail1 = this->tailSegments[1];
  //       tail2 = this->headSegments[1]->end->outg;
  //     } else {
  //       head2 = this->tailSegments[1]->beg->inco;
  //       tail1 = this->tailSegments[0];
  //       tail2 = this->headSegments[1]->end->outg;
  //     }
  //   } else if (this->actPops[1] > 0) {
  //     head1 = this->headSegments[1];

  //     // true head segment of second worm
  //     // (the tails are not necessarily in order with the heads)
  //     if (this->tailSegments[0]->beg->inco->pop[actComp] > this->tailSegments[0]->pop[actComp]) {
  //       head2 = this->tailSegments[0]->beg->inco;
  //       tail1 = this->tailSegments[1];
  //       tail2 = this->headSegments[0]->end->outg;
  //     } else {
  //       head2 = this->tailSegments[1]->beg->inco;
  //       tail1 = this->tailSegments[0];
  //       tail2 = this->headSegments[0]->end->outg;
  //     }
  //   } else {
  //     // both heads carry negative population
  //     head1 = this->tailSegments[0]->beg->inco;

  //     // true head segment of second worm
  //     head2 = this->tailSegments[1]->beg->inco;

  //     tail1 = this->headSegments[0]->end->outg;
  //     tail2 = this->headSegments[1]->end->outg;
  //   }

  //   long _t1 = head1->end->t > tail1->beg->t ? head1->end->t - tail1->beg->t : head1->end->t + tMax - tail1->beg->t;
  //   long _t2 = head2->end->t > tail2->beg->t ? head2->end->t - tail2->beg->t : head2->end->t + tMax - tail2->beg->t;
  //   if (_t1 + _t2 != (long) this->numParticles[actComp]) {
  //     _t1 = head2->end->t > tail1->beg->t ? head2->end->t - tail1->beg->t : head2->end->t + tMax - tail1->beg->t;
  //     _t2 = head1->end->t > tail2->beg->t ? head1->end->t - tail2->beg->t : head1->end->t + tMax - tail2->beg->t;
  //   }

  //   ////
  //   ////
  //   ////
  //   if (_t1 + _t2 != (long) this->numParticles[actComp]) {
  //     cout << "Worm::trySample: ERROR: _t1 + _t2 != (long) this->numParticles[actComp]" << endl;
  //     cout << _t1 + _t2 << " != " << this->numParticles[actComp] << endl;
  //   }
  //   if (min(_t1, _t2) != min(t1, t2) || max(_t1, _t2) != max(t1, t2)) {
  //     cout << "Worm::trySample: ERROR: min(_t1, _t2) != min(t1, t2) || max(_t1, _t2) != max(t1, t2)" << endl;
  //     cout << min(_t1, _t2) * this->beta / tMax << " vs " << min(t1, t2) * this->beta / tMax << endl;
  //     cout << max(_t1, _t2) * this->beta / tMax << " vs " << max(t1, t2) * this->beta / tMax << endl;
  //     // exit(1);
  //   }
  // }







  // smallest possible separation
  const auto _T = positiveModulo<long, long>(t2_mid - t1_mid, tMax);
  const auto T  = min(_T, tMax - _T);

  const unsigned i0 = (T  / ((double) 0.5 * tMax)) * fpc_hist_N;
  // i1 ≥ i2
  const unsigned i1 = (max(t1, t2) / (double) this->_maxNsDiff[actComp]) * fpc_hist_N;
  const unsigned i2 = (min(t1, t2) / (double) this->_maxNsDiff[actComp]) * fpc_hist_N;
  const unsigned i = fpc_hist_N*(fpc_hist_N*i0 + i1) + i2;

  if (i1 + i2 > fpc_hist_N - 1) {
    cout << "Worm::trySample: i1 + i2 > fpc_hist_N - 1" << endl;
    exit(1);
  }

  ////
  //// compute fermionic sign
  ////
  int fermionicExchangeSign = 1;
  if (this->hasFermionicExchange) {
    cout << "Worm::trySample: implement fermionic exchange sign for fpc" << endl;
    exit(EXIT_SUCCESS);
  }

  ////
  //// time ordering
  ////
  int fermionicTimeOrderingSign = 1;
  if ( ! this->isBosonic[actComp]) {
    // fpc = <T[a a a^† a^†]> = ± <T[a a^† a a^†]>
    // if there is no overlap in time between the worms we must add a minus sign
    if (2*T > t1 + t2) fermionicTimeOrderingSign = -1;
  }

  #pragma message("compute fermionic sign but with a trick: if all fermionic spicies have less (Ns) than two particles -> 1")
  this->fpc_hist[actComp][i]       += this->sign * fermionicExchangeSign * fermionicTimeOrderingSign;
  this->fpc_hist_count[actComp][i] += 1;

}