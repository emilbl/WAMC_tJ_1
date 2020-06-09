#include "../Worm.h"

using namespace std;


bool Worm::sampleDDC () {
  // PROFILING
  // profiler("Worm::sampleDDC");

  if (allowMultiCompWorm) {
    cout << "Worm::sampleDDC: Implement \"allowMultiCompWorm\"" << endl;
    exit(EXIT_SUCCESS);
  }


  ////
  //// both worms must be of the same species
  ////
  if (this->numWorms != 2) return false;
  if (this->actComps[0] != this->actComps[1]) return false;
  const auto actComp = this->actComps[0];

  ////
  //// find the actual head segments heads
  ////
  auto [hs1, hs2, ts1, ts2] = this->findHeadsAndTails(actComp);

  // DEBUG
  if (debug) {
    if ( ! allowAntiWorm) {
      if ( (hs1 != this->headSegments[0] && hs1 != this->headSegments[1]) ||
           (hs2 != this->headSegments[0] && hs2 != this->headSegments[1])    ) {
        cout << "Worm::sampleDDC: ERROR: unexpected worm head" << endl;
        if (shutItDown) this->shutDown();
      }
      if ( (ts1 != this->tailSegments[0] && ts1 != this->tailSegments[1]) ||
           (ts2 != this->tailSegments[0] && ts2 != this->tailSegments[1])    ) {
        cout << "Worm::sampleDDC: ERROR: unexpected worm tail" << endl;
        if (shutItDown) this->shutDown();
      }
    }
    if (hs1 == hs2) {
      cout << "Worm::sampleDDC: ERROR: equal worm heads" << endl;
      if (shutItDown) this->shutDown();
    }
    if (ts1 == ts2) {
      cout << "Worm::sampleDDC: ERROR: equal worm tails" << endl;
      if (shutItDown) this->shutDown();
    }
  }

  ////
  //// must be pairwise located on same sites
  ////
  if ( (hs1->siteIndex != ts1->siteIndex || hs2->siteIndex != ts2->siteIndex) &&
       (hs1->siteIndex != ts2->siteIndex || hs2->siteIndex != ts1->siteIndex)    ) return false;

  // group spatially
  if (hs1->siteIndex == ts2->siteIndex) swap(ts1, ts2);

  // DEBUG
  if (debug) {
    if (hs1->siteIndex != ts1->siteIndex) {
      cout << "Worm::sampleDDC: ERROR: mismatch of site index of first head-tail pair" << endl;
      if (shutItDown) this->shutDown();
    }
    if (hs2->siteIndex != ts2->siteIndex) {
      cout << "Worm::sampleDDC: ERROR: mismatch of site index of second head-tail pair" << endl;
      if (shutItDown) this->shutDown();
    }
  }


  ////
  //// trace one worm to
  //// 1. figure if there is exchange of world lines
  //// 2. compute overlap
  ////
  auto seg = hs1;
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

  // exchange
  const int exchangeSign = ( ! this->isBosonic[actComp] && seg == ts2) ? -1 : 1;


  // worm1
  const long t1_end = hs1->end->t;
  const long t1_beg = seg->beg->t < t1_end ? seg->beg->t : seg->beg->t - tMax;
  const long t1     = t1_end - t1_beg;

  // worm2
  const long t2 = (long) this->numParticles[actComp] - t1;
  const long t2_end = hs2->end->t;
  const long t2_beg = t2_end - t2;


  // smallest possible separation
  const auto t1_mid = t1_beg + 0.5 * t1;
  const auto t2_mid = t2_beg + 0.5 * t2;
  const auto _T = positiveModulo<long, long>(t2_mid - t1_mid, tMax);
  const auto T  = min(_T, tMax - _T);

  // cout << "----------------" << endl;
  // cout << "s = " << exchangeSign << endl;
  // cout << t1_beg / (double) tMax << " -> " << t1_end / (double) tMax << " : " << t1 / (double) tMax << endl;
  // cout << t2_beg / (double) tMax << " -> " << t2_end / (double) tMax << " : " << t2 / (double) tMax << endl;


  // must be overlapping
  if (2*T > t1 + t2) return false;

  int orderSign;
  long overlap;
  if (max(t1, t2) - min(t1, t2) > 2*T) {
    ////
    ////            ∙──────────────────────────∙
    ////                             ∙────∙
    ////
    orderSign = 1;
    overlap   = min(t1, t2);
  } else {
    ////
    ////           ∙───────────────∙
    ////                       ∙───────────────∙
    ////
    orderSign = -1;
    overlap   = 0.5 * (t1 + t2) - T;
  }

  // DEBUG
  if (debug) {
    if (overlap < 0) {
      cout << "Worm::sampleDDC: ERROR: negative overlap" << endl;
      if (shutItDown) this->shutDown();
    }
    if (overlap > min(t1, t2)) {
      cout << "Worm::sampleDDC: ERROR: overlap larger than the shortest worm" << endl;
      if (shutItDown) this->shutDown();
    }
  }


  ////
  //// compute spatial and temporal indices
  ////
  const unsigned spatial_i = this->spatialIndexDifference(hs1->siteIndex, hs2->siteIndex);

  const long maxOverlap = this->_maxNsDiff[actComp] / 2;
  const unsigned overlap_i = overlap * ddc_hist_N_o / maxOverlap;

  const long maxCombinedWorm = this->_maxNsDiff[actComp];
  const unsigned diff_i = (this->numParticles[actComp] - 2*overlap) * ddc_hist_N_d / maxCombinedWorm;

  ////
  //// increase
  ////
  const unsigned i = overlap_i * ddc_hist_N_d * this->numTransl + diff_i * this->numTransl + spatial_i;
  this->ddc_hist[actComp][i]       += this->sign * exchangeSign * orderSign;
  this->ddc_hist_count[actComp][i] += 1;


  return true;
}