#include "Worm.h"

using namespace std;

void Worm::checkWormConfiguration () const {

  ////
  //// check head and tail
  ////
  if (this->numWorms != this->headSegments.size()) {
    cout << settings::cout::enterRed << "Worm::checkWormConfiguration: this->numWorms != this->headSegments.size()" << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }
  if (this->numWorms != this->tailSegments.size()) {
    cout << settings::cout::enterRed << "Worm::checkWormConfiguration: this->numWorms != this->tailSegments.size()" << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }
  if (this->numWorms != this->headSegmentIndices.size()) {
    cout << settings::cout::enterRed << "Worm::checkWormConfiguration: this->numWorms != this->headSegmentIndices.size()" << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }
  if (this->numWorms != this->tailSegmentIndices.size()) {
    cout << settings::cout::enterRed << "Worm::checkWormConfiguration: this->numWorms != this->tailSegmentIndices.size()" << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) this->shutDown();
  }

  for (unsigned w = 0; w < this->numWorms; w++) {
    // population of the worm
    array<int, numComps> wormPop = {};
    if (allowMultiCompWorm) {
      wormPop = this->wormPop;
    } else if (this->headSegments[w]) {
      wormPop[this->actComps[w]] = this->actPops[w];
    }


    if ( ! allowAntiWorm) {
      // make sure that there are no worms with negative populations
      if (this->actPops[w] <=0 ) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: this->actPops[" << w << "] = " << this->actPops[w] << " <= 0 even though allowAntiWorm == false" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
    }

    // check head
    if (this->sites[this->headSegments[w]->siteIndex].size() <= this->headSegmentIndices[w] || this->headSegments[w] != this->sites[this->headSegments[w]->siteIndex][this->headSegmentIndices[w]]) {
      cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: mismatch between this->headSegmentIndices[" << w << "] and this->headSegments[" << w << "]" << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) this->shutDown();
    }
    if ( ! this->headSegments[w]->end) {
      cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: ! this->headSegments[" << w << "]->end" << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) this->shutDown();
    }
    if (wormPop != Worm::diff<int, unsigned, unsigned>(this->headSegments[w]->pop, this->headSegments[w]->end->outg->pop)) {
      cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: pop imbalance at head" << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) this->shutDown();
    }

    // check tail
    if (this->sites[this->tailSegments[w]->siteIndex].size() <= this->tailSegmentIndices[w] || this->tailSegments[w] != this->sites[this->tailSegments[w]->siteIndex][this->tailSegmentIndices[w]]) {
      cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: mismatch between this->tailSegmentIndices[" << w << "] and this->tailSegments[" << w << "]" << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) this->shutDown();
    }
    if ( ! this->tailSegments[w]->beg) {
      cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: ! this->tailSegments[" << w << "]->beg" << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) this->shutDown();
    }
    if (wormPop != Worm::diff<int, unsigned, unsigned>(this->tailSegments[w]->pop, this->tailSegments[w]->beg->inco->pop)) {
      cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: pop imbalance at tail" << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) this->shutDown();
    }
  }

  ////
  ////
  ////

  // zero population difference
  array<int, numComps> zeroPopDiff = {};

  // the vectors should be consistent with the linked lists
  for (unsigned i = 0; i < this->numSites; i++) {
    unsigned numSegments = this->sites[i].size();

    if ( ! numSegments) {
      cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: numSegments == 0 for lattice site index = " << i << "" << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) this->shutDown();
    }

    if ( ! numSegments%2) {
      cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: an even number of segments for lattice site index = " << i << "" << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) this->shutDown();
    }

    if (this->sites[i][0]->pop != this->sites[i].back()->pop) {
      cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: population imbalance at time boundary for lattice site index = " << i << "" << settings::cout::resetStyle << endl;
      if (settings::mode::shutItDown) this->shutDown();
    }

    for (unsigned j = 0; j < numSegments; j++) {
      auto seg = this->sites[i][j];

      // segment-node connections
      if (j == 0 && seg->beg) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: first segment: segment->beg != nullptr" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
      if (j == numSegments - 1 && seg->end) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: last segment: segment->end != nullptr" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
      if (j > 0 && ! seg->beg) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": segment->beg == nullptr" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
      if (j > 0 && seg->beg->inco != this->sites[i][j - 1]) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": segment->beg->inco incorrect" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
      if (j < numSegments - 1 && ! seg->end->outg) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": segment->end->outg == nullptr" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
      if (j < numSegments - 1 && seg->end->outg != this->sites[i][j + 1]) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": segment->end->outg incorrect" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }

      // node-node connections
      if (seg->beg && seg->beg->conn && seg->beg->conn->conn != seg->beg) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": seg->beg->conn->conn != seg->beg" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
      if (seg->beg && ! seg->beg->conn) {
        const bool isAHead = find(this->headSegments.begin(), this->headSegments.end(), seg->beg->inco) != this->headSegments.end();
        const bool isATail = find(this->tailSegments.begin(), this->tailSegments.end(), seg) != this->tailSegments.end();
        if ( ! isAHead && ! isATail) {
          cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": node without connection" << settings::cout::resetStyle << endl;
          if (settings::mode::shutItDown) this->shutDown();
        }
      }

      // make sure that the distances of the connected nodes are equal but with opposite sign
      if (seg->beg && seg->beg->conn) {
        vector<double> negDist(this->numDims);
        transform(seg->beg->dist.begin(),
                  seg->beg->dist.end(),
                  negDist.begin(),
                  negate<double>());
        if (seg->beg->conn->dist != negDist) {
          cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": the distances do not agree: -seg->beg->dist=" << negDist << " != " << seg->beg->conn->dist << "=seg->beg->conn->dist" << settings::cout::resetStyle << endl;
          if (settings::mode::shutItDown) this->shutDown();
        }
      }

      // make sure that the distance is correct
      if (seg->beg && seg->beg->conn) {
        unsigned I = seg->beg->conn->inco->siteIndex,
                 J = seg->siteIndex;
        vector<double> Ri, Rj, dist(this->numDims);
        this->lattice.getR(J, Rj);
        this->lattice.getR(I, Ri);

        // in case of a boundary crossing
        if (this->lattice.boundaryCrossed(J, I)) {
          this->lattice.getGhostSiteR(J, I, Ri);
        }

        transform(Ri.begin(), Ri.end(), Rj.begin(), dist.begin(), minus<double>());
        if (seg->beg->dist != dist) {
          cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": the distances do not agree: seg->beg->dist=" << seg->beg->dist << " != " << dist << "=dist" << settings::cout::resetStyle << endl;
          if (settings::mode::shutItDown) this->shutDown();
        }
      }

      // make sure that that the times of connected nodes are the same
      if (seg->beg && seg->beg->conn && seg->beg->t != seg->beg->conn->t) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " <<  j << ": seg->beg->t != seg->beg->conn->t" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }

      // times and lengths
      if (seg->beg && (seg->beg->t <= 0 || (long) seg->beg->t >= tMax)) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " <<  j << ": invalid seg->beg->t" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
      if (j > 0 && j < numSegments - 1 && seg->end->t <= seg->beg->t) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": seg->end->t <= seg->beg->t" <<  settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }

      // population
      if (
        seg->end &&
        seg->end->conn &&
        Worm::diff<int, unsigned, unsigned>(seg->pop, seg->end->outg->pop) != Worm::diff<int, unsigned, unsigned>(seg->end->conn->outg->pop, seg->end->conn->inco->pop)
      ) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": population imbalance seg->end" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
      if (seg->beg && zeroPopDiff == Worm::diff<int, unsigned, unsigned>(seg->pop, seg->beg->inco->pop)) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " <<  j - 1 << " -> " << j << ": no population difference" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
      const auto nextSeg = seg->end ? seg->end->outg : this->sites[i][0];
      const auto popDiff = Worm::diff<int, unsigned, unsigned>(seg->pop, nextSeg->pop);
      for (unsigned a = 0; a < numComps; a++) {
        if (abs(popDiff[a]) > 1) {
          cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: abs(popDiff[a]) > 1 at site " << i << ", segment " <<  j << " -> " << j + 1 << settings::cout::resetStyle << endl;
          if (settings::mode::shutItDown) this->shutDown();
        }
      }


      bool suspiciousPop = false;
      for (unsigned k = 0; k < numComps; k++) {
        if (seg->pop[k] > UINT_MAX/2) {
          suspiciousPop = true;
          break;
        }
      }
      if (suspiciousPop) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": suspicious population " << seg->pop << "" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }




      ////
      //// pop and pop difference valid
      ////
      if (hasInvalidPop(seg->pop)) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": invalid population " << seg->pop << "" << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }
      if (j < numSegments - 1 && hasInvalidPopDiff(seg->pop, seg->end->outg->pop)) {
        cout << settings::cout::enterRed << "Worm::checkWormConfiguration: ERROR: site " << i << ", segment " << j << ": invalid population difference " << seg->pop << " -> " << seg->end->outg->pop << settings::cout::resetStyle << endl;
        if (settings::mode::shutItDown) this->shutDown();
      }

    }
  }
}
