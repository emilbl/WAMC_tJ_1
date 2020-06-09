#include "../Worm.h"

using namespace std;


double Worm::getInsertProposalProbability (
  const vector<unsigned> & activeWormComponents
) const {
  if (activeWormComponents.size() == 0) {
    // should insert
    return 1;
  } else if (activeWormComponents.size() == this->maxNumWorms) {
    // must not insert another worm
    return 0;
  } else {
    if (modelType == tJ) {
      // means we are trying to insert a second worm
      if (activeWormComponents[0] == 0) {
        // has a spin worm
        return 0.1;
      } else {
        // has a hole worm
        return 0.1;
      }

    } else {
      cout << "Worm::insertProposalWeight: implement model!" << endl;
      this->shutDown();
    }
  }


  cout << "Worm::insertProposalWeight: this should never occur" << endl;
  exit(EXIT_SUCCESS);
}