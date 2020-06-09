#ifndef SGN_H
#define SGN_H

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif