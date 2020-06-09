#ifndef COUT_VECTOR_H
#define COUT_VECTOR_H

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <map>
#include <utility>

namespace std {

  template<typename T>
  ostream& operator<< (ostream & out, const vector<T> & v) {
    out << "[";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
      out << v[i];
      if (i != last)
        out << ", ";
    }
    out << "]";
    return out;
  }

  template<typename T, size_t size>
  ostream& operator<< (ostream & out, const array<T, size> & v) {
    out << "[";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
      out << v[i];
      if (i != last)
        out << ", ";
    }
    out << "]";
    return out;
  }

  template<typename T1, typename T2>
  ostream& operator<< (ostream & out, const map<T1, T2> & m) {
    out << "[" << endl;
    for (auto it = m.begin(); it != m.end(); it++) {
      out << "\t" << it->first << ": " << it->second;
      if (next(it, 1) != m.end())
        out << ",";
      out << endl;
    }
    out << "]";
    return out;
  }

  template<typename T1, typename T2>
  ostream& operator<< (ostream & out, const pair<T1, T2> & p) {
    out << "(" << p.first << ", " << p.second << ")";
    return out;
  }


}


#endif