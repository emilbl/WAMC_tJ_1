#ifndef HASH_VECTOR_H
#define HASH_VECTOR_H

#include <functional>
#include <vector>

using namespace std;

// namespace std {

//   template<class T> 
//   struct hash<std::vector<T> > {
//     auto operator() (const std::vector<T> & key) const {
//       std::hash<T> hasher;
//       size_t result = 0;
//       for (size_t i = 0; i < key.size(); ++i) {
//         result = result * 31 + hasher(key[i]); // ??
//       }
//       return result;
//     }
//   };

// }

// // struct vectorHasher {
// //   template<class T>
// //   size_t operator()(const vector<T> & key) const {
// //     hash<T> hasher;
// //     size_t result = 0;
// //     for (size_t i = 0; i < key.size(); ++i) {
// //       result = result * 31 + hasher(key[i]); // ??
// //     }
// //     return result;
// //   }
// // };

#endif