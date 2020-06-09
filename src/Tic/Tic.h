#ifndef TIC_H
#define TIC_H

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

class Tic {
  private:
    std::chrono::high_resolution_clock::time_point start;

  public:

  private:

  public:
    Tic ();

    double toc () const;

    void reset ();
};

#endif