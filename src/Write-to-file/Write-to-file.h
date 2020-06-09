#ifndef WRITE_TO_FILE_H
#define WRITE_TO_FILE_H

// #include <mpi.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <array>

#include "../settings/settings.h"

class WriteToFile {
  public:

  private:
    const std::string outputDir;

  public:
    WriteToFile (const std::string &);

    void mkdir (const std::string &) const;

    // void removeContents (const std::string &) const;

    template<typename T>
    void aVector (const typename std::vector<T>::const_iterator &,
                  const typename std::vector<T>::size_type &,
                  const std::string &,
                  const bool) const;

    void aString (const std::string &,
                  const std::string &) const;

  private:

};

#endif