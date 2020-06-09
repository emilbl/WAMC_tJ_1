#include "Write-to-file.h"

using namespace std;


WriteToFile::WriteToFile (
  const string & outputDir
) :
  outputDir{outputDir}
{
  this->mkdir(outputDir);
}


void WriteToFile::mkdir (
  const string & folderName
) const {
  #pragma GCC diagnostic push                        // disables warning temporarily
  #pragma GCC diagnostic ignored "-Wunused-result"   //

  system(("mkdir -p " + folderName).c_str());

  #pragma GCC diagnostic pop                         // enables warning once again
}

// void WriteToFile::removeContents (
//   const string & folderName
// ) {
//   #pragma GCC diagnostic push                        // disables warning temporarily
//   #pragma GCC diagnostic ignored "-Wunused-result"   //

//   system(("rm -r " + pathToRoot + folderName).c_str());

//   #pragma GCC diagnostic pop                         // enables warning once again
// }


template<typename T>
void WriteToFile::aVector (
  const typename vector<T>::const_iterator & beg,
  const typename vector<T>::size_type & size,
  const string & fileName,
  const bool append
) const {
  // open file
  FILE * pFile;
  pFile = fopen((this->outputDir + "/" + fileName).c_str(), append ? "ab" : "wb");

  // save vector
  fwrite(&*beg,
         sizeof(*beg),
         size,
         pFile);

  // close file
  fclose (pFile);
}
template void WriteToFile::aVector<int> (
  const vector<int>::const_iterator &,
  const vector<int>::size_type &,
  const string &,
  const bool
) const;
template void WriteToFile::aVector<unsigned> (
  const vector<unsigned>::const_iterator &,
  const vector<unsigned>::size_type &,
  const string &,
  const bool
) const;
template void WriteToFile::aVector<long> (
  const vector<long>::const_iterator &,
  const vector<long>::size_type &,
  const string &,
  const bool
) const;
template void WriteToFile::aVector<double> (
  const vector<double>::const_iterator &,
  const vector<double>::size_type &,
  const string &,
  const bool
) const;


void WriteToFile::aString (
  const string & theString,
  const string & fileName
) const {
  ofstream myfile;

  myfile.open(this->outputDir + "/" + fileName,
              ios::out | ios::binary);

  if (myfile.is_open()) {

    myfile << theString;
    myfile.close();

  } else {
    cout << settings::cout::enterRed
         << "Write-to-file::binary ERROR: unable to open file  ->  EXIT"
         << settings::cout::resetStyle << endl;
    if (settings::mode::shutItDown) exit (EXIT_FAILURE);
  }
}



// template<typename T>
// void WriteToFile::binary (
//   string pathToFile,
//   vector<T> & data
// ) {
//   ofstream myfile;
//   myfile.open(pathToFile, ios::out | ios::binary);

//   if (myfile.is_open()) {

//     for (unsigned i = 0; i < data.size(); i++) {
//       myfile << data[i] << (i < data.size() - 1 ? "\n" : "");
//     }
//     myfile.close();

//   } else {
//     cout << settings::cout::enterRed << "Write-to-file::binary ERROR: unable to open file  ->  EXIT" << settings::cout::resetStyle << endl;
//     if (settings::mode::shutItDown) exit (EXIT_FAILURE);
//   }

// }
// template void WriteToFile::binary<double> (string, vector<double> &);


// template<typename T>
// void WriteToFile::binaryAppend (
//   string pathToFile,
//   vector<T> & data
// ) {
//   ofstream myfile;
//   myfile.open(pathToFile, ios::out | ios::app | ios::binary);

//   if (myfile.is_open()) {

//     myfile << endl;
//     for (unsigned i = 0; i < data.size(); i++) {
//       myfile << data[i] << (i < data.size() - 1 ? "\n" : "");
//     }
//     myfile.close();

//   } else {
//     cout << settings::cout::enterRed << "Write-to-file::binaryAppend ERROR: unable to open file  ->  EXIT" << settings::cout::resetStyle << endl;
//     if (settings::mode::shutItDown) exit (EXIT_FAILURE);
//   }

// }
// template void WriteToFile::binaryAppend<double> (string, vector<double> &);