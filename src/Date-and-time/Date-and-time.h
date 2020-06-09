#ifndef DATE_AND_TIME_H
#define DATE_AND_TIME_H

#include <iostream>
#include <string>
#include <ctime>

struct DateAndTime {

  static std::string getDate ();

  static std::string getTime ();

  static std::string getBoth (const std::string & = "_");

};

#endif