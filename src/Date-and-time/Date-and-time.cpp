#include "Date-and-time.h"

using namespace std;


string DateAndTime::getDate () {
  char buffer[80];
  time_t rawtime = time(nullptr);
  strftime(buffer, 80, "%y-%m-%d", localtime(&rawtime));

  return buffer;
}


string DateAndTime::getTime () {
  char buffer[80];
  time_t rawtime = time(nullptr);
  strftime(buffer, 80, "%H:%M:%S", localtime(&rawtime));

  return buffer;
}


string DateAndTime::getBoth (
  const string & delimiter
) {
  return DateAndTime::getDate() + delimiter + DateAndTime::getTime();
}
