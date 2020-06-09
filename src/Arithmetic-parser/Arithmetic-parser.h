#ifndef ARITHMETIC_PARSER_H
#define ARITHMETIC_PARSER_H

#include <iostream>
#include <algorithm>
#include <string>

#include "../settings/settings.h"

class ArithmeticParser {
  public:

  private:
    char * expr_string;

  public:
    ArithmeticParser () {};

    template<typename T>
    T parse (std::string _expr) {
      // remove all spaces
      std::string::iterator end_pos = std::remove(_expr.begin(), _expr.end(), ' ');
      _expr.erase(end_pos, _expr.end());


      this->expr_string = &_expr[0];

      return expression<T>();
    }

  private:

    char peek () { return *this->expr_string; }
    char get () { return *this->expr_string++; }


    template<typename T>
    T number () {
      T result = get() - '0';
      while (peek() >= '0' && peek() <= '9') {
        result = 10 * result + get() - '0';
      }
      return result;
    }

    template<typename T>
    T factor () {
      if (peek() >= '0' && peek() <= '9')
        return number<T>();
      else if (peek() == '(') {
        get(); // '('
        T result = expression<T>();
        get(); // ')'
        return result;
      } else if (peek() == '-') {
        get();
        return -factor<T>();
      }

      return (T) 0; // error
    }

    template<typename T>
    T term () {
      T result = factor<T>();
      while (peek() == '*' || peek() == '/')
        if (get() == '*')
          result *= factor<T>();
        else
          result /= factor<T>();

      return result;
    }

    template<typename T>
    T expression () {
      T result = term<T>();
      while (peek() == '+' || peek() == '-')
        if (get() == '+')
          result += term<T>();
        else
          result -= term<T>();
      return result;
    }

};

#endif