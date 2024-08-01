#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <functional>
#include <stdexcept>

/* Take out all uses of format to make it compatible with earlier c++ compilers.
#ifdef USE_FMT
#include <fmt/core.h>
using fmt::format;
#else
//#include <format>
//using std::format;
#endif
*/

/**************************************************************************************************/
/* class ParseArgs is a command line argument or switch parser.                                   */
/* Usage                                                                                          */
/* 1) Define the variables that you want the switches to modify. They can be 1 of 4 base types:   */
/*    bool for boolean switches, int64_t for integer switches,                                    */
/*    double for real switches and std::string_view for text switches                             */
/*    They can also be 1 of 4 vector types which are std::vector<base types>                      */
/* 2) Create an instance of PlusArgs (i.e. PlusArgs plusArgs;)                                    */
/* 3) For each desired switch, call the function                                                  */
/*    addArg(const char* arg, const char* desc,  <type>* const valPtr, bool req = false)          */
/*       char* arg is the switch or argument name on the command line.                            */
/*       char* desc is the description that will be printed with the                              */
/*       <type> valPtr is a pointer to a user variable to be modified by the command line arg     */
/*       bool req is an optional flag which indicates this argument or switch is required.        */
/* 3a) Optionally add a lambda that performs some check on the value of the argument. The lambda  */
/*     is of type void and should take as input, a pointer to the variable provided in step 3. If */
/*     the check fails, the lambda should throw any std::exception with text explaining the error.*/
/*     check fails, the lambda should throw any std::exception with text explaining the error.    */
/*     addCheck(const char* argName, std::function<void(<type>*)> check)                          */
/*       char* arg is the switch or argument name on the command line.                            */
/*       std::function<void(ArgValue_t*)> check is a reference to a lambda that does the check.   */
/* 4) Call ParseArgs.parse(const int argc, const char* argv[], int start = 1)                     */
/*       int argc is the number is char* in argv                                                  */
/*       char* argv[] is an array of pointers to char string whish are the args.                  */
/*       int start the index of argv that parsing should start on                                 */
/*       if parse() returns true then no errors we found.                                         */
/*       if parse() returns false, then the program should clean up and exit().                   */
/*                                                                                                */
/* Regardless of how the arguments types defined, the user can enter the arguments in 2 formats:  */
/* 1. the command line switch and value separated by a space                                      */
/* 2. the command line switch and value separated by an '='                                       */
/* However, a boolean switch without a value will set the boolean variable to true.               */
/* Also there is one predefined switch -h which will print the list of user defined switches      */
/* with the provided descriptions.                                                                */
/* The integer and real values can be an equation that will be evaluated, and the result will     */
/* copied copied to the variable.  The equation is evaluated using standard order and precedence  */
/* rules. Operators include + - * / % and ^.  Parentheses can be used to change calculation order.*/
/* For the integer  type, the numbers can be entered as binary, octal, or hex by prefixing the    */
/* number with '0b', '0' or '0x' respectively. Otherwise, the number is interpreted as decimal.   */
/*                                                                                                */
/* Example:                                                                                       */
/*        :                                                                                       */
/*        :                                                                                       */
/* #include "ParseArgs.h"                                                                         */
/*                                                                                                */
/* int main(int argc, char* argv[]) {                                                             */
/*                                                                                                */
/*   //Define variables to be modified by the cmd line args                                       */
/*   bool aval = false;                                                                           */
/*   bool bval = true;                                                                            */
/*   int64_t ival = 0;                                                                            */
/*   double dval = 0.0;                                                                           */
/*   std::string_view sVal = "";                                                                  */
/*   std::vector<bool> vbVals{ false, false };                                                    */
/*                                                                                                */
/*   // Create an instance of ParseArgs and add the desired cmd line args or switches             */
/*   // Optionally add a lambda to perform a check on the parsed value.                           */
/*   ParseArgs parseArgs;                                                                         */
/*   parseArgs.addArg("-a", "Flag only arg", &aval);                                              */
/*   parseArgs.addArg("-b", "Boolean arg", &bval);                                                */
/*   parseArgs.addArg("-i", "Integer arg", &ival);                                                */
/*   parseArgs.addArg("-r", "Real arg", &dval);                                                   */
/*   parseArgs.addCheck([](double* const value) {                                                 */
/*     if (*value > 10.0 || *value < -10.0)                                                       */
/*       throw std::range_error("Value must be in the range of -10 to +10 inclusive");            */
/*     });                                                                                        */
/*   parseArgs.addArg("-s", "String arg", &sVal, true);                                           */
/*   parseArgs.addArg("-vb", "Vector of bools arg", &vbVals);                                     */
/*                                                                                                */
/*   // Call parse with argc and argv.  Exit if there was an error                                */
/*   if (!parseArgs.parse(argc, argv)) exit(0);                                                   */
/*                                                                                                */
/*   // print the current values of all the arguments                                             */
/*   std::cout << parseArgs.getValuesString();                                                    */
/*                                                                                                */
/*   // Since arg parsing is done, the data in parseArgs can optionally be deleted to save space. */
/*   parseArgs.clear();                                                                           */
/*   return 0;                                                                                    */
/* }                                                                                              */
/*                                                                                                */
/* Compilation of ParsArgs.h requires c++20                                                       */
/*                                                                                                */
/* The above example would parse the next 4 lines with without error.                             */
/* -a -b=false -i 5 "-r 6.7 -s abc -vb true,true                                                  */
/* -a -b=false -i=5 "-r=6.7 -s=abc -vb=true,true                                                  */
/* -a -b=false -i=3+2 "-r=(10+3.4)/2 -s=abc -vb=true,true                                         */
/* -a -b=false -i=0b101 "-r=(10+3.4)/2 -s=abc -vb=true,true                                       */
/**************************************************************************************************/


class ParseEq {
public:
  static double parseRealEquation(const std::string_view& in, size_t& strIdx, size_t precedence, int64_t& parenCnt) {
    double currentValue = 0;
    std::string digits;

    // skip over leading spaces
    while (strIdx < in.length() && in[strIdx] == ' ') {
      strIdx++;
    }
    // then check for open paren and push 
    if (strIdx < in.length() && in[strIdx] == '(') {
      strIdx++;
      parenCnt++;
      currentValue = handleParen(in, strIdx, 0, parenCnt);
    }
    else { // else must be a number
      digits = in.substr(strIdx);
      size_t numDigits = 0;
      currentValue = std::stod(digits, &numDigits);
      strIdx += numDigits;
    }

    // next process an operator.
    while (strIdx < in.length()) {

      if (in[strIdx] == '+' || in[strIdx] == '-') {
        if (precedence >= 1) return currentValue;
        char op = in[strIdx++];
        double returnValue = parseRealEquation(in, strIdx, 1, parenCnt);
        if (op == '+') currentValue += returnValue;
        else currentValue -= returnValue;
        continue;
      }

      if (in[strIdx] == '*' || in[strIdx] == '/') {
        if (precedence >= 2) return currentValue;
        char op = in[strIdx++];
        double returnValue = parseRealEquation(in, strIdx, 2, parenCnt);
        if (op == '*') currentValue *= returnValue;
        else if (op == '/') currentValue /= returnValue;
        continue;
      }

      if (in[strIdx] == '^') {
        if (precedence >= 3) return currentValue;
        char op = in[strIdx++];
        double returnValue = parseRealEquation(in, strIdx, 2, parenCnt);
        if (op == '^') currentValue = std::pow(currentValue, returnValue);
        continue;
      }

      if (in[strIdx] == ')') {
        if (parenCnt <= 0)
          throw std::invalid_argument("mismatched parentheses");
        return currentValue;
      }

      if (in[strIdx] == ' ') {
        strIdx++; // skip over spaces
        continue;
      }

      throw std::invalid_argument("unexpected character");
      strIdx++;
    }
    return(currentValue);
  }

  static double handleParen(const std::string_view& in, size_t& strIdx, size_t precedence, int64_t& parenCnt) {
    double currentValue = parseRealEquation(in, strIdx, 0, parenCnt);
    if (strIdx < in.length() && in[strIdx] == ')') parenCnt--;
    strIdx++;
    return currentValue;
  }

  static int64_t parseIntEquation(const std::string_view& in, size_t& strIdx, size_t precedence, int64_t& parenCnt) {
    int64_t currentValue = 0;
    std::string digits;

    // skip over leading spaces
    while (strIdx < in.length() && in[strIdx] == ' ') {
      strIdx++;
    }
    // then check for open paren and push 
    if (strIdx < in.length() && in[strIdx] == '(') { 
      strIdx++;
      parenCnt++;
      currentValue = handleIntParen(in, strIdx, 0, parenCnt);
    }
    else { // else must be a number
      int radix= 0;
      if (strIdx < in.length() && in[strIdx] == '0') { // is radix other than 10?
        strIdx++;
        if (strIdx < in.length() && (in[strIdx] == 'b' || in[strIdx] == 'B')) { //is radix binary?
          radix = 2;
          strIdx++;
        }
        else if (strIdx < in.length() && (in[strIdx] == 'x' || in[strIdx] == 'X')) { //is radix hex? 
          radix = 16;
          strIdx++;
        }
        else { // then must be octal
          radix = 8;
        }
      }
      else { // then radix must be 10
        radix = 10;
      }
      digits = in.substr(strIdx);
      size_t numDigits = 0;
      currentValue = std::stoll(digits, &numDigits, radix);
      strIdx += numDigits;
    }

    // next process an operator.
    while (strIdx < in.length()) {

      if (in[strIdx] == '+' || in[strIdx] == '-') {
        if (precedence >= 1) return currentValue;
        char op = in[strIdx++];
        int64_t returnValue = parseIntEquation(in, strIdx, 1, parenCnt);
        if (op == '+') currentValue += returnValue;
        else currentValue -= returnValue;
        continue;
      }

      if (in[strIdx] == '*' || in[strIdx] == '/' || in[strIdx] == '%') {
        if (precedence >= 2) return currentValue;
        char op = in[strIdx++];
        int64_t returnValue = parseIntEquation(in, strIdx, 2, parenCnt);
        if (op == '*') currentValue *= returnValue;
        else if (op == '/') currentValue /= returnValue;
        else if (op == '%') currentValue %= returnValue;
        continue;
      }

      if (in[strIdx] == '^') {
        if (precedence >= 3) return currentValue;
        char op = in[strIdx++];
        int64_t returnValue = parseIntEquation(in, strIdx, 2, parenCnt);
        if (op == '^') currentValue = (int64_t)std::pow(currentValue, returnValue);
        continue;
      }

      if (in[strIdx] == ')') {
        if (parenCnt <= 0)
          throw std::invalid_argument("mismatched parentheses");
        return currentValue;
      }

      if (in[strIdx] == ' ') {
        strIdx++; // skip over spaces
        continue;
      }

      throw std::invalid_argument("unexpected character");
      strIdx++;
    }
    return(currentValue);
  }

  static int64_t handleIntParen(const std::string_view& in, size_t& strIdx, size_t precedence, int64_t& parenCnt) {
    int64_t currentValue = parseIntEquation(in, strIdx, 0, parenCnt);
    if (strIdx < in.length() && in[strIdx] == ')') parenCnt--;
    strIdx++;
    return currentValue;
  }

};


class ParseArgs {

  /*** ArgBase is base struct for structs that hold information about the specific arguments ***/
  struct ArgBase {
    const char* argName;      // name of the argument or switch
    const char* description;  // description of the argument or switch
    const bool required;      // boolean indicating whether this arg or switch is required
    bool found;               // Boolean indicating if this switch was found during parsing

    ArgBase(const char* arg, const char* desc, bool req) : argName{ arg }, description{ desc }, required{ req }, found{ false } {}

    virtual void  parse(std::string_view s) = 0;
    virtual std::string toString() = 0;
    virtual std::string type() = 0;
    virtual ~ArgBase() {};
  };

  /**** ArgBool holds information about specific single boolean arguments ****/
  struct ArgBool : ArgBase {
    bool* const value;  // pointer to the value
    std::function<void(bool*)> check = nullptr; // pointer to a potential value check

    ArgBool(const char* arg, const char* desc, bool* valPtr, bool req) : ArgBase(arg, desc, req), value{ valPtr } {} // constructor

    void parse(std::string_view s) override {  // parse for the specific type
      if (s == "")
        *value = true;
      else if (s == "true")
        *value = true;
      else if (s == "false")
        *value = false;
      else throw std::invalid_argument("Must be 'true' or 'false'");
      if (check != nullptr) check(value);
    }

    std::string toString() override { return *value ? "true" : "false"; }  // conversion to a string

    std::string type() override { return "boolean"; }  // string indication the type.
  };

  /**** ArgInt holds information about specific integer arguments ****/
  struct ArgInt : ArgBase {
    int64_t* const value;
    std::function<void(int64_t*)> check = nullptr;

    ArgInt(const char* arg, const char* desc, int64_t* const valPtr, bool req) : ArgBase(arg, desc, req), value{ valPtr } {
    }

    void parse(std::string_view s) override {
      if (s == "")
        throw std::invalid_argument("Must include an integer value");
      else {
        size_t strIdx = 0;
        int64_t parenCount = 0;
        *value = ParseEq::parseIntEquation(std::string(s), strIdx, 0, parenCount);
        if (parenCount != 0)
          throw std::invalid_argument("mismatched parentheses");
        else if (strIdx != s.length())
          throw std::invalid_argument("parse error");
        if (check != nullptr) check(value);
      }
    }

    std::string toString() override { return std::to_string(*value); }

    std::string type()  override { return "integer number"; }
  };

  /**** ArgReal holds information about specific real number arguments ****/
  struct ArgReal : ArgBase {
    double* const value;
    std::function<void(double*)> check = nullptr;

    ArgReal(const char* arg, const char* desc, double* const valPtr, bool req) : ArgBase(arg, desc, req), value{ valPtr } {}

    void parse(std::string_view s) override {
      if (s == "")
        throw std::invalid_argument(std::string("Must include a real number value"));
      else {
        size_t strIdx = 0;
        int64_t parenCount = 0;
        *value = ParseEq::parseRealEquation(s, strIdx, 0, parenCount);
        if (parenCount != 0)
          throw std::invalid_argument("mismatched parentheses");
        else if (strIdx != s.length())
          throw std::invalid_argument("parse error");
        if (check != nullptr) check(value);
      }
    }

    std::string toString() override { return std::to_string(*value); }

    std::string type() override { return "real number"; }
  };

  /**** ArgText holds information about specific text arguments ****/
  struct ArgText : ArgBase {
    std::string_view* const value;
    std::function<void(std::string_view*)> check = nullptr;

    ArgText(const char* arg, const char* desc, std::string_view* const valPtr, bool req) : ArgBase(arg, desc, req), value{ valPtr } {}

    void parse(std::string_view s) override {
      if (s == "")
        throw std::invalid_argument("Must include a text value");
      else
        *value = s;

      if (check != nullptr) check(value);
    }

    std::string toString() override { return '"' + std::string(*value) + '"'; }

    std::string type() override { return "text"; }
  };

public:

  /** These addArg member functions build and modify an ArgBool object **/
  void addArg(const char* arg, const char* desc, bool* const valPtr, bool req = false) {
    addArgT<ArgBool, bool>(arg, desc, valPtr, req);
  }
  void addCheck(std::function<void(bool*)> check) {
    addCheckT<ArgBool, bool>(check);
  }
  void addCheck(const char* argName, std::function<void(bool*)> check) {
    addCheckT<ArgBool, bool>(argName, check);
  }

  /** These addArg member functions build and modify an ArgInt objects **/
  void addArg(const char* arg, const char* desc, int64_t* const valPtr, bool req = false) {
    addArgT<ArgInt, int64_t>(arg, desc, valPtr, req);
  }
  void addCheck(std::function<void(int64_t*)> check) {
    addCheckT<ArgInt, int64_t>(check);
  }
  void addCheck(const char* argName, std::function<void(int64_t*)> check) {
    addCheckT<ArgInt, int64_t>(argName, check);
  }

  /** These addArg member functions build and modify an ArgReal objects **/
  void addArg(const char* arg, const char* desc, double* const valPtr, bool req = false) {
    addArgT<ArgReal, double>(arg, desc, valPtr, req);
  }
  void addCheck(std::function<void(double*)> check) {
    addCheckT<ArgReal, double>(check);
  }
  void addCheck(const char* argName, std::function<void(double*)> check) {
    addCheckT<ArgReal, double>(argName, check);
  }

  /** These addArg member functions build and modify an ArgText objects **/
  void addArg(const char* arg, const char* desc, std::string_view* const valPtr, bool req = false) {
    addArgT<ArgText, std::string_view>(arg, desc, valPtr, req);
  }
  void addCheck(std::function<void(std::string_view*)> check) {
    addCheckT<ArgText, std::string_view>(check);
  }
  void addCheck(const char* argName, std::function<void(std::string_view*)> check) {
    addCheckT<ArgText, std::string_view>(argName, check);
  }

#ifndef REMOVE_VECTOR_ARGUMENT_CODE
private:

  /**** ArgBoolArray holds information about specific boolean array arguments ****/
  struct ArgBoolArray : ArgBase {
    std::vector<bool>* const values;
    std::function<void(std::vector<bool>*)> check = nullptr; // pointer to a potential value check

    ArgBoolArray(const char* arg, const char* desc, std::vector<bool>* valPtr, bool req) : ArgBase(arg, desc, req), values{ valPtr } {}

    void parse(std::string_view s) override {
      if (s == "")
        throw std::invalid_argument("Must include a list of boolean values");
      else {
        auto pVals = split(s, ',');
        size_t len = values->size();
        if (pVals.size() != len)
          throw std::invalid_argument("Number of values doesn't match array size of " + std::to_string(len));
        for (int idx = 0; idx < len; idx++) {
          if (pVals[idx] == "true")
            values->at(idx) = true;
          else if (pVals[idx] == "false")
            values->at(idx) = false;
          else throw std::invalid_argument("Boolean values must be 'true' or 'false'");
        }
      }
      if (check != nullptr) check(values);
    }

    std::string toString() override {
      std::string tmp;
      for (int i = 0; i < values->size(); i++) {
        if (i != 0)  tmp += ',';
        tmp += values->at(i) ? "true" : "false";
      }
      return tmp;
    }

    std::string type() override { return "boolean array"; }
  };

  /**** ArgIntArray holds information about specific integer array arguments ****/
  struct ArgIntArray : ArgBase {
    std::vector<int64_t>* const values;
    std::function<void(std::vector<int64_t>*)> check = nullptr; // pointer to a potential value check

    ArgIntArray(const char* arg, const char* desc, std::vector<int64_t>* valPtr, bool req) : ArgBase(arg, desc, req), values{ valPtr } {}

    void parse(std::string_view s) override {
      if (s == "")
        throw std::invalid_argument("Must include a list of integer values");
      else {
        auto pVals = split(s, ',');
        size_t len = values->size();
        if (pVals.size() != len)
          throw std::invalid_argument("Number of values doesn't match array size of " + std::to_string(len));
        for (int idx = 0; idx < len; idx++) {
          size_t last_c = 0;
          values->at(idx) = std::stoll(std::string(pVals[idx]), &last_c);
          if (last_c != pVals[idx].size())
            throw std::invalid_argument("Must include integer number values");
        }
      }
      if (check != nullptr) check(values);
    }

    std::string toString() override {
      std::string tmp;
      for (int i = 0; i < values->size(); i++) {
        if (i != 0)  tmp += ',';
        tmp += std::to_string(values->at(i));
      }
      return tmp;
    }

    std::string type() override { return "integer array"; }
  };

  /**** ArgRealArray holds information about specific integer array arguments ****/
  struct ArgRealArray : ArgBase {
    std::vector<double>* const values;
    std::function<void(std::vector<double>*)> check = nullptr; // pointer to a potential value check

    ArgRealArray(const char* arg, const char* desc, std::vector<double>* valPtr, bool req) : ArgBase(arg, desc, req), values{ valPtr } {}

    void parse(std::string_view s) override {
      if (s == "")
        throw std::invalid_argument("Must include a list of real values");
      else {
        auto pVals = split(s, ',');
        size_t len = values->size();
        if (pVals.size() != len)
          throw std::invalid_argument("Number of values doesn't match array size of " + std::to_string(len));
        for (int idx = 0; idx < len; idx++) {
          size_t last_c = 0;
          values->at(idx) = std::stod(std::string(pVals[idx]), &last_c);
          if (last_c != pVals[idx].size())
            throw std::invalid_argument("Must include a real number values");
        }
      }
      if (check != nullptr) check(values);
    }

    std::string toString() override {
      std::string tmp;
      for (int i = 0; i < values->size(); i++) {
        if (i != 0)  tmp += ',';
        tmp += std::to_string(values->at(i));
      }
      return tmp;
    }

    std::string type() override { return "real array"; }
  };

  /**** ArgTextArray holds information about specific text array arguments ****/
  struct ArgTextArray : ArgBase {
    std::vector<std::string_view>* const values;
    std::vector<std::string> localValues;
    std::function<void(std::vector<std::string_view>*)> check = nullptr; // pointer to a potential value check

    ArgTextArray(const char* arg, const char* desc, std::vector<std::string_view>* valPtr, bool req) : ArgBase(arg, desc, req), values{ valPtr } {}

    void parse(std::string_view s) override {
      if (s == "")
        throw std::invalid_argument("Must include a list of text values");
      else {
        auto pVals = split(s, ',');
        size_t len = values->size();
        if (pVals.size() != len)
          throw std::invalid_argument("Number of values doesn't match array size of " + std::to_string(len));
        for (int idx = 0; idx < len; idx++) {
          values->at(idx) = pVals[idx];
        }
      }
      if (check != nullptr) check(values);
    }

    std::string toString() override {
      std::string tmp;
      for (int i = 0; i < values->size(); i++) {
        if (i != 0)  tmp += ',';
        tmp += '"' + std::string(values->at(i)) + '"';
      }
      return tmp;
    }

    std::string type() override { return "text array"; }
  };

public:
  /** These addArg member functions build and modify an ArgBoolArray objects **/
  void addArg(const char* arg, const char* desc, std::vector<bool>* const valPtr, bool req = false) {
    addArgT<ArgBoolArray, std::vector<bool>>(arg, desc, valPtr, req);
  }
  void addCheck(std::function<void(std::vector<bool>*)> check) {
    addCheckT<ArgBoolArray, std::vector<bool>>(check);
  }
  void addCheck(const char* argName, std::function<void(std::vector<bool>*)> check) {
    addCheckT<ArgBoolArray, std::vector<bool>>(argName, check);
  }

  /** These addArg member functions build and modify an ArgIntArray objects **/
  void addArg(const char* arg, const char* desc, std::vector<int64_t>* const valPtr, bool req = false) {
    addArgT<ArgIntArray, std::vector<int64_t>>(arg, desc, valPtr, req);
  }
  void addCheck(std::function<void(std::vector<int64_t>*)> check) {
    addCheckT<ArgIntArray, std::vector<int64_t>>(check);
  }
  void addCheck(const char* argName, std::function<void(std::vector<int64_t>*)> check) {
    addCheckT<ArgIntArray, std::vector<int64_t>>(argName, check);
  }

  /** These addArg member functions build and modify an ArgRealArray objects **/
  void addArg(const char* arg, const char* desc, std::vector<double>* const valPtr, bool req = false) {
    addArgT<ArgRealArray, std::vector<double>>(arg, desc, valPtr, req);
  }
  void addCheck(std::function<void(std::vector<double>*)> check) {
    addCheckT<ArgRealArray, std::vector<double>>(check);
  }
  void addCheck(const char* argName, std::function<void(std::vector<double>*)> check) {
    addCheckT<ArgRealArray, std::vector<double>>(argName, check);
  }

  /** These addArg member functions build and modify an ArgBoolArray objects **/
  void addArg(const char* arg, const char* desc, std::vector<std::string_view>* const valPtr, bool req = false) {
    addArgT<ArgTextArray, std::vector<std::string_view>>(arg, desc, valPtr, req);
  }
  void addCheck(std::function<void(std::vector<std::string_view>*)> check) {
    addCheckT<ArgTextArray, std::vector<std::string_view>>(check);
  }
  void addCheck(const char* argName, std::function<void(std::vector<std::string_view>*)> check) {
    addCheckT<ArgTextArray, std::vector<std::string_view>>(argName, check);
  }

#endif // REMOVE_VECTOR_ARGUMENT_CODE

private:

  /** split is a static utility member function that splits a string with a number of separators into separate strings.
  *** string in is the input string
  *** char sep is the separator character about which the string will be separated
  *** returned is a vector of individual separated strings **/

  static std::vector<std::string_view> split(std::string_view in, char sep) {
    std::vector<std::string_view> out;
    std::string_view word;
    size_t start_idx = 0;
    while (true) {
      size_t end_idx = in.find(sep, start_idx);
      word = in.substr(start_idx, end_idx - start_idx);
      out.push_back(word);
      if (end_idx == std::string::npos) break;
      start_idx = end_idx + 1;
    }
    return out;
  }


  /** The following template addArg member function builds a specific ArgValue_t object, each for a different arg type
  *** Specific versions of this function for each allowed type are built in the public section of this class
  *** Parameters:
  *** char* arg is the switch or argument name on the command line.
  *** char* desc is the description that will be printed with the
  *** ArgValue_t valPtr is a pointer to a user variable to be modified by the command line arg
  *** bool req = false is a flag which indicates these argment or switch is required.
  **/
  template<class ArgBase_t, class ArgValue_t>
  void addArgT(const char* arg, const char* desc, ArgValue_t* const valPtr, bool req) {
    auto* tempArgBase = new ArgBase_t(arg, desc, valPtr, req);
    std::pair< std::string, ArgBase* > argEntry{ arg, (ArgBase*)tempArgBase };
    argsMap.insert(argEntry);
    lastMapEntry = argEntry;
  }

  /** The following template addCheck member function adds a "check" lambda of function to the last
  *** ArgBase_t object added to the argMap..
  *** Specific versions of this function for each allowed type are built in the public section of this class
  *** Parameters:
  *** std::function<void(ArgValue_t*)> check is a reference to a lambda that does the check on the args value.
  ***/
  template<class ArgBase_t, class ArgValue_t>
  void addCheckT(std::function<void(ArgValue_t*)> check) {
    if (typeid(ArgBase_t) == typeid(*lastMapEntry.second)) {
      static_cast<ArgBase_t*>(lastMapEntry.second)->check = check;
    }
    else {
      std::cout << "addCheck: check is the wrong type for " << lastMapEntry.first << std::endl;
    }
  }

  /** The following template addCheck member function adds a "check" lambda of function to an ArgBase_t object
  *** with the name argName.
  *** Specific versions of this function for each allowed type are built in the public section of this class
  *** Parameters:
  *** char* arg is the switch or argument name on the command line.
  *** std::function<void(ArgValue_t*)> check is a reference to a lambda that does the check on the args value.
  ***/
  template<class ArgBase_t, class ArgValue_t>
  void addCheckT(const char* argName, std::function<void(ArgValue_t*)> check) {
    auto match = argsMap.find(std::string(argName));
    ArgBase* ptr = match->second;
    if (match != argsMap.end()) {
      // parse the value for this arg and the mark it found
      if (typeid(ArgBase_t) == typeid(*ptr)) {
        static_cast<ArgBase_t*>(match->second)->check = check;
      }
      else {
        std::cout << "addCheck: check is the wrong type for " << lastMapEntry.first << std::endl;
      }
    }
    else {
      // print error if arg name is not in the argMap.
      std::cout << "addCheck does not recognise argument name " << argName << ".\n";
    }
  }

  // argsMap hold and provides the lookup of user defined arguments.
  std::map<std::string, ArgBase*> argsMap;
  // this holds the last mam entry for use by functions that add additional info the Arg* objects 
  std::pair<std::string, ArgBase*> lastMapEntry;
  std::string_view printHelpPretext = "";

public:

  /** The parse member function parses the standard C command line arguments and modifies the variables
  *** int argc is the number of args
  *** char* argv[] array of pointers to char strings holding the actual arguments
  *** return is a bool that is true if no errors found (program should continue) or
                               false if errors found or it the -h arg was present (program should be terminated) **/
  bool parse(const int argc, const char* argv[], int start = 1) {
    std::string_view argName;
    std::string_view argValue;
    // loop through the arguments
    for (int i = start; i < argc; i++) {
      std::string_view s(argv[i]);
      if (s == "-h") {
        printHelp();
        return false;
      }
      try {
        // pars out the argument name and the value text
        // first check for the equals format and split the argument
        size_t eqi = s.find('=');
        if (eqi != std::string::npos) {
          argName = s.substr(0, eqi);
          argValue = s.substr(eqi + 1);
        }
        // if not the equals case, see if the next argv is a valid value.
        else {
          argName = s;
          // first check that there is a next argv
          if (i >= argc - 1) argValue = "";
          else {
            // if so, then check to see if it contains an equals indicating its the next arg, not a value;
            std::string_view sp1(argv[i + 1]);
            if (sp1.find('=') != std::string::npos) argValue = "";
            else {
              // if no equals then check to see if it matches an arg name;
              if (argsMap.find(std::string(sp1)) != argsMap.end()) argValue = "";
              else {
                // The next argv could be a valid value to send to be processed.
                argValue = sp1;
                i = i + 1;
              }
            }
          }
        }
        // look up the  argument name or switch in the argMap 
        auto match = argsMap.find(std::string(argName));
        if (match != argsMap.end()) {
          // parse the value for this arg and the mark it found
          match->second->parse(argValue);
          match->second->found = true;
        }
        else {
          // print error if arg name is not in the argMap.
          std::cout << "Command line argument error on " << s << std::endl;
          std::cout << argName << " is not a recognised argment\n";
          printHelp(false);
          return false;
        }
      }
      catch (const std::exception& e) {
        // catch any errors found in parsing.
        std::cerr << "Command line argument error on " << s << ".\n";
        std::string newWhat = e.what();
        if (newWhat.find("stoll") != std::string::npos || newWhat.find("stod") != std::string::npos)
          newWhat = "expected a value or open parenthasis";
        std::cerr << newWhat << std::endl;
        printHelp(false);
        return false;
      }
    }
    // after parsing is done check to tee the the required argument were present.
    for (auto arg : argsMap) {
      if (arg.second->required && !arg.second->found) {
        std::cout << "Required command line argument " << arg.first << " was not present.\n";
        printHelp(false);
        return false;
      }
    }
    return true;
  }

  void addHelpPretext(std::string_view pt) {
    printHelpPretext = pt;
  }

  /** The printHelp member function prints the list of user defined arguments with type and descriptions.
  *** It is called by parse() when -h is found in the arguments or when a parse error is identified.
  *** bool printDefaults : if true prints the default value for each of the args.**/
  void printHelp(bool printDefaults = true) {
    // set up a lambda to add text strings to a specific width
    auto pw = [](std::string in, size_t w) -> std::string {
      std::string out = in;
      while (out.length() < w) out += ' ';
      return out;
      };
    int aTypeWidth = 15;
    // scan for maximum argment name (switch) width
    size_t aNameWidth = 3;
    for (auto arg : argsMap) {
      if (arg.first.size() > aNameWidth)  aNameWidth = arg.first.size();
    }
    // start by printing the pretext if available
    if (printHelpPretext != "")
      std::cout << printHelpPretext << std::endl;
    // now print the table of defined args, header first
    std::cout << "Argument format: arg value or arg=value\n";
    //std::cout << format("{:{}} : {:{}} {}\n", "arg", aNameWidth, "value type", aTypeWidth, "Description");
    std::cout << pw("arg", aNameWidth) << " : " << pw("value type", aTypeWidth) << " " << "Description" << std::endl;
    std::cout << "----------------------------------------------------------------------------------------\n";
    //std::cout << format("{:{}} : {:{}} {}\n", "-h", aNameWidth, " ", aTypeWidth, "Prints this help message");
    std::cout << pw("-h", aNameWidth) << " : " << pw(" ", aTypeWidth) << " " << "Prints this help message" << std::endl;
    // then each defined arg
    for (auto arg : argsMap) {
      //std::cout << format("{:{}} : {:{}} {}\n", arg.second->argName, aNameWidth, arg.second->type(), aTypeWidth, arg.second->description);
      std::cout << pw(arg.second->argName, aNameWidth) << " : " << pw(arg.second->type(), aTypeWidth) << " " << arg.second->description << std::endl;
      //if (printDefaults) std::cout << format("{:{}}    default = {}\n", "", aNameWidth + aTypeWidth, arg.second->toString());
      if (printDefaults) std::cout << pw("", aNameWidth + aTypeWidth) << "    default = " << arg.second->toString() << std::endl;
    }
  }

  /** The getValuesString member function returns the command line string for all the args with their current settings. **/
  std::string getValuesString() {
    std::string out;
    for (auto arg : argsMap) {
      out += std::string(arg.second->argName) + "=" + arg.second->toString() + " ";
    }
    return out;
  }

  /**  The clear member function deletes all the data stored in the args map. **/
  void clear() {
    // first delete all the arg* objects.
    for (auto arg : argsMap) {
      delete arg.second;
    }
    argsMap.clear();
  }

  // destructor
  ~ParseArgs() {
    clear();
  }
};
