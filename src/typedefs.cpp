#pragma once
#define BOOST_DISABLE_ASSERTS 

#include "boost\cstdint.hpp"

typedef float float32;
typedef double float64;

typedef int8_t int8;
typedef uint8_t uint8;

typedef int16_t int16;
typedef uint16_t uint16;



enum Action
   {
      Ignore,
      Overwrite,
      Store
   };



#include <boost/python/exception_translator.hpp>
#include <exception>

struct my_exception : std::exception
{
	std::string s;
	my_exception(std::string ss) : s(ss) {}
	const char* what() const throw() { return s.c_str(); }
};

void translate(my_exception const& e)
{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_RuntimeError, e.what());
}
