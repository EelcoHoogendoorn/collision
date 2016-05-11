
#pragma once

#include <string>
#include <exception>

#include <python.h>
#include <boost/python.hpp>
#include <boost/python/exception_translator.hpp>


struct python_exception : std::exception
{
	std::string s;
	python_exception(std::string ss) : s(ss) {}
	const char* what() const throw() { return s.c_str(); }
};

void translate(python_exception const& e)
{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_RuntimeError, e.what());
}
