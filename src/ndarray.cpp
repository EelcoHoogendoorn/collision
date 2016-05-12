/*
defines numpy boost array types used
*/
#pragma once

#include <python.h>

#include "typedefs.cpp"
#include "numpy_boost_python.hpp"

//using namespace boost::python;
template <int N, typename T>
using ndarray = numpy_boost<T, N>;


// ugly hack apparently required to init the numpy C API
#if PY_MAJOR_VERSION >= 3
int
#else
void
#endif
init_numpy()
{
import_array();
}
