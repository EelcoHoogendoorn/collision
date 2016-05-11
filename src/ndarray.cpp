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

//typedef numpy_boost<int, 1> int_1;
//typedef numpy_boost<int, 2> int_2;
//
//typedef numpy_boost<float, 1> float_1;
//typedef numpy_boost<float, 2> float_2;

//typedef numpy_boost<float3, 1> float3_1;


//typedef numpy_boost<double, 3> double_3;



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
