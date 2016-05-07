/*
defines numpy boost array types used
*/
#pragma once

#include "typedefs.cpp"
#include "numpy_boost_python.hpp"

//using namespace boost::python;

typedef numpy_boost<int, 1> int_1;
typedef numpy_boost<int, 2> int_2;

typedef numpy_boost<float, 1> float_1;
typedef numpy_boost<float, 2> float_2;

//typedef numpy_boost<float3, 1> float3_1;


typedef numpy_boost<double, 3> double_3;

