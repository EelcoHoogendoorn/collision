#pragma once

#include <array>
#include <limits>
#include <iostream>
#include <functional>
#include <algorithm>
#include <iterator>

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/join.hpp>
#include <boost/assign.hpp>

#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/adjacent_filtered.hpp>
#include <boost/range/adaptor/uniqued.hpp>

#include <boost/range/algorithm/set_algorithm.hpp>
//#include <boost/range/adaptors.hpp>       // somehow gives a link error?


#include "typedefs.cpp"
#include "exception.cpp"

#include "numpy_eigen/array.cpp"
#include "numpy_boost/ndarray.cpp"


using namespace boost;
using namespace boost::adaptors;
using namespace boost::assign;
