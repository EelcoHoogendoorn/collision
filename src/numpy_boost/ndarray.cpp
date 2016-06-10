/*
defines numpy boost array types used
*/
#pragma once

#include <vector>

#include <boost/array.hpp>
#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>

#include <python.h>

#include "numpy_boost_python.hpp"


// make ndarray subclass of numpy_boost?
//template <typename T, int N=1>
//using ndarray = numpy_boost<T, N>;
// macro is backwardsw compatible with older version of cpp
#define ndarray numpy_boost


// extension method for ndarray; should be a static method on class, really
// copy range of known size to ndarray of the same type
template<typename range_t>
ndarray<typename boost::range_value<range_t>::type>
ndarray_from_range(const range_t input) {
    typedef typename boost::range_value<range_t>::type element_t;
    const boost::array<int, 1> shape = {{ boost::distance(input) }};
    ndarray<element_t> output(shape);
    boost::copy(input, output.begin());
    return output;
}

// copy iterable range of unknown size to ndarray of the same type, via a vector temporary
template<typename range_t>
ndarray<typename boost::range_value<range_t>::type>
ndarray_from_iterable(const range_t input) {
    typedef typename boost::range_value<range_t>::type element_t;
    std::vector<element_t> tmp;
    for (const element_t e : input)
        tmp.push_back(e);
    return ndarray_from_range(tmp);
}


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
