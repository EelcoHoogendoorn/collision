#pragma once

#include <limits>
#include <iostream>
#include <functional>
#include <algorithm>

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/join.hpp>

#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/adjacent_filtered.hpp>
//#include <boost/range/adaptors.hpp>       // somehow gives a link error?


#include "../typedefs.cpp"
#include "../numpy_eigen/array.cpp"
#include "../numpy_boost/ndarray.cpp"
#include "../numpy_boost/exception.cpp"


using namespace boost;
using namespace boost::adaptors;


class ScopedGILRelease {
public:
    inline ScopedGILRelease() {
        m_thread_state = PyEval_SaveThread();
    }

    inline ~ScopedGILRelease() {
        PyEval_RestoreThread(m_thread_state);
        m_thread_state = NULL;
    }

private:
    PyThreadState* m_thread_state;
};
