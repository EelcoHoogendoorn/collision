#include <iostream>
#include "test.cpp"
#include "linalg.cpp"
#include "ndarray.cpp"
#include <python.h>

/*
boost python interface definition
*/

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

using namespace boost::python;
BOOST_PYTHON_MODULE(collision)
{
	//initialize numpy support
	init_numpy();
	//register array types employed; needed to avoid runtime error
	numpy_boost_python_register_type<int, 2>();

	//global functions
	def("main", main);


}
