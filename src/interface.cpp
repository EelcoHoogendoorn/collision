/*
boost python interface definition
*/

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "exception.cpp"
#include "interaction_map.cpp"
#include "typedefs.cpp"


typedef VertexGridHash<int16, float32, 2> Grid2d;
typedef VertexGridHash<int16, float32, 3> Grid3d;


using namespace boost::python;
BOOST_PYTHON_MODULE(Collision)
{
	//initialize numpy support
	init_numpy();
	
	//register array types employed; needed to avoid runtime error
	numpy_boost_python_register_type<int16, 1>();
	numpy_boost_python_register_type<int16, 2>();
	numpy_boost_python_register_type<int32, 1>();
	numpy_boost_python_register_type<int32, 2>();
	numpy_boost_python_register_type<float32, 1>();
	numpy_boost_python_register_type<float32, 2>();


	class_<Grid2d>("Grid2d", init<ndarray<2, float32>, float32>())
//		.add_property("cell_ids",   &Grid2d::get_cell_ids,    &Grid2d::set_cell_ids)
		.add_property("permutation",&Grid2d::get_permutation, &Grid2d::set_permutation)
		.add_property("pivots",     &Grid2d::get_pivots,      &Grid2d::set_pivots)
		.def_readonly("n_buckets",  &Grid2d::n_buckets)
//		.def("unit_test", &Grid2::unit_test)
		;

	class_<Grid3d>("Grid3d", init<ndarray<2, float32>, float32>())
//		.add_property("cell_ids",   &Grid3d::get_cell_ids,    &Grid3d::set_cell_ids)
		.add_property("permutation",&Grid3d::get_permutation, &Grid3d::set_permutation)
		.add_property("pivots",     &Grid3d::get_pivots,      &Grid3d::set_pivots)
		.def_readonly("n_buckets",  &Grid3d::n_buckets)
//		.def("unit_test", &Grid3::unit_test)
		;

	register_exception_translator<python_exception>(&translate);

}
