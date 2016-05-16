/*
boost python interface definition
*/

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#define BOOST_DISABLE_ASSERTS

#include "typedefs.cpp"
#include "numpy_boost/exception.cpp"
#include "interaction_map.cpp"
#include "triangle_mesh.cpp"
#include "collision_info.cpp"

typedef PointGrid<float32, int16, 2> Grid2d;
typedef PointGrid<float32, int16, 3> Grid3d;

typedef TriangleMesh<float32> Mesh;

typedef CollisionInfo<Grid3d, Mesh> Info;



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


	class_<Grid2d>("Grid2d", init<ndarray<float32, 2>, float32>())
		.add_property("cells",      &Grid2d::get_cells,       &Grid2d::set_cells)
		.add_property("permutation",&Grid2d::get_permutation, &Grid2d::set_permutation)
		.add_property("pivots",     &Grid2d::get_pivots,      &Grid2d::set_pivots)
		.def_readonly("n_buckets",  &Grid2d::n_buckets)
//		.def("unit_test", &Grid2::unit_test)
		;

	class_<Grid3d>("Grid3d", init<ndarray<float32, 2>, float32>())
		.add_property("cells",      &Grid3d::get_cells,       &Grid3d::set_cells)
		.add_property("permutation",&Grid3d::get_permutation, &Grid3d::set_permutation)
		.add_property("pivots",     &Grid3d::get_pivots,      &Grid3d::set_pivots)
		.def_readonly("n_buckets",  &Grid3d::n_buckets)
//		.def("unit_test", &Grid3::unit_test)
		;

	class_<Mesh>("Mesh", init<ndarray<float32, 2>, ndarray<float32, 2>, ndarray<int32, 2>, float32>())
		.add_property("boxes",      &Mesh::get_boxes, &Mesh::set_boxes)
		;

	class_<Info>("Info", init<Grid3d&, Mesh&, bool>())
		.add_property("depth",      &Info::get_depth, &Info::set_depth)
		.add_property("triangle",   &Info::get_triangle, &Info::set_triangle)
		.add_property("bary",       &Info::get_triangle, &Info::set_triangle)
		.add_property("normal",     &Info::get_triangle, &Info::set_triangle)
		;

	register_exception_translator<python_exception>(&translate);

}
