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

typedef GridSpec<float32, int64, 2> Spec2d;
typedef GridSpec<float32, int64, 3> Spec3d;

typedef PointGrid<Spec2d> Grid2d;
typedef PointGrid<Spec3d> Grid3d;

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
//	numpy_boost_python_register_type<Grid3d::index_t, 2>();
	numpy_boost_python_register_type<int64, 1>();
	numpy_boost_python_register_type<int64, 2>();
	numpy_boost_python_register_type<int32, 2>();
	numpy_boost_python_register_type<float32, 1>();
	numpy_boost_python_register_type<float32, 2>();


	class_<Spec2d>("Spec2d", init<ndarray<float32, 2>, float32>())
	    .def("stencil", &Spec2d::compute_offsets)
		;
	class_<Grid2d>("Grid2d", init<Spec2d, ndarray<float32, 2>, ndarray<Grid2d::index_t>>())
		.add_property("cells",      &Grid2d::get_cells,       &Grid2d::set_cells)
		.add_property("permutation",&Grid2d::get_permutation, &Grid2d::set_permutation)
		.add_property("pivots",     &Grid2d::get_pivots,      &Grid2d::set_pivots)
		.def_readonly("n_buckets",  &Grid2d::n_buckets)
//        .def("get_pairs", &Grid2d::get_pairs)
//		.def("unit_test", &Grid2::unit_test)
		;

	class_<Spec3d>("Spec3d", init<ndarray<float32, 2>, float32>())
	    .def("stencil", &Spec3d::compute_offsets)
		;
	class_<Grid3d>("Grid3d", init<Spec3d, ndarray<float32, 2>, ndarray<Grid3d::index_t>>())
		.add_property("cells",      &Grid3d::get_cells,       &Grid3d::set_cells)
		.add_property("permutation",&Grid3d::get_permutation, &Grid3d::set_permutation)
		.add_property("pivots",     &Grid3d::get_pivots,      &Grid3d::set_pivots)
		.def_readonly("n_buckets",  &Grid3d::n_buckets)
//        .def("get_pairs", &Grid3d::get_pairs)
        .def("update", &Grid3d::update)
//		.def("unit_test", &Grid3::unit_test)
		;

	class_<Mesh>("Mesh", init<ndarray<float32, 2>, ndarray<float32, 2>, ndarray<int32, 2>, float32, float32>())
		.add_property("boxes",      &Mesh::get_boxes, &Mesh::set_boxes)
		;

	class_<Info>("Info", init<Grid3d&, Mesh&, bool>())
		.add_property("depth",      &Info::get_depth,       &Info::set_depth)
		.add_property("triangle",   &Info::get_triangle,    &Info::set_triangle)
		.add_property("bary",       &Info::get_bary,        &Info::set_bary)
		.add_property("normal",     &Info::get_normal,      &Info::set_normal)
		;

	register_exception_translator<python_exception>(&translate);

}
