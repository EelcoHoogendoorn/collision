/*
boost python interface definition
*/

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#define BOOST_DISABLE_ASSERTS

#include "typedefs.cpp"
#include "numpy_boost/exception.cpp"
#include "grid/grid_spec.cpp"
#include "grid/point_grid.cpp"
//#include "triangle_mesh.cpp"
//#include "collision_info.cpp"

using namespace boost;

typedef float32_t   real_t;         // real type of coordinate space
typedef int64_t     fixed_t;        // representation of cell coordinates; 32 bits would likely suffice
typedef int32_t     index_t;        // indexing type; we never expect to have more than 4 billion of anything

typedef GridSpec<real_t, fixed_t, index_t, 2> Spec2d;
typedef GridSpec<real_t, fixed_t, index_t, 3> Spec3d;

typedef PointGrid<Spec2d> Grid2d;
typedef PointGrid<Spec3d> Grid3d;

//typedef TriangleMesh<float32> Mesh;
//
//typedef CollisionInfo<Grid3d, Mesh> Info;



using namespace boost::python;
BOOST_PYTHON_MODULE(Collision)
{
	//initialize numpy support
	init_numpy();
	
	//register array types employed; needed to avoid runtime error
	numpy_boost_python_register_type<int16_t, 1>();
	numpy_boost_python_register_type<int16_t, 2>();
	numpy_boost_python_register_type<int32_t, 1>();
	numpy_boost_python_register_type<int32_t, 2>();
	numpy_boost_python_register_type<int64_t, 1>();
	numpy_boost_python_register_type<int64_t, 2>();

	numpy_boost_python_register_type<real_t, 1>();
	numpy_boost_python_register_type<real_t, 2>();
	numpy_boost_python_register_type<real_t, 3>();


	class_<Spec2d>("Spec2d", init<ndarray<real_t, 2>, real_t>())
//	    .def("stencil", &Spec2d::compute_offsets)
		;
	class_<Grid2d>("Grid2d", init<Spec2d, ndarray<real_t, 2>, ndarray<index_t>>())
		.def("cells",       &Grid2d::get_cells)
		.def("permutation", &Grid2d::get_permutation)
        .def("pairs",       &Grid2d::get_pairs)
//		.add_property("pivots",     &Grid2d::get_pivots,      &Grid2d::set_pivots)
//		.def_readonly("n_buckets",  &Grid2d::n_buckets)
//	    .def("stencil", &Spec2d::compute_offsets)
//		.def("unit_test", &Grid2::unit_test)
		;

//	class_<Spec3d>("Spec3d", init<ndarray<float32, 2>, float32>())
//	    .def("stencil", &Spec3d::compute_offsets)
//		;
//	class_<Grid3d>("Grid3d", init<Spec3d, ndarray<float32, 2>, ndarray<Grid3d::index_t>>())
//		.add_property("cells",      &Grid3d::get_cells,       &Grid3d::set_cells)
//		.add_property("permutation",&Grid3d::get_permutation, &Grid3d::set_permutation)
//		.add_property("pivots",     &Grid3d::get_pivots,      &Grid3d::set_pivots)
//		.def_readonly("n_buckets",  &Grid3d::n_buckets)
//        .def("get_pairs", &Grid3d::get_pairs)
//        .def("update", &Grid3d::update)
////		.def("unit_test", &Grid3::unit_test)
//		;

//	class_<Mesh>("Mesh", init<ndarray<float32, 2>, ndarray<float32, 2>, ndarray<int32, 2>, float32, float32>())
//		.add_property("boxes",      &Mesh::get_boxes, &Mesh::set_boxes)
//		;
//
//	class_<Info>("Info", init<Grid3d&, Mesh&, bool>())
//		.add_property("depth",      &Info::get_depth,       &Info::set_depth)
//		.add_property("triangle",   &Info::get_triangle,    &Info::set_triangle)
//		.add_property("bary",       &Info::get_bary,        &Info::set_bary)
//		.add_property("normal",     &Info::get_normal,      &Info::set_normal)
//		;

	register_exception_translator<python_exception>(&translate);

}
