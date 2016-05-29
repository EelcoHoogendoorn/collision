/*
boost python interface definition
*/

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#define BOOST_DISABLE_ASSERTS

#include "typedefs.cpp"
#include "exception.cpp"

#include "grid/grid_spec.cpp"
#include "grid/point_grid.cpp"
#include "grid/box_grid.cpp"
//#include "triangle_mesh.cpp"
//#include "collision_info.cpp"

using namespace boost;

typedef float32_t   real_t;         // real type of coordinate space
typedef int64_t     fixed_t;        // representation of cell coordinates; 32 bits would give only 10 bits in 3d
typedef int32_t     index_t;        // indexing type; we never expect to have more than 4 billion of anything

typedef GridSpec<real_t, fixed_t, index_t, 2> Spec2d;
typedef GridSpec<real_t, fixed_t, index_t, 3> Spec3d;

typedef PointGrid<Spec2d> PointGrid2d;
typedef PointGrid<Spec3d> PointGrid3d;

typedef BoxGrid<Spec2d> BoxGrid2d;
typedef BoxGrid<Spec3d> BoxGrid3d;

//typedef TriangleMesh<float32> Mesh;
//typedef CollisionInfo<Grid3d, Mesh> Info;



using namespace boost::python;
BOOST_PYTHON_MODULE(Collision)
{
    // init GIL control
    PyEval_InitThreads();

	// initialize numpy support
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
	    .def("stencil",             &Spec2d::compute_offsets)
		;
	class_<PointGrid2d>("PointGrid2d", init<Spec2d, ndarray<real_t, 2>, ndarray<index_t>>())
		.def("cells",               &PointGrid2d::get_cells)
		.def("permutation",         &PointGrid2d::get_permutation)
        .def("intersect_self",      &PointGrid2d::get_pairs)
        .def("unique_keys",         &PointGrid2d::get_unique_keys)
		;
	class_<BoxGrid2d>("BoxGrid2d", init<Spec2d, ndarray<real_t, 2>>())
		.def("permutation",		    &BoxGrid2d::get_permutation)
		.def("object_id",		    &BoxGrid2d::get_object_id)
		.def("intersect_self",	    &BoxGrid2d::intersect_self)
		.def("intersect_points",    &BoxGrid2d::intersect_points)
        .def("unique_keys",         &BoxGrid2d::get_unique_keys)
		;

	class_<Spec3d>("Spec3d", init<ndarray<real_t, 2>, real_t>())
	    .def("stencil",             &Spec3d::compute_offsets)
		;
	class_<PointGrid3d>("PointGrid3d", init<Spec3d, ndarray<real_t, 2>, ndarray<index_t>>())
		.def("cells",               &PointGrid3d::get_cells)
		.def("permutation",         &PointGrid3d::get_permutation)
        .def("intersect_self",      &PointGrid3d::get_pairs)
        .def("update",              &PointGrid3d::update)
		;
	class_<BoxGrid3d>("BoxGrid3d", init<Spec3d, ndarray<real_t, 2>>())
		.def("permutation",		    &BoxGrid3d::get_permutation)
		.def("object_id",		    &BoxGrid3d::get_object_id)
		.def("intersect_self",	    &BoxGrid3d::intersect_self)
		.def("intersect_points",    &BoxGrid3d::intersect_points)

//        .def("intersect_points",   &BoxGrid3d::intersect_points)
		;

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
