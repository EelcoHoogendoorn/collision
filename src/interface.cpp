/*
boost python interface definition
*/

#include "exception.cpp"
#include "interaction_map.cpp"

typedef VertexGridHash<int3> Grid;


using namespace boost::python;
BOOST_PYTHON_MODULE(Collision)
{
	//initialize numpy support
	init_numpy();
	//register array types employed; needed to avoid runtime error
	numpy_boost_python_register_type<int, 1>();
	numpy_boost_python_register_type<int, 2>();
	numpy_boost_python_register_type<float, 1>();
	numpy_boost_python_register_type<float, 2>();
//	numpy_boost_python_register_type<int3, 1>();


	class_<Grid>("VertexGridHash", init<ndarray<2, float>, float>())
		.add_property(
		"indices", &Grid::get_indices, &Grid::set_indices)
		.add_property(
		"pivots", &Grid::get_pivots, &Grid::set_pivots)
		
//		.def_readonly("vertices", &Grid::n_vertices)
		
//		.def("unit_test", &Grid::unit_test)
		;


//	class_<TriangleMesh>("TriangleMesh", init<float_2, float_2, int_2, float>())
//		.def("boundingbox", &TriangleMesh::boundingbox)
//		.add_property(
//		"position", &TriangleMesh::get_position, &TriangleMesh::set_position)
//		.add_property(
//		"incidence", &TriangleMesh::get_incidence, &TriangleMesh::set_incidence)
//		.add_property(
//		"bbmin", &TriangleMesh::get_bbmin, &TriangleMesh::set_bbmin)
//		.add_property(
//		"bbmax", &TriangleMesh::get_bbmax, &TriangleMesh::set_bbmax)
//
//		.def_readonly("vertices", &TriangleMesh::vertices)
//		;
//
//
//	class_<CollisionInfo>("CollisionInfo", init<VertexGridHash&, TriangleMesh&, const bool>())
//		.add_property(
//		"bary", &CollisionInfo::get_bary, &CollisionInfo::set_bary)
//		.add_property(
//		"normal", &CollisionInfo::get_normal, &CollisionInfo::set_normal)
//		.add_property(
//		"depth", &CollisionInfo::get_depth, &CollisionInfo::set_depth)
//		.add_property(
//		"triangle", &CollisionInfo::get_triangle, &CollisionInfo::set_triangle)
//
//		//.def_readonly("vertexgrid", &CollisionInfo::vg)
//		//.def_readonly("trianglemesh", &CollisionInfo::tm);
//
//		.def("vertexgrid",&get_vertexgrid,return_value_policy<reference_existing_object>())
//		.def("trianglemesh",&get_trianglemesh,return_value_policy<reference_existing_object>())
//
//		.def("unit_test", &CollisionInfo::unit_test)
//		;

	
	register_exception_translator<python_exception>(&translate);

}
