#include <iostream>
#include "interactions.cpp"
#include "linalg.cpp"
#include "ndarray.cpp"

/*
interface definition
boost python ftw
*/



using namespace boost::python;

BOOST_PYTHON_MODULE(collision)
{
	//initialize numpy support
	import_array();
	//register array types employed; needed to avoid runtime error
	numpy_boost_python_register_type<int, 1>();
	numpy_boost_python_register_type<int, 2>();
	numpy_boost_python_register_type<float, 1>();
	numpy_boost_python_register_type<float, 2>();

	//global functions
	def("triangles_versus_points", triangles_versus_points);
	def("membrane_membrane", membrane_membrane);
	def("membrane_cortex", membrane_cortex);
	def("cortex_cortex", cortex_cortex);
	def("cortex_cortex_density", cortex_cortex_density);
	def("cortex_cortex_double_density", cortex_cortex_double_density);
//	def("cortex_cortex_membrane_density", cortex_cortex_membrane_density);

	def("map_vertex_to_triangle", map_vertex_to_triangle);
	def("map_triangle_to_vertex", map_triangle_to_vertex);
	
	def("map_merge", map_merge);


	class_<VertexGridHash>("VertexGridHash", init<float_2, float_2, float_2, float_2, float>())

		.add_property(
		"position", &VertexGridHash::get_position, &VertexGridHash::set_position)
		.add_property(
		"cell_id", &VertexGridHash::get_cell_id, &VertexGridHash::set_cell_id)
		.add_property(
		"indices", &VertexGridHash::get_indices, &VertexGridHash::set_indices)
		.add_property(
		"pivots", &VertexGridHash::get_pivots, &VertexGridHash::set_pivots)
		
		.def_readonly("vertices", &VertexGridHash::vertices)
		
		.def("unit_test", &VertexGridHash::unit_test)
		;


	class_<TriangleMesh>("TriangleMesh", init<float_2, float_2, float_2, float_2, int_2, float>())
		.def("boundingbox", &TriangleMesh::boundingbox)

		.add_property(
		"position", &TriangleMesh::get_position, &TriangleMesh::set_position)
		.add_property(
		"incidence", &TriangleMesh::get_incidence, &TriangleMesh::set_incidence)
		.add_property(
		"bbmin", &TriangleMesh::get_bbmin, &TriangleMesh::set_bbmin)
		.add_property(
		"bbmax", &TriangleMesh::get_bbmax, &TriangleMesh::set_bbmax)

		.def_readonly("vertices", &TriangleMesh::vertices)
		;


	class_<CollisionInfo>("CollisionInfo", init<VertexGridHash&, TriangleMesh&, const bool>())
		.add_property(
		"bary", &CollisionInfo::get_bary, &CollisionInfo::set_bary)
		.add_property(
		"normal", &CollisionInfo::get_normal, &CollisionInfo::set_normal)
		.add_property(
		"depth", &CollisionInfo::get_depth, &CollisionInfo::set_depth)
		.add_property(
		"triangle", &CollisionInfo::get_triangle, &CollisionInfo::set_triangle)

		//.def_readonly("vertexgrid", &CollisionInfo::vg)
		//.def_readonly("trianglemesh", &CollisionInfo::tm);

		.def("vertexgrid",&get_vertexgrid,return_value_policy<reference_existing_object>())
		.def("trianglemesh",&get_trianglemesh,return_value_policy<reference_existing_object>())

		.def("unit_test", &CollisionInfo::unit_test)
		;

	
	register_exception_translator<my_exception>(&translate);


}
