#pragma once

#include <limits>
#include <iostream>
#include <algorithm>

#include <boost/range.hpp>


#include "typedefs.cpp"
#include "numpy_eigen/array.cpp"
#include "numpy_boost/ndarray.cpp"
//#include "interaction_map.cpp"


/*
class to facilitate interaction of a triangle mesh with a pointcloud
only handles precompute of bounding boxes really,
but these may be needed many times per timestep
should we iteratively deform the mesh along its normals, implementing normal calcs here would make sense too
*/

template <typename real_type>
class TriangleMesh {
public:

	typedef int32 index_type;
	typedef Eigen::Array<index_type, 1, 3> triangle_type;
	typedef Eigen::Array<real_type, 1, 3> vector_type;
	typedef Eigen::Array<real_type, 2, 3> box_type;

	const index_type    n_triangles;
	const index_type    n_vertices;
	const real_type     inner;
	const real_type     outer;

	const ndarray<vector_type>      position;
	const ndarray<vector_type>      normal;
	const ndarray<triangle_type>    triangles;
	const ndarray<box_type>         boxes;


	explicit TriangleMesh(
		ndarray<real_type,  2> position,
		ndarray<real_type,  2> normal,
		ndarray<index_type, 2> triangles,
		real_type inner,
		real_type outer
	) :
		position    (position.view<vector_type>()),
		normal      (normal.view<vector_type>()),
		triangles   (triangles.view<triangle_type>()),
		inner       (inner),
		outer       (outer),
		n_vertices  (position.size()),
		n_triangles (triangles.size()),
		boxes       (init_boxes())
	{
	}

	//compute bounding box for each triangle
	auto init_boxes() const
	{
		auto inf = std::numeric_limits<real_type>::infinity();

		ndarray<box_type> boxes(ndarray<real_type, 2>({ n_triangles, 6 }).view<box_type>());

		for (auto t : boost::irange(0, n_triangles))	//loop over all triangles
		{
			box_type& box = boxes[t];
			box.row(0).fill(+inf);
			box.row(1).fill(-inf);

			auto triangle = triangles[t];
			for (auto i : boost::irange(0, 3))		//loop over all vertices incident to the triangle
			{
				auto v = triangle(i);

				auto i = position[v] - normal[v] * inner;
				auto o = position[v] + normal[v] * outer;

				box.row(0) = box.row(0).min(i).min(o);
				box.row(1) = box.row(1).max(i).max(o);
			}
		}
		return boxes;
	}

    // map the vertices of this mesh to a point grid; extent can be reused
//    Grid3d get_point_grid() const {
//        return Grid3d(vertices);
//    }

	auto get_boxes() { return boxes.unview<real_type>(); }
	void set_boxes(ndarray<real_type, 2> b) {  }

};
