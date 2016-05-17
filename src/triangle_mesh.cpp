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

template <typename real_t>
class TriangleMesh {
public:

	typedef int32_t                     index_t;
	typedef Eigen::Array<index_t, 1, 3> triangle_t;
	typedef Eigen::Array<real_t,  1, 3> vector_t;
	typedef Eigen::Array<real_t,  2, 3> box_t;

	const index_t    n_triangles;
	const index_t    n_vertices;
	const real_t     inner;
	const real_t     outer;

	const ndarray<vector_t>      position;
	const ndarray<vector_t>      normal;
	const ndarray<triangle_t>    triangles;
	const ndarray<box_t>         boxes;


	explicit TriangleMesh(
		ndarray<real_t,  2> position,
		ndarray<real_t,  2> normal,
		ndarray<index_t, 2> triangles,
		real_t inner,
		real_t outer
	) :
		position    (position.view<vector_t>()),
		normal      (normal.view<vector_t>()),
		triangles   (triangles.view<triangle_t>()),
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
		auto inf = std::numeric_limits<real_t>::infinity();

		ndarray<box_t> boxes(ndarray<real_t, 2>({ n_triangles, 6 }).view<box_t>());

		for (auto t : boost::irange(0, n_triangles))	//loop over all triangles
		{
			box_t& box = boxes[t];
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

	auto get_boxes() { return boxes.unview<real_t>(); }
	void set_boxes(ndarray<real_t, 2> b) {  }

};
