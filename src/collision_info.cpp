#pragma once

#include <iostream>
#include <math.h>
#include <algorithm>
#include <tuple>

#include <boost/range/irange.hpp>

#include "typedefs.cpp"
#include "numpy_eigen/array.cpp"
#include "numpy_boost/ndarray.cpp"
#include "interaction_map.cpp"
#include "triangle_mesh.cpp"


enum Action
{
	Ignore,
	Overwrite,
	Store
};



/*
class to handle interactions between a triangle mesh and a pointcloud
make this a method of mesh?
*/
template<typename grid_t, typename mesh_t>
class CollisionInfo {

public:
	typedef int32 index_t;
	typedef float32 real_t;
	typedef Eigen::Array<real_t, 1, 3> real3;
	typedef Eigen::Array<real_t, 3, 3> real33;
	typedef Eigen::Matrix<real_t, 3, 1> _real3;
	typedef Eigen::Matrix<real_t, 3, 3> _real33;
	typedef Eigen::Array<index_t, 1, 3> triangle_t;

	const grid_t& grid;
	const mesh_t& mesh;
	CollisionInfo& info;

	// contact properties
	ndarray<real3>      bary;       // bary position of contact
	ndarray<real3>      normal;     // normal of contact
	ndarray<real_t>     depth;      // depth in triangle
	ndarray<index_t>    triangle;   // index of matching triangle

	index_t count;		//number of collisions found
	const bool self_intersect;

	// property interface
	auto get_bary()     { return bary.unview<real_t>(); }
	auto get_normal()   { return normal.unview<real_t>(); }
	auto get_depth()    { return depth; }
	auto get_triangle() { return triangle; }
	void set_bary       (ndarray<real_t, 2> bary) {}
	void set_normal     (ndarray<real_t, 2> normal) {}
	void set_depth      (ndarray<real_t>    depth) {}
	void set_triangle   (ndarray<index_t>   triangle) {}


	//construct from vertex count alone
	explicit CollisionInfo(grid_t& grid, mesh_t& mesh, const bool self_intersect) :
		grid    (grid),
		mesh    (mesh),
		info    (*this),
		bary    (ndarray<real_t, 2>({ grid.n_points, 3 }).view<real3>()),
		normal  (ndarray<real_t, 2>({ grid.n_points, 3 }).view<real3>()),
		depth   ({ grid.n_points }),
		triangle({ grid.n_points }),
		count   (-1),
		self_intersect(self_intersect)
	{
		fill(triangle, -1);     // mark all points as non-intersecting
		triangles_versus_points();
}


	/*
	compute barycentric coords of origin-centered triangle along normal direction d
	*/
	static inline _real3 area_project(const _real33& s, const _real3& d) {
		const _real3 r(
			d.dot(s.col(1).cross(s.col(2))),
			d.dot(s.col(2).cross(s.col(0))),
			d.dot(s.col(0).cross(s.col(1)))
		);
		return r / r.sum();
	}

	/*
	iterative intersection algorithm on a point versus a triangle swept along its triangle normals
	returns bary, normal and depth
	tvp are triangle vertex positions relative to origin point
	tvn are triangle vertex normals
	*/
	static auto triangle_point_test(const _real33& tvp, const _real33& tvn) {
		auto getnormal = [&](const _real3& bary)   { return tvn * bary;};					//get normal, given a bary
		auto getbary   = [&](const _real3& normal) { return area_project(tvp, normal);};	//get bary, given a normal
		auto iterate   = [&](const _real3& bary)   { return getbary(getnormal(bary));};		//one iteration is a composition of the two

		const real_t pou = 1.0 / 3;          // start in the center of the triangle
		const _real3 init(pou, pou, pou);

		const _real3 bary   = iterate(iterate(init));		// two fixed point iterations seem to suffice
		const _real3 normal = getnormal(bary).normalized();	// need normal to compute depth
		const  real_t  depth  = (tvp * bary).dot(normal);		// compute depth along normal

		return std::make_tuple(bary.transpose(), normal.transpose(), depth);
	}


	//apply body to each contact pair
	//body signature is (vertexindex, triangleindex) -> void
	template <class F>
	void for_each_contact(const F& body) const {
		for (const index_t p : boost::irange(0, vg.n_points)) {
			const index_t t = triangle[p];
			if (t == -1) continue;	//skip unpaired points
			body(p, t);
		}
	}

	//wrap bounding box iterator, adding in triangle-awareness
	template <class F>
	void for_each_vertex_in_triangle(index_t t, F& body) const {
		grid.for_each_vertex_in_bounding_box(mesh.boxes[t], body);
	}

	//wrap bounding box iterator, adding in triangle-awareness
	template <class F>
	void for_each_vertex_in_triangle_naive(index_t t, F& body) const {
		grid.for_each_vertex_in_bounding_box_naive(mesh.boxes[t], body);
	}


	/*
	iteration over triangle-vertex pairs, in order to initialize the collisioninfo object
	this encapsulates all the mutations and boilerplate
	*/
	template <class F>
	void for_each_pair(const F& body) {
		info.count = 0;			//initialize number of interactions

		for (const index_t t : boost::irange(0, mesh.n_triangles)) {		//loop over triangles
			//some mutable optimizations
			bool loaded = false;
			triangle_t tvi;
			real33 tvp, tvn;

			//loop over all vertices in bounding box of triangle;
			info.for_each_vertex_in_triangle(t, [&](const index_t v) {
				if (!loaded) {
					//read triangle into local mem matrix for efficient geometry computations
					//only needed if we hit a vertex at all
					tvi = mesh.triangles[t];
					for (const index_t i : boost::irange(0, 3)) {
						const index_t c = tvi[i];
						tvp.col(i) = mesh.position[c];
						tvn.col(i) = mesh.normal[c];
					}
					loaded = true;
				}

				//skip vertex-conflicts on self intersections
				if (info.self_intersect)
				    if ((tvi == v).any()) return;

				//call the body; this tells us what to think of this triangle-vertex pair
				const auto ret = body(t, v, tvp, tvn);
				const Action action = std::get<0>(ret);

				if (action == Action::Ignore) return;
				if (action == Action::Store) info.count += 1;			//register novel interaction

				//store the result
				std::tie(info.bary    [v],
					     info.normal  [v],
					     info.depth   [v]) = std::get<1>(ret);		//ret is a tuple of tuples...
				         info.triangle[v] = t;
			});
		}
	}


	/* collide vertexgrid and trianglemesh, building up collision info */
	void triangles_versus_points() {
		for_each_pair([&](	//loop over all nearby triangle-vertex pairs
			const index_t t, const index_t v,
			const real33& tvp, const real33& tvn					//triangle vertex positions and normals
		) {
			//intersect translated triangle, and unpack results
			const auto  intersection(triangle_point_test(tvp.colwise() - grid.position[v].transpose(), tvn));
			const real3 bary(std::get<0>(intersection));
			const real_t  d   (std::get<2>(intersection));

			//decide what to do with the intersection result?
			const Action action((
				d > +mesh.inner ||	//check intersection result for bounds
				d < -mesh.outer ||
				(bary < 0).any()) ?
				    Action::Ignore :			            // if out of bounds, ignore
				    info.triangle[v] == -1 ?		        // check for novelty
        				Action::Store :			            // if novel, store
        				d*d < info.depth[v]*info.depth[v] ? // check for superiority
            				Action::Overwrite :		        // if superior, overwrite
            				Action::Ignore);

			return std::make_tuple(action, intersection);			//return result and what to do with it
		});
	}


	void unit_test() {
		//test triangle iteration and lookup
		for (auto t : boost::irange(0, mesh.n_triangles)) {
			std::vector<int> clever;
			for_each_vertex_in_triangle(t, [&](auto v) {clever.push_back(v);});
			boost::sort(clever);

			std::vector<int> naive;
			for_each_vertex_in_triangle_naive(t, [&](auto v) {naive.push_back(v);});
			boost::sort(naive);

			if (!boost::equal(clever, naive)) {
				boost::copy(
					clever,
					std::ostream_iterator<int>(std::cout, ","));
				boost::copy(
					naive,
					std::ostream_iterator<int>(std::cout, ","));

				throw python_exception("bug in triangle iteration detected!");
			}
		}
	}

};


//// dont remember why this indicetion for the accessors was necessary..
//const VertexGridHash& get_vertexgrid  (CollisionInfo& ci){return ci.vg;}
//const TriangleMesh&   get_trianglemesh(CollisionInfo& ci){return ci.tm;}
