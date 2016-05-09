#pragma once

#include <boost/range/irange.hpp>
#include <math.h>
#include <algorithm>
#include <tuple>

#include "linalg.cpp"
#include "ndarray.cpp"
#include "interaction_map.cpp"
#include "triangle_mesh.cpp"


/*
compute barycentric coords of origin-centered triangle along normal direction d
*/
inline _float3 area_project(const _float33& s, const _float3& d)
{
	const _float3 r(
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
std::tuple<const float3, const float3, const float> triangle_point_test(
	const _float33& tvp, const _float33& tvn)
{
	const auto getnormal = [&](const _float3& bary)  { return tvn * bary;};						//get normal, given a bary
	const auto getbary   = [&](const _float3& normal){ return area_project(tvp, normal);};		//get bary, given a normal
	const auto iterate   = [&](const _float3& bary)  { return getbary(getnormal(bary));};		//one iteration is a composition of the two

	const  float  pou = 1.0/3;
	const _float3 init(pou,pou,pou);

	const _float3 bary   = iterate(iterate(init));			//two fixed point iterations seem to suffice
	const _float3 normal = getnormal(bary).normalized();	//need normal to compute depth
	const  float  depth  = (tvp * bary).dot(normal);		//compute depth along normal

	return std::make_tuple(bary, normal, depth);
}


/*
class to handle interactions between a triangle mesh and a pointcloud
*/
class CollisionInfo {

public:
	VertexGridHash& vg;
	TriangleMesh& tm;

	//precomputed collision properties
	float_2 bary;
	float_2 normal;
	float_1 depth;
	int_1 triangle;

	int count;		//number of collisions found
	const bool self_intersect;

    // property interface
	float_2 get_bary(){return this->bary;}
	void set_bary(float_2 bary){this->bary = bary;}
	float_2 get_normal(){return this->normal;}
	void set_normal(float_2 normal){this->normal = normal;}
	float_1 get_depth(){return this->depth;}
	void set_depth(float_1 depth){this->depth = depth;}
	int_1 get_triangle(){return this->triangle;}
	void set_triangle(int_1 triangle){this->triangle = triangle;}


	//construct from vertex count alone
	CollisionInfo(VertexGridHash& vg, TriangleMesh& tm, const bool self_intersect):
		vg(vg), tm(tm),
		bary    ({vg.vertices, 3}),
		normal  ({vg.vertices, 3}),
		depth   ({vg.vertices}),
		triangle({vg.vertices}),
		count   (-1),
		self_intersect(self_intersect)
	{
		fill(triangle, -1);
	}

	//apply body to each contact pair
	//body signature is (vertexindex, triangleindex) -> void
	template <class F>
	void for_each_contact(const F& body)
	{
		for (const int v: boost::irange(0, vg.vertices))
		{
			const int t = triangle[v];
			if (t==-1) continue;	//skip unpaired vertices
			body(v, t);
		}
	}

	//wrap bounding box iterator, adding in triangle-awareness
	template <class F>
	void for_each_vertex_in_triangle(const int t, const F& body)
	{
		vg.for_each_vertex_in_bounding_box(
			tm.bbmin.range<const float3>()[t],
			tm.bbmax.range<const float3>()[t],
			[&](const int v){body(v);}
		);
	}
	//wrap bounding box iterator, adding in triangle-awareness
	template <class F>
	void for_each_vertex_in_triangle_naive(const int t, const F& body)
	{
		vg.for_each_vertex_in_bounding_box_naive(
			tm.bbmin.range<const float3>()[t],
			tm.bbmax.range<const float3>()[t],
			[&](const int v){body(v);}
		);
	}

	/*
	iteration over triangle-vertex pairs, in order to initialize the collisioninfo object
	this encapsulates all the mutations and boilerplate
	*/
	template <class F>
	void for_each_pair(const F& body)
	{
		CollisionInfo& ci = *this;
		auto& tm = ci.tm;

		auto tm_incidence	= tm.incidence	.range<const int3>();;
		auto tm_position	= tm.position	.range<const float3>();
		auto tm_normal		= tm.normal		.range<const float3>();

		auto ci_bary		= ci.bary		.range<float3>();		//mutable output arrays
		auto ci_normal		= ci.normal		.range<float3>();
		auto ci_depth		= ci.depth		.range<float>();
		auto ci_triangle	= ci.triangle	.range<int>();

		ci.count = 0;			//initialize number of interactions

		for (const int t: boost::irange(0, tm.triangles))		//loop over triangles
		{
			//some mutable optimizations
			bool loaded = false;
			int3 tvi;
			float33 tvp, tvn;

			//loop over all vertices in bounding box of triangle;
			ci.for_each_vertex_in_triangle(t, [&](const int v)
			{
				if (!loaded)
				{
					//read triangle into local mem matrix for efficient geometry computations
					//only needed if we hit a vertex at all
					tvi = tm_incidence[t];
					for (const int i: boost::irange(0,3))
					{
						const int v = tvi[i];
						tvp.col(i) = tm_position[v];
						tvn.col(i) = tm_normal[v];
					}
					loaded = true;
				}

				//skip vertex-conflicts on self intersections
				if (ci.self_intersect) if ((tvi==v).any()) return;

				//call the body; this tells us what to think of this triangle-vertex pair
				const auto ret = body(t, v, tvp, tvn);
				const Action action = std::get<0>(ret);

				if (action==Action::Ignore) return;
				if (action==Action::Store) ci.count += 1;			//register novel interaction

				//store the result
				ci_triangle[v] = t;
				std::tie(ci_bary  [v],
						 ci_normal[v],
						 ci_depth [v]) = std::get<1>(ret);		//ret is a tuple of tuples...
			});
		}
	}


    /* collide vertexgrid and trianglemesh */
    void triangles_versus_points()
    {
        auto vg_position	= vg.position	.range<const float3>();

        auto ci_depth		= depth		    .range<const float>();
        auto ci_triangle	= triangle	    .range<const int>();

        for_each_pair([&]	//loop over all nearby triangle-vertex pairs
        (
            const int t, const int v,
            const float33& tvp, const float33& tvn					//triangle vertex positions and normals
        )
        {
            //intersect translated triangle, and unpack results
            const auto   intersection = triangle_point_test(tvp.colwise()-vg_position[v], tvn);
            const float3 bary		  = std::get<0>(intersection);
            const float  depth        = std::get<2>(intersection);

            //decide what to do with the intersection result?
            const Action action =
                depth > +2*tm.thickness ||							//check intersection result for bounds
                depth < -2*tm.thickness ||
                (bary < 0).any() ?
                    Action::Ignore:									//if out of bounds, ignore
                    ci_triangle[v] == -1 ?							//check for novelty
                        Action::Store :								//if novel, store
                        depth*depth < ci_depth[v]*ci_depth[v] ?		//check for superiority
                            Action::Overwrite:						//if superior, overwrite
                            Action::Ignore;

            return std::make_tuple(action, intersection);			//return result and what to do with it
        });
    }


	void unit_test()
	{
		//test triangle iteration and lookup
		for (const int t: boost::irange(0, tm.triangles))
		{
			std::vector<int> clever;
			for_each_vertex_in_triangle      (t, [&](const int v){clever.push_back(v);});
			boost::sort(clever);

			std::vector<int> naive;
			for_each_vertex_in_triangle_naive(t, [&](const int v){naive.push_back(v);});
			boost::sort(naive);

			if (!boost::equal(clever, naive))
			{
				boost::copy(
					clever,
					std::ostream_iterator<int>(std::cout, ","));
				boost::copy(
					naive,
					std::ostream_iterator<int>(std::cout, ","));

				throw my_exception("bug in triangle iteration detected!");
			}
		}
	}

};


// dont remember why this indicetion for the accessors was necessary..
const VertexGridHash& get_vertexgrid  (CollisionInfo& ci){return ci.vg;}
const TriangleMesh&   get_trianglemesh(CollisionInfo& ci){return ci.tm;}
