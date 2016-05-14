#pragma once

#include <limits>

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/any_range.hpp>

#include "linalg.cpp"
#include "ndarray.cpp"


/*
class to facilitate interaction of a triangle mesh with a pointcloud
only handles precompute of bounding boxes really,
but these may be needed many times per timestep
should we iteratively deform the mesh along its normals, implementing normal calcs here would make sense too
*/

template <typename real_type>
class TriangleMesh{
public:

    typedef Eigen::Array<real_type, 2, 3> vector_type;
    typedef Eigen::Array<real_type, 2, 3> box_type;

	float_2 position;
	float_2 normal;
	int_2 incidence;

	float_2 bbmin;
	float_2 bbmax;

	const int   triangles;
	const int   vertices;
	const float thickness;

	TriangleMesh(float_2 position, float_2 normal, int_2 incidence, float thickness):
		position(position), normal(normal), incidence(incidence), thickness(thickness),
			bbmin(incidence.shape()), bbmax(incidence.shape()),
			vertices(position.shape()[0]),
			triangles(incidence.shape()[0])
	{
		boundingbox();
	}

	//compute bounding boxes for each triangle
	void boundingbox()
	{
		const auto _position	= position	.range<const float3>();
		const auto _normal		= normal	.range<const float3>();
		const auto _incidence   = incidence	.range<const int3>();

		//could precompute a temp inner/outer array; not sure it is worthwhile tho

		const auto _bbmin		= bbmin		.range<float3>();
		const auto _bbmax		= bbmax		.range<float3>();

		const float inf = std::numeric_limits<float>::infinity();
		boost::fill(_bbmin, float3(+inf, +inf, +inf));
		boost::fill(_bbmax, float3(-inf, -inf, -inf));

		for (const int t: boost::irange(0, triangles))	//loop over all triangles
		{
			const int3 I = _incidence[t];
			for (const int i: boost::irange(0,3))		//loop over all vertices incident to the triangle
			{
				const int v = I[i];

				const float3 _i = _position[v] - _normal[v] * thickness * 2;
				const float3 _o = _position[v] + _normal[v] * thickness * 2;

				_bbmin[t] = _bbmin[t].min(_i);	_bbmin[t] = _bbmin[t].min(_o);
				_bbmax[t] = _bbmax[t].max(_i);	_bbmax[t] = _bbmax[t].max(_o);
			}
		}
	}


	int_2 get_incidence(){return this->incidence;}
	void set_incidence(int_2 incidence){this->incidence = incidence;}
	float_2 get_position(){return this->position;}
	void set_position(float_2 position){this->position = position;}
	float_2 get_bbmin(){return this->bbmin;}
	void set_bbmin(float_2 bbmin){this->bbmin = bbmin;}
	float_2 get_bbmax(){return this->bbmax;}
	void set_bbmax(float_2 bbmax){this->bbmax = bbmax;}

};
