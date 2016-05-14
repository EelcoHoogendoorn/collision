#pragma once

#include <math>
#include <algorithm>
#include <tuple>

#include <boost/range/irange.hpp>

#include "interaction_map.cpp"
#include "triangle_mesh.cpp"
#include "collision_info.cpp"


/*
this file contains all functions that are directly relevant to the behavior of the simulation
exactly how we find all nearby vertices and triangles in linear time is interesting, 
but an implementation detail as far as the bio-physics/mathematics is concerned

this module strives for maximum 'correctness', without sacrificing any performance
mutable state is for rookies ;)
*/


/*
add contact forces between membranes
*/
void membrane_membrane(
	CollisionInfo& ci, 
	const float_1 modulation,
	const float scaling, 
	const float damping)
{
	const auto vg_normal		= ci.vg.normal.range<const float3>();

	ci.for_each_contact_react([&](  //loop over each colliding triangle-vertex pair, returning force vector
		const int v, const int t, 
		const float depth, const float3& normal, const float3& rv)		
	{
		const float d = depth / ci.tm.thickness;
		//quadratic with two roots of positive curvature
		const auto quadratic = [d](const float l, const float r) {return (d-l)*(d-r);};	

		//force curve
		const float curve = 
			ci.self_intersect ?
				d > 0 ?
					-quadratic(2,3):		//internal repulsion
					d > -1 ?
						quadratic(-1,-2):	//regular behavior
						0					//no self-stickiness. may not be so bad with bigger spacing though
				:							//else if no self intersection
				d > 0 ?
					0:						//shouldnt come this far. but maybe we should keep forcing if we do anyway?
					quadratic(-1,-2);		//regular behavior

		const float angle = 1;//std::max(0.0f, -(vg_normal[v]*normal).sum());
		const float mod = modulation[v];

		return								//return force vector between this vertex and triangle
			normal * angle * scaling * mod * curve
			- rv * damping * mod;
	});
}

/*
membrane-cortex interaction
quite a cheap function; cost is in ci
*/
void membrane_cortex(
	CollisionInfo& ci, 
	const float scaling,
	const float damping)
{
	ci.for_each_contact_react([&](
		const int v, const int t,
		const float depth, const float3& normal, const float3& rv)
	{
		const float d = depth / ci.tm.thickness;
		const float curve = 1+abs(d);
		const float rvn = (rv*normal).sum();
		return normal * 
			(curve * depth * scaling 
			- rvn*damping);
	});
}


/*
interactions between cortical fluid particles 
*/
void cortex_cortex(
	VertexGridHash& vg, 
	const float scale, 
	const float damping)
{
	const float ls2 = vg.lengthscale*vg.lengthscale;

	vg.for_each_vertex_pair([&](
		const int i, const int j, 
		const float3& rp, const float3& rv)
	{
		const float d2 = (rp*rp).sum() + 1e-12;
		const float d = sqrt(d2);
		const float3 normal = rp/d;
		return 
			normal * (sqrt(ls2)-d) * scale 
			- rv * damping;
	});
}

/*
translate pressure field into force contributions between particles
*/
void cortex_cortex_density(
	VertexGridHash& vg, 
	const float_1 density, 
	const float scale, 
	const float damping)
{
	const float ls2 = vg.lengthscale*vg.lengthscale;

	vg.for_each_vertex_pair([&](
		const int i, const int j, 
		const float3& rp, const float3& rv)
	{
		const float d2 = (rp*rp).sum() + 1e-12;

		const float d = sqrt(d2);
		const float dd = d/vg.lengthscale;
		const float3 normal = rp/d;

		const float _density = (density[i]+density[j])/2;
		const float curve = (1-dd*dd);

		return
			  normal * curve * _density * scale 
			- rv * curve * damping;
	});
}
void cortex_cortex_double_density(
	VertexGridHash& vg, 
	const float_1 near_density, 
	const float_1 far_density, 
	const float scale, 
	const float damping)
{
	const float ls2 = vg.lengthscale*vg.lengthscale;

	vg.for_each_vertex_pair([&](
		const int i, const int j, 
		const float3& rp, const float3& rv)
	{
		const float d2 = (rp*rp).sum() + 1e-12;

		const float d = sqrt(d2);
		const float dd = d/vg.lengthscale;
		const float3 normal = rp/d;

		const float q = 1-dd;

		const float _near_density = (near_density[i]+near_density[j])/2;
		const float _far_density  = (far_density [i]+far_density [j])/2;
		return
			  normal * (_near_density * q + _far_density * q*q) * scale
			- rv * q * damping;
	});
}

/*
put in seperate module?
map from vertices to triangles' vertices; typical use case; map from cortex to membrane vertices
do we need metrics anywhere?
*/
void map_vertex_to_triangle(CollisionInfo& ci, float_1 vg_density, float_1 tm_density)
{
	auto& vg = ci.vg;
	auto& tm = ci.tm;

	const auto ci_bary		= ci.bary		.range<const float3>();
	const auto tm_incidence	= tm.incidence	.range<const int3>();

	ci.for_each_contact([&](const int v, const int t)
	{
		const int3   tvi  = tm_incidence[t];
		const float3 bary = ci_bary[v];
		for (const int i: boost::irange(0,3))
			tm_density[tvi[i]] +=  vg_density[v] * bary[i];
	});
}


/*
distribution projected on vertices is linear interpolation over the triangles
that is, we measure a concentration, not an amount
*/
void map_triangle_to_vertex(CollisionInfo& ci, float_1 tm_density, float_1 vg_density)
{
	auto& vg = ci.vg;
	auto& tm = ci.tm;

	const auto ci_bary		= ci.bary		.range<const float3>();
	const auto tm_incidence	= tm.incidence	.range<const int3>();

	ci.for_each_contact([&](const int v, const int t)
	{
		const int3 tvi = tm_incidence[t];
		const float3 bary = ci_bary[v];

		for (const int i: boost::irange(0,3))
			vg_density[v] += tm_density[tvi[i]] * bary[i];
	});
}

/*
function to calc permutation operator that performs cortex merging
*/
int map_merge(
	VertexGridHash& vg,
	int_1 row_id,
	int_1 col_id,
	const float limit)
{
	//allocate temp pair array
	int_1 pair_id(vg.vertices);
	boost::fill(pair_id, -1);
	int counter = 0;				//ugh; just when you think mutable state was banished form this file
	auto pair = [&](const int i, const int j){pair_id[i] = pair_id[j] = counter++;};	//pair up two vertices
	auto paired = [&](const int v){return pair_id[v]!=-1;};							//predicate to check pairing state

	const float ls2 = vg.lengthscale*vg.lengthscale;

	//look for pairs
	vg.for_each_vertex_pair([&](
		const int i, const int j, 
		const float3& rp, const float3& rv)
	{
		const float d2 = (rp*rp).sum() + 1e-12;
		const float d = sqrt(d2);
		if (d<limit*vg.lengthscale)
			if (!paired(i) && !paired(j))
				pair(i,j);

		return float3(0,0,0);		//use non-reactant looping instead?
	});

	//read availability of support from ci here?
	for (const int v: boost::irange(0,vg.vertices))
	{
		col_id[v] = v;
		row_id[v] = paired(v) ? pair_id[v] : counter++;
	}
	return counter;		//number of remapped vertices
}