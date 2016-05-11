#pragma once

#include <limits>
#include <functional>

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/combine.hpp>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/adjacent_filtered.hpp>
//#include <boost/range/adaptors.hpp>

#include "linalg.cpp"
#include "ndarray.cpp"
#include "maps.cpp"

using namespace boost;
using namespace boost::adaptors;


/*
vertex grid spatial collision map
this provides O(1) spatial lookups with no requirements at all on the point layout, 
the datastructures are easy to construct and avoid dynamic allocation or cache-unfriendly structures

since we use this for surfaces embedded in 3d space, the virtual voxel grid is quite sparse
this makes a dense grid less attractive; plus the part of the datastructure that will require random access
fits snugly in L1 cache

the only way to make this faster would be to actively reorder the input points, to exploit temporal coherency in the lexsort
*/

typedef int3 cell_type;

class VertexGridHash {
    /*
    this datastructure allows for coarse/braod collision queries
    it is one of the most simple datastructures to implement and debug,
    and is chosen here with future GPU-porting in mind, as it maps well to parallel architectures
    */

public:
	const ndarray<1, float3> position;       //vertex positions
	const int n_vertices;   //number of vertices
	const float lengthscale;//size of a virtual voxel, which equals the maximum interaction range of a vertex-vertex pair query

protected:
	float3 pmin, pmax;	    //maximum extent of pointcloud; used to map coordinates to positive integers
	const int3 size;	    //number of virtual buckets in each direction; used to prevent out-of-bound lookup

	const ndarray<1, cell_type> cell_id;	        //the cell coordinates a vertex resides in
	const ndarray<1, int      > indices;	        //index array mapping the vertices to lexographically sorted order
	int_1 pivots;		    //boundaries between buckets in vertices as viewed under indices
	const int n_cells;	    //number of cells

	const HashMap<cell_type, int> cell_map; // maps int3 cell coordinates to bucket index

public:
	//interface methods
//	float_2 get_position() const {return this->position;}
//	void set_position(float_2 position){this->position = position;}
//	int_2 get_cell_id() const {return this->cell_id;}
//	void set_cell_id(int_2 cell_id){this->cell_id = cell_id;}
	int_1 get_indices() const {return this->indices;}
	void set_indices(int_1 indices){}
	int_1 get_pivots() const {return this->pivots;}
	void set_pivots(int_1 pivots){this->pivots = pivots;}


	explicit VertexGridHash(const float_2 position, const float lengthscale):
		position(position.view<float3>()),
		n_vertices(position.size()),
		lengthscale(lengthscale),
		size(measure()),
		cell_id(compute_cells()),
		indices(indexing()),
		pivots({n_vertices+1}),
		n_cells(find_cells()),
		cell_map(
		    irange(0, n_cells) | transformed([&](const int i){return cell_from_index(i);}),
		    irange(0, n_cells)
		)
	{
	    // empty constructor; noice
	}

protected:
	//map a global coord into the grid local coords
	inline float3 transform(const float3& v) const {
		return (v - pmin) / lengthscale;
	}
	inline int3 cell_from_local_position(const float3& v) const {
		return (v - 0.5).cast<int>();	// defacto floor
	}
	inline int3 cell_from_position(const float3& v) const {
		return cell_from_local_position(transform(v));
	}
	//convert bucket index into cell coords
	int3 cell_from_index(const int b) const {
		return cell_id[indices[pivots[b]]];
	}


	//determine extents of data, return size
	int3 measure()
	{
		pmin.fill(+std::numeric_limits<float>::infinity()); 
		pmax.fill(-std::numeric_limits<float>::infinity());
		for (const float3 p: position)
		{
			pmin = pmin.min(p);
			pmax = pmax.max(p);
		}
		return transform(pmax - pmin).cast<int>() + 1;  // compute size
	}

	ndarray<1, cell_type> compute_cells() const
	{
		//determine grid cells
	    ndarray<1, cell_type> cidx({n_vertices});
		for (const int v: irange(0, n_vertices))
			cidx[v] = cell_from_position(position[v]);
	    return cidx;
	}
	//finds the index vector that puts the vertices in a lexographically sorted order
	ndarray<1, int> indexing() const
	{
		//create index array, based on lexographical ordering
	    ndarray<1, int> idx({n_vertices});
		boost::copy(irange(0, n_vertices), idx.begin());
		boost::sort(
	        idx,
			[&](const int l, const int r) {
				const int3 cl = cell_id[l];
				const int3 cr = cell_id[r];
				return 
					cl[0]!=cr[0] ?
						cl[0]<cr[0]:
						cl[1]!=cr[1] ?
							cl[1]<cr[1]:
							cl[2]<cr[2];
        });
        return idx;
	}
	//divide the sorted vertices into buckets, containing vertices in the same virtual voxel
	int find_cells()
	{
	    indexing();

		int np = 0;		//number of pivots
		const auto add_pivot = [&](const int b) {pivots[np] = b; np += 1;};
//		const auto _cell_id  = cell_id.view<const cell_type>();

		add_pivot(0);

//        auto res = indices.range<const int>()
//                    | transformed([&](const int i){return _cell_id[i];})
//                    | indexed(1)
//                    | adjacent_filtered([](auto a, auto b){return (a.value() != b.value()).any();})
//                    | transformed([](auto i){return i.index();});
//        for (const int i : res)
//			add_pivot(i);

		for (const int i: irange(1, n_vertices))
			if ((cell_id[indices[i]] != cell_id[indices[i-1]]).any())		//if different from the previous one in sorted space
				add_pivot(i);
		add_pivot(n_vertices);

		if (np == n_vertices)
		    throw python_exception("every vertex is in its own cell; lengthscale probably needs to go way up");

//        int_1::extent_gen extents;
//        pivots.view(boost::extents[np - 1]);
		return np - 1;
	}


protected:
	auto indices_from_bucket(const int b) const {
		return (b == -1) ? irange(0, 0) : irange(pivots[b], pivots[b+1]);
    }
	auto vertices_from_cell(const int3& cell) const {
		return indices_from_bucket(cell_map.read(cell))
		        | transformed([&](const int i){return indices[i];});
    }
public:
    // public traversal interface; what this class is all about
	template <class F>
	void for_each_vertex_in_cell(const int3& cell, const F& body) const {
		for (const int i: indices_from_bucket(cell_map.read(cell)))
			body(indices[i]);
	}

	//loop over each occupied cell in the grid
	template <class F>
	void for_each_cell(const F& body) const {
		for (const int b: irange(0, n_cells))
			body(cell_from_index(b));
	}

	//loop over bounding box and apply body
	template <class F>
	void for_each_vertex_in_bounding_box(const float3& gmin, const float3& gmax, const F& body) const
	{
		const auto in_box = [&](const int v)
		{
			const float3 vp = position[v];
			return !((vp < gmin).any() || (vp > gmax).any());
		};

		const float3 lmin = transform(gmin);
		const float3 lmax = transform(gmax);

		//intersected volume is not positive; bail
		if ((lmin.max(float3(0,0,0)) > lmax.min(size.cast<float>())).any()) return;				

		//compute local cell coords; constrain to [0-size)
		const int3 lb =  cell_from_local_position(lmin)   .max(int3(0,0,0));
		const int3 ub = (cell_from_local_position(lmax)+1).min(size       );

		//loop over all cells that intersect with bb
		for (const int x: irange(lb[0],ub[0]))
			for (const int y: irange(lb[1],ub[1]))
				for (const int z: irange(lb[2],ub[2]))
					for (const int v: vertices_from_cell(int3(x,y,z)))
						if (in_box(v))
						    body(v);
	}
	//for unit testing purposes
	template <class F>
	void for_each_vertex_in_bounding_box_naive(const float3& gmin, const float3& gmax, const F& body) const
	{
		const auto in_box = [&](const int v)
		{
			const float3 vp = position[v];
			return !((vp<gmin).any() || (vp>gmax).any());
		};

		if ((gmin>pmax).any() || (gmax<pmin).any()) return;

		for (const int v: irange(0, n_vertices))
			if (in_box(v))
			    body(v);
	}

};