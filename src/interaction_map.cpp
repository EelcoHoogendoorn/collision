#pragma once

#include <limits>

#include <boost/range.hpp>

#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/adjacent_filtered.hpp>
//#include <boost/range/adaptors.hpp>       // somehow gives a link error?

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

template<typename cell_type>
class VertexGridHash {
    /*
    this datastructure allows for coarse/braod collision queries
    it is one of the most simple datastructures to implement and debug,
    and is chosen here with future GPU-porting in mind, as it maps well to parallel architectures
    */

public:
	const ndarray<1, const float3> position;    // vertex positions
	const int n_vertices;                       // number of vertices
	const float lengthscale;                    // size of a virtual voxel, which equals the maximum interaction range of a vertex-vertex pair query

protected:
    const float23 extents;                      // maximum extent of pointcloud; used to map coordinates to positive integers
	const cell_type size;                       // number of virtual buckets in each direction; used to prevent out-of-bound lookup

	const ndarray<1, const cell_type> cell_id;	// the cell coordinates a vertex resides in
	const ndarray<1,       int      > indices;	// index array mapping the vertices to lexographically sorted order
	      ndarray<1,       int      > pivots;	// boundaries between buckets in vertices as viewed under indices
	const int n_cells;	                        // number of cells

	const HashMap<cell_type, int> bucket_from_cell; // maps cell coordinates to bucket index

public:
	//interface methods
	ndarray<1, int> get_indices() const {return this->indices;}
	void set_indices(ndarray<1, int> indices) {int a=3;}
	ndarray<1, int> get_pivots() const {return this->pivots;}
	void set_pivots(ndarray<1, int> pivots){this->pivots = pivots;}

    // grid constructor
	explicit VertexGridHash(const ndarray<2, float> position, const float lengthscale):
		position    (position.view<const float3>()),
		n_vertices  (position.size()),
		lengthscale (lengthscale),
		extents     (init_extents()),
		size        (init_size()),
		cell_id     (init_cells()),
		indices     (init_indices()),
		pivots      ({n_vertices+1}),
		n_cells     (init_pivots()),
		bucket_from_cell(       // create a map to invert the cell_from_bucket function
		    irange(0, n_cells) | transformed([&](const int i){return cell_from_bucket(i);}),
		    irange(0, n_cells)
		)
	{
	    // empty constructor; noice
	}

protected:
	//map a global coord into the grid local coords
	inline float3 transform(const float3& v) const {
		return (v - extents.row(0)) / lengthscale;
	}
	inline cell_type cell_from_local_position(const float3& v) const {
		return (v - 0.5).cast<int>();	// defacto floor
	}
	inline cell_type cell_from_position(const float3& v) const {
		return cell_from_local_position(transform(v));
	}
	//convert bucket index into cell coords
	inline cell_type cell_from_bucket(const int b) const {
		return cell_id[indices[pivots[b]]];
	}


	//determine extents of data, return size
	float23 init_extents() const
	{
	    float23 ext;
		ext.row(0).fill(+std::numeric_limits<float>::infinity());
		ext.row(1).fill(-std::numeric_limits<float>::infinity());
		for (const float3 p: position)
		{
			ext.row(0) = ext.row(0).min(p);
			ext.row(1) = ext.row(1).max(p);
		}
		return ext;
	}

	cell_type init_size() const {
	    return transform(extents.row(1) - extents.row(0)).cast<int>() + 1;
	}

	// determine grid cells
	ndarray<1, const cell_type> init_cells() const
	{
		// silly indirection, because we cannot yet allocate custom type
	    ndarray<2, cell_type::Scalar> cidxbase({n_vertices, cell_type::SizeAtCompileTime});
	    ndarray<1, cell_type> cidx = cidxbase.view<cell_type>();
		for (const int v: irange(0, n_vertices))
			cidx[v] = cell_from_position(position[v]);
	    return cidxbase.view<const cell_type>();
	}
	//finds the index vector that puts the vertices in a lexographically sorted order
	ndarray<1, int> init_indices() const
	{
		//create index array, based on lexographical ordering
	    ndarray<1, int> idx({n_vertices});
		boost::copy(irange(0, n_vertices), idx.begin());
		boost::sort(
	        idx,
			[&](const int l, const int r) {
				const cell_type cl = cell_id[l];
				const cell_type cr = cell_id[r];
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
	int init_pivots()
	{
		int np = 0;		//number of pivots
		const auto add_pivot = [&](const int b) {pivots[np] = b; np += 1;};

        auto res = indices
                    | transformed([&](const int i){return cell_id[i];})
                    | indexed(1)
                    | adjacent_filtered([](const auto a, const auto b){return (a.value() != b.value()).any();})
                    | transformed([](const auto i){return i.index();});

		add_pivot(0);
        for (const int i : res)
			add_pivot(i);

		if (np >= n_vertices)
		    throw python_exception("every vertex is in its own cell; lengthscale probably needs to go way up");

//        int_1::extent_gen extents;
//        pivots.view(boost::extents[np - 1]);
		return np - 1;
	}


protected:
	auto indices_from_bucket(const int b) const {
		return (b == -1) ? irange(0, 0) : irange(pivots[b], pivots[b+1]);
    }
    auto indices_from_cell(const cell_type& cell) const {
        return indices_from_bucket(bucket_from_cell[cell]);
    }
	auto vertices_from_cell(const cell_type& cell) const {
		return indices_from_cell(cell)
		        | transformed([&](const int i){return indices[i];});
    }

public:
    // public traversal interface; what this class is all about
	template <class F>
	void for_each_vertex_in_cell(const cell_type& cell, const F& body) const {
		for (const int i: indices_from_cell(cell))
			body(indices[i]);
	}

	//loop over each occupied cell in the grid
	template <class F>
	void for_each_cell(const F& body) const {
		for (const int b: irange(0, n_cells))
			body(cell_from_bucket(b));
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
		const cell_type lb =  cell_from_local_position(lmin)   .max(cell_type(0,0,0));
		const cell_type ub = (cell_from_local_position(lmax)+1).min(size            );

		//loop over all cells that intersect with bb
		for (const int x: irange(lb[0],ub[0]))
			for (const int y: irange(lb[1],ub[1]))
				for (const int z: irange(lb[2],ub[2]))
					for (const int v: vertices_from_cell(cell_type(x, y, z)))
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
			return !((vp < gmin).any() || (vp > gmax).any());
		};

		if ((gmin > extents.row(1)).any() || (gmax < extents.row(0)).any()) return;

		for (const int v: irange(0, n_vertices))
			if (in_box(v))
			    body(v);
	}

};