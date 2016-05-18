#pragma once

#include <limits>
#include <iostream>
#include <functional>
#include <algorithm>

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>

#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/adjacent_filtered.hpp>
//#include <boost/range/adaptors.hpp>       // somehow gives a link error?

#include "typedefs.cpp"
#include "numpy_eigen/array.cpp"
#include "numpy_boost/ndarray.cpp"
#include "numpy_boost/exception.cpp"
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

template<typename real_t, typename fixed_t, int NDim>
class PointGrid {
	/*
	this datastructure allows for coarse/braod collision queries
	it is one of the most simple datastructures to implement and debug,
	and is chosen here with future GPU-porting in mind, as it maps well to parallel architectures
	*/

public:
    typedef real_t                          real_t;     // expose as public
	typedef int32_t                         index_t;       // 32 bit int is a fine size type; 4 billion points isnt very likely
	typedef int64_t                         hash_t;

	typedef Eigen::Array<real_t,  2, NDim>	box_t;
	typedef Eigen::Array<real_t,  1, NDim>	vector_t;
	typedef Eigen::Array<fixed_t, 1, NDim>	cell_t;

	const ndarray<vector_t>      position;    // positions
	const index_t                n_points;    // number of points
	const real_t                 lengthscale; // size of a virtual voxel

	const box_t                  extents;     // maximum extent of pointcloud; used to map coordinates to positive integers
	const cell_t                 shape;       // number of virtual buckets in each direction; used to prevent out-of-bound lookup
	const cell_t                 strides;     // for lex-ranking cells
public:
	const ndarray<fixed_t>       cell_id;     // the cell coordinates a vertex resides in
	const ndarray<index_t>       permutation; // index array mapping the vertices to lexographically sorted order
	const ndarray<index_t>       pivots;	  // boundaries between buckets of cells as viewed under permutation
	const index_t                n_buckets;   // number of cells
	const ndarray<fixed_t>       stencil;	  // relative jumps to perform on grid to scan the entire stencil

	const HashMap<fixed_t, index_t, index_t> bucket_from_cell; // maps cell coordinates to bucket index

public:
	//interface methods
	auto get_cells()        const { return cell_id; }
	auto get_permutation()  const { return permutation; }
	auto get_pivots()       const { return pivots; }
	void set_cells          (ndarray<fixed_t> cells)        {}
	void set_permutation    (ndarray<index_t> permutation)  {}
	void set_pivots         (ndarray<index_t> pivots)       {}

	// constructor
	explicit PointGrid(ndarray<real_t, 2> position, real_t lengthscale) :
		position	(position.view<vector_t>()),
		n_points	(position.size()),
		lengthscale	(lengthscale),
		extents		(init_extents()),
		shape		(init_shape()),
		strides		(init_strides()),
		cell_id		(init_cells()),
		permutation	(init_permutation()),
		pivots		(init_pivots()),
		n_buckets	(pivots.size() - 1),
		bucket_from_cell(       // create a map to invert the cell_from_bucket function
			boost::combine(
				irange(0, n_buckets) | transformed([&](index_t b) {return cell_from_bucket(b);}),
				irange(0, n_buckets)
			)
		)
	{   // empty constructor; noice
	}

	explicit PointGrid(ndarray<real_t, 2> position, real_t lengthscale, ndarray<index_t, 2> stencil) :
		position	(position.view<vector_t>()),
		n_points	(position.size()),
		lengthscale	(lengthscale),
		extents		(init_extents()),
		shape		(init_shape()),
		strides		(init_strides()),
		cell_id		(init_cells()),
		permutation	(init_permutation()),
		pivots		(init_pivots()),
		n_buckets	(pivots.size() - 1),
		stencil     (init_stencil),
		bucket_from_cell(       // create a map to invert the cell_from_bucket function
			boost::combine(
				irange(0, n_buckets) | transformed([&](index_t b) {return cell_from_bucket(b);}),
				irange(0, n_buckets)
			)
		)
	{   // empty constructor; noice
	}

private:
	//determine extents of data
	auto init_extents() const {
		real_t inf = std::numeric_limits<real_t>::infinity();
		box_t extents;
		extents.row(0).fill(+inf);
		extents.row(1).fill(-inf);
		for (vector_t p : position) {
			extents.row(0) = extents.row(0).min(p);
			extents.row(1) = extents.row(1).max(p);
		}
		return extents;
	}
	// integer shape of the domain
	cell_t init_shape() const {      // interestingly, using auto as return type fails spectacularly
		return transform(extents.row(1) - extents.row(0)).cast<fixed_t>() + 1;	// use +0.5 before cast?
	}
	// find strides for efficient lexsort
	auto init_strides() const {
		//		boost::partial_sum(shape.cast<int>(), begin(strides), std::multiplies<int>());   // doesnt work somehow
		cell_t strides;
		strides(0) = 1;
		for (auto i : irange(1, NDim))
			strides(i) = strides(i - 1) * shape(i - 1);
		return strides;
	}
	// determine grid cells
	auto init_cells() const {
		auto cell_id = ndarray<fixed_t>({ n_points });
		for (index_t v : irange(0, n_points))
			cell_id[v] = hash_from_cell(cell_from_position(position[v]));
		return cell_id;
	}
	// finds the index vector that puts the vertices in a lexographically sorted order
	auto init_permutation() const {
		ndarray<index_t> permutation({ n_points });
		// init with initial order; 0 to n
		boost::copy(irange(0, n_points), permutation.begin());
		// branching-free lex sorting ftw
		auto _cell_id = cell_id.range();
		auto lex = [&](index_t l, index_t r) {return _cell_id[r] > _cell_id[l];};
        // wow, casting permutation to raw range yield factor 3 performance in sorted case
		boost::sort(permutation.range(), lex);
		return permutation;
	}
	//divide the sorted vertices into buckets, containing vertices in the same virtual voxel
	auto init_pivots() const {
		// allocate array of size n_points, becuase it plays nicely with the rest of our numpy mempool
		ndarray<index_t> pivots({ n_points });

		index_t n_pivots = 0;		//number of pivots
		auto add_pivot = [&](index_t p) {pivots[n_pivots++] = p;};

		auto res = permutation
			| transformed([&](auto i) {return cell_id[i];})
			| indexed(0)
			| adjacent_filtered([](auto a, auto b) {return a.value() != b.value();})
			| transformed([](auto i) {return i.index();});

		for (index_t p : res)
			add_pivot(p);
		if (n_pivots == n_points)
			throw python_exception("every vertex is in its own cell; lengthscale probably needs to go way up");
		add_pivot(n_points);

		return pivots.resize(n_pivots);
	}
	// initialize the stencil of hash offsets
	auto init_stencil(ndarray<index_t, 2> stencil) const {
	    offsets = ndarray<cell_t>({stencil.size()});
        index_t n_offsets = 0;

	    auto push_back = [&](auto hash){offsets[n_offsets++] = hash;};

	    for (index_t o : irange(0, stencil.size())) {
	        auto hash = hash_from_cell(stencil[o]);
	        if (hash > 0) push_back(hash);
	    }
	    return offsets.resize(n_offsets);
	}


protected:
	//map a global coord into the grid local coords
	inline vector_t transform(const vector_t& v) const {
		return (v - extents.row(0)) / lengthscale;
	}
	inline cell_t cell_from_local_position(const vector_t& v) const {
		return v.cast<fixed_t>();	// we want to round towards zero; surprised that we do not need a -0.5 for that..
	}
	inline cell_t cell_from_position(const vector_t& v) const {
		return cell_from_local_position(transform(v));
	}
	//convert bucket index into cell coords
	inline fixed_t cell_from_bucket(index_t b) const {
		return cell_id[permutation[pivots[b]]];
	}
	inline fixed_t hash_from_cell(cell_t cell) const {
	    return (cell * strides).sum();
	}

protected:
	auto indices_from_bucket(index_t b) const {
		return (b == -1) ? irange(0, 0) : irange(pivots[b], pivots[b + 1]);
	}
	auto indices_from_cell(fixed_t cell) const {
		return indices_from_bucket(bucket_from_cell[cell]);
	}
	auto vertices_from_cell(fixed_t cell) const {
		return indices_from_cell(cell)
			| transformed([&](index_t i) {return permutation[i];});
	}

public:
	// public traversal interface; what this class is all about
	template <class F>
	void for_each_vertex_in_cell(fixed_t cell, const F& body) const {
		for (index_t v : vertices_from_cell(cell))
			body(v);
	}

	//loop over each occupied cell in the grid
	template <class F>
	void for_each_cell(const F& body) const {
		for (index_t b : irange(0, n_buckets))
			body(cell_from_bucket(b));
	}

	//loop over bounding box and apply body
	template <class F>
	void for_each_vertex_in_bounding_box(const box_t& box, const F& body) const
	{
	    const vector_t gmin = box.row(0);
	    const vector_t gmax = box.row(1);

		const auto in_box = [&](index_t v)
		{
			const vector_t& vp = position[v];
			return !((vp < gmin).any() || (vp > gmax).any());
		};

		const vector_t lmin = transform(gmin);
		const vector_t lmax = transform(gmax);

		//intersected volume is not positive; bail
		if ((lmin.max(vector_t(0, 0, 0)) > lmax.min(shape.cast<real_t>())).any()) return;

		//compute local cell coords; constrain to [0-shape)
		const cell_t lb =  cell_from_local_position(lmin).max(cell_t(0, 0, 0));
		const cell_t ub = (cell_from_local_position(lmax) + 1).min(shape);

		//loop over all cells that intersect with bb
		for (auto x : irange(lb[0], ub[0]))
			for (auto y : irange(lb[1], ub[1]))
				for (auto z : irange(lb[2], ub[2]))
					for (index_t v : vertices_from_cell(hash_from_cell(cell_t(x, y, z))))
						if (in_box(v))
							body(v);
	}
	//for unit testing purposes
	template <class F>
	void for_each_vertex_in_bounding_box_naive(const box_t& box, const F& body) const
	{
		const auto in_box = [&](auto v)
		{
			const vector_t& vp = position[v];
			return !((vp < box.row(0)).any() || (vp > box.row(1)).any());
		};

		if ((box.row(0) > extents.row(1)).any() || (box.row(1) < extents.row(0)).any()) return;

		for (auto v : irange(0, n_points))
			if (in_box(v))
				body(v);
	}

	//symmetric iteration over all vertex pairs, summing reaction forces
	template <class F>
	void for_each_vertex_pair(const F& body) const
	{
		const real_t ls2 = lengthscale*lengthscale;

		//this function wraps sign conventions regarding relative position and force
		const auto wrapper = [&](const index_t vi, const index_t vj){
			const vector_t rp = position[vi] - position[vj];
			const real_t d2 = (rp*rp).sum();
			if (d2 > ls2) return;
			body(vi,vj,d2);  // store index pair
		};

		//loop over all buckets
		for_each_cell([&](const fixed_t ci) {
            const auto bi = vertices_from_cell(ci);
			//interaction within bucket
			for (index_t vi : bi)
				for (index_t vj : bi)
					if (vi==vj) break; else
						wrapper(vi, vj);
			//loop over all neighboring buckets
			for (index_t o : stencil) {
				const auto bj = vertices_from_cell(ci + o);
				for (const index_t vj : bj)		//loop over other guy first; he might be empty, giving early exit
					for (const index_t vi : bi)
						wrapper(vi, vj);
			}
		});
	}
    // compute [n, 2] array of all paris within length_scale distance
	auto get_pairs() const
	{
	    typedef Eigen::Array<index_t, 1, 2> pair_t;
	    std::vector<pair_t> pairs;
	    for_each_vertex_pair([&](index_t i, index_t j, real_t d2) {
	        pairs.push_back(pair_t(i, j));
	    });
	    index_t n_pairs(pairs.size());
	    ndarray<pair_t> _pairs({n_pairs});
        boost::copy(pairs, _pairs.begin());
        return _pairs.unview<index_t>();
	}
};
