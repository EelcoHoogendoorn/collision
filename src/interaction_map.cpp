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
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/adjacent_filtered.hpp>
//#include <boost/range/adaptors.hpp>       // somehow gives a link error?

#include "typedefs.cpp"
#include "numpy_eigen/array.cpp"
#include "numpy_boost/ndarray.cpp"
#include "numpy_boost/exception.cpp"
#include "maps.cpp"


using namespace boost;
using namespace boost::adaptors;


template<typename range_t>
auto ndarray_from_range(const range_t input) {
    typedef typename range_value<range_t>::type element_t;
    std::vector<element_t> tmp;
    for (element_t e : input)
        tmp.push_back(e);
    ndarray<element_t> output({(int32)tmp.size()});
    copy(tmp, output.begin());
    return output;
}



template<typename real_t, typename fixed_t, int NDim>
class GridSpec {
    /* helper class which stores the defining properties of a grid
    */

public:
	typedef fixed_t                         fixed_t;     // expose as public
	typedef real_t                          real_t;     // expose as public

	typedef GridSpec<real_t, fixed_t, NDim> self_t;
	typedef int32_t                         index_t;       // 32 bit int is a fine size type; 4 billion points isnt very likely
	typedef int64_t                         hash_t;

	typedef Eigen::Array<real_t, 2, NDim>	box_t;
	typedef Eigen::Array<real_t, 1, NDim>	vector_t;
	typedef Eigen::Array<fixed_t, 1, NDim>	cell_t;

	const real_t scale;    // size of a virtual voxel
	const box_t  box;      // maximum extent of pointcloud; used to map coordinates to positive integers
	const cell_t shape;    // number of virtual buckets in each direction; used to prevent out-of-bound lookup
	const cell_t strides;  // for lex-ranking cells

	GridSpec(ndarray<real_t, 2> position, real_t scale) :
		scale	(scale),
		box		(init_box(position)),
		shape	(init_shape()),
		strides	(init_strides())
//		stencil	(init_stencil(stencil))
	{
	}

private:
	//determine bounding box from point cloud positions
	auto init_box(ndarray<real_t, 2> position) const {
		real_t inf = std::numeric_limits<real_t>::infinity();
		box_t box;
		box.row(0).fill(+inf);
		box.row(1).fill(-inf);
		for (vector_t p : position.view<vector_t>()) {
			box.row(0) = box.row(0).min(p);
			box.row(1) = box.row(1).max(p);
		}
		return box;
	}
	// integer shape of the domain
	cell_t init_shape() const {      // interestingly, using auto as return type fails spectacularly
		return transform(box.row(1) - box.row(0)).cast<fixed_t>() + 1;	// use +0.5 before cast?
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

public:
	//map a global coord into the grid local coords
	inline vector_t transform(const vector_t& v) const {
		return (v - box.row(0)) / scale;
	}
	inline cell_t cell_from_local_position(const vector_t& v) const {
		return v.cast<fixed_t>();	// we want to round towards zero; surprised that we do not need a -0.5 for that..
	}
	inline cell_t cell_from_position(const vector_t& v) const {
		return cell_from_local_position(transform(v));
	}
	inline fixed_t hash_from_cell(cell_t cell) const {
		return (cell * strides).sum();
	}

	// initialize the stencil of hash offsets
	ndarray<hash_t> compute_offsets(ndarray<fixed_t, 2> stencil) const {
        auto arr = ndarray_from_range(
            stencil.view<cell_t>()
                | transformed([&](cell_t c){return hash_from_cell(c);})
                | filtered([&](hash_t h){return h > 0;})
                );
        boost::sort(arr);
        return arr;
	}

};



template<typename spec_t>
class PointGrid {
	/*
	provide spatial lookup in O(1) time for n-dimensional point clouds
	*/

public:
    typedef PointGrid<spec_t>				self_t;
    typedef typename spec_t::real_t         real_t;
	typedef typename spec_t::index_t		index_t;
	typedef typename spec_t::hash_t			hash_t;
	typedef typename spec_t::fixed_t		fixed_t;

	typedef typename spec_t::box_t			box_t;
	typedef typename spec_t::vector_t		vector_t;
	typedef typename spec_t::cell_t			cell_t;

	const spec_t				 spec;
	const ndarray<vector_t>      position;    // positions
	const index_t                n_points;    // number of points

public:
	const ndarray<fixed_t>       cell_id;     // the cell coordinates a vertex resides in
	const ndarray<index_t>       permutation; // index array mapping the vertices to lexographically sorted order
	const ndarray<index_t>       pivots;	  // boundaries between buckets of cells as viewed under permutation
	const index_t                n_buckets;   // number of cells
	const ndarray<index_t>       offsets;	  // determines stencil

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
	explicit PointGrid(spec_t spec, ndarray<real_t, 2> position) :
		spec(spec),
		position(position.view<vector_t>()),
		n_points(position.size()),
		cell_id(init_cells()),
		permutation(init_permutation()),
		pivots(init_pivots()),
		n_buckets(pivots.size() - 1),
		bucket_from_cell(       // create a map to invert the cell_from_bucket function
			boost::combine(
				irange(0, n_buckets) | transformed([&](index_t b) {return cell_from_bucket(b);}),
				irange(0, n_buckets)
			)
		)
	{
	}
	// constructor
	explicit PointGrid(spec_t spec, ndarray<real_t, 2> position, ndarray<index_t> offsets) :
		spec(spec),
		position(position.view<vector_t>()),
		n_points(position.size()),
		cell_id(init_cells()),
		permutation(init_permutation()),
		pivots(init_pivots()),
		n_buckets(pivots.size() - 1),
		offsets(offsets),
		bucket_from_cell(       // create a map to invert the cell_from_bucket function
			boost::combine(
				irange(0, n_buckets) | transformed([&](index_t b) {return cell_from_bucket(b);}),
				irange(0, n_buckets)
			)
		)
	{
	}

	// create a new pointgrid, using the permutation of existing pointgrid as initial guess
	self_t update(const ndarray<real_t, 2> position) {
	    return self_t(spec, position, permutation, offsets);
	}
	explicit PointGrid(spec_t spec, ndarray<real_t, 2> position, ndarray<index_t> permutation, ndarray<index_t> offsets) :
		spec		(spec),
		position	(position.view<vector_t>()),
		n_points	(position.size()),
		cell_id		(init_cells()),
		permutation	(init_permutation(permutation)),
		pivots		(init_pivots()),
		n_buckets	(pivots.size() - 1),
		bucket_from_cell(       // create a map to invert the cell_from_bucket function
			boost::combine(
				irange(0, n_buckets) | transformed([&](index_t b) {return cell_from_bucket(b);}),
				irange(0, n_buckets)
			)
		)
	{
	}

private:
	// determine grid cells
	auto init_cells() const {
		auto cell_id = ndarray<fixed_t>({ n_points });
		for (index_t v : irange(0, n_points))
			cell_id[v] = spec.hash_from_cell(spec.cell_from_position(position[v]));
		return cell_id;
	}
	// finds the index vector that puts the vertices in a lexographically sorted order
	auto init_permutation() const {
	    return init_permutation(irange(0, n_points));
	}
	template<typename range_t>
	auto init_permutation(const range_t initial_permutation) const {
		ndarray<index_t> _permutation({ n_points });
		// init with initial order; 0 to n
		boost::copy(initial_permutation, _permutation.begin());
		// branching-free lex sorting ftw
		auto _cell_id = cell_id.range();
		auto lex = [&](index_t l, index_t r) {return _cell_id[r] > _cell_id[l];};
        // wow, casting permutation to raw range yield factor 3 performance in sorted case
		boost::sort(_permutation.range(), lex);
		return _permutation;
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


protected:
	//convert bucket index into cell coords
	inline fixed_t cell_from_bucket(index_t b) const {
		return cell_id[permutation[pivots[b]]];
	}
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

	//loop over bounding box and apply body; this should be moved to subclass really
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

		const vector_t lmin = spec.transform(gmin);
		const vector_t lmax = spec.transform(gmax);

		//intersected volume is not positive; bail
		if ((lmin.max(vector_t(0, 0, 0)) > lmax.min(spec.shape.cast<real_t>())).any()) return;

		//compute local cell coords; constrain to [0-shape)
		const cell_t lb =  spec.cell_from_local_position(lmin).max(cell_t(0, 0, 0));
		const cell_t ub = (spec.cell_from_local_position(lmax) + 1).min(spec.shape);

		//loop over all cells that intersect with bb
		for (auto x : irange(lb[0], ub[0]))
			for (auto y : irange(lb[1], ub[1]))
				for (auto z : irange(lb[2], ub[2]))
					for (index_t v : vertices_from_cell(spec.hash_from_cell(cell_t(x, y, z))))
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

		if ((box.row(0) > spec.box.row(1)).any() || (box.row(1) < spec.box.row(0)).any()) return;

		for (auto v : irange(0, n_points))
			if (in_box(v))
				body(v);
	}

	//symmetric iteration over all vertex pairs, summing reaction forces
	template <class F>
	void for_each_vertex_pair(const F& body) const
	{
		const real_t ls2 = spec.scale*spec.scale;

		//this function wraps sign conventions regarding relative position and force
		const auto wrapper = [&](const index_t vi, const index_t vj){
			const vector_t rp = position[vi] - position[vj];
			const real_t d2 = (rp*rp).sum();
			if (d2 > ls2) return;
			body(vi, vj, d2);  // store index pair
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
			for (index_t o : offsets) {
				const auto bj = vertices_from_cell(ci + o);
				for (const index_t vj : bj)		//loop over other guy first; he might be empty, giving early exit
					for (const index_t vi : bi)
						wrapper(vi, vj);
			}
		});
	}
    // compute [n, 2] array of all paris within length_scale distance
	ndarray<index_t, 2> get_pairs() const
	{
	    typedef Eigen::Array<index_t, 1, 2> pair_t;

	    std::vector<pair_t> pairs;
	    for_each_vertex_pair([&](index_t i, index_t j, real_t d2) {
	        pairs.push_back(pair_t(i, j));
	    });
	    index_t n_pairs(pairs.size());
//	    ndarray<pair_t> _pairs({n_pairs});
//        boost::copy(pairs, _pairs.begin());
//        ndarray<index_t, 2> output = _pairs.unview<index_t>();
        ndarray<index_t, 2> output({ n_pairs, 2});
        for (index_t i : irange(0, n_pairs))
        {
            output[i][0] = pairs[i][0];
            output[i][1] = pairs[i][1];
        }
        return output;
	}

};
