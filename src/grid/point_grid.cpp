#pragma once

#include <functional>
#include <algorithm>

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>

#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
//#include <boost/range/adaptors.hpp>       // somehow gives a link error?

#include "typedefs.cpp"
#include "numpy_eigen/array.cpp"
#include "numpy_boost/ndarray.cpp"
#include "numpy_boost/exception.cpp"


using namespace boost;
using namespace boost::adaptors;


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
	const ndarray<fixed_t>       cell_id;     // the hash of a cell a point resides in
	const ndarray<index_t>       offsets;	  // determines stencil
	const SparseGrid<spec_t>     grid;        // defines buckets

public:
	//interface methods
	auto get_cells()        const { return cell_id; }
	void set_cells          (ndarray<fixed_t> cells)        {}


	// constructor
	explicit PointGrid(spec_t spec, ndarray<real_t, 2> position) :
		spec        (spec),
		position    (position.view<vector_t>()),
		n_points    (position.size()),
		cell_id     (init_cells()),
		grid        (spec, cell_id)
	{
	}
	// constructor with nonzero stencil, for self intersection
	explicit PointGrid(spec_t spec, ndarray<real_t, 2> position, ndarray<index_t> offsets) :
		spec        (spec),
		position    (position.view<vector_t>()),
		n_points    (position.size()),
		cell_id     (init_cells()),
		offsets     (offsets),
		grid        (spec, cell_id)
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
		offsets     (offsets),
		grid        (spec, cell_id, permutation)
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


public:
	// initialize the stencil of hash offsets
	ndarray<hash_t> compute_offsets(ndarray<fixed_t, 2> stencil) const {
        auto arr = ndarray_from_range(
            stencil.view<cell_t>()
                | transformed([&](cell_t c){return spec.hash_from_cell(c);})
                | filtered([&](hash_t h){return h > 0;})
                );
        boost::sort(arr);
        return arr;
	}


	// public traversal interface; what this class is all about
	template <class F>
	void for_each_point_in_cell(fixed_t cell, const F& body) const {
		for (index_t p : grid.objects_from_key(cell))
			body(p);
	}

	//loop over each occupied cell in the grid
	template <class F>
	void for_each_cell(const F& body) const {
		for (index_t b : irange(0, n_buckets))
			body(grid.key_from_bucket(b));
	}


	// symmetric iteration over all point pairs
	template <class F>
	void for_each_pair(const F& body) const
	{
		const real_t ls2 = spec.scale*spec.scale;

		const auto wrapper = [&](const index_t i, const index_t j){
			const vector_t rp = position[i] - position[j];
			const real_t d2 = (rp*rp).sum();
			if (d2 > ls2) return;
			body(i, j, d2);
		};

		//loop over all buckets
		for_each_cell([&](const fixed_t ci) {
            const auto bi = grid.objects_from_key(ci);
			//interaction within bucket
			for (index_t pi : bi)
				for (index_t pj : bi)
					if (pi==pj) break; else
						wrapper(pi, pj);
			//loop over all neighboring buckets
			for (index_t o : offsets) {
				const auto bj = grid.objects_from_key(ci + o);
				for (const index_t pj : bj)		//loop over other guy first; he might be empty, giving early exit
					for (const index_t pi : bi)
						wrapper(pi, pj);
			}
		});
	}
    // compute [n, 2] array of all pairs within length_scale distance
	ndarray<index_t, 2> get_pairs() const
	{
	    typedef erow<index_t, 2> pair_t;

	    std::vector<pair_t> pairs;
	    for_each_pair([&](index_t i, index_t j, real_t d2) {
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

	// self intersection
	auto intersect() const {
	}

	auto intersect(const PointGrid& other) const {
	}

	auto intersect(const ObjectGrid& other) const {
	}
};


