#pragma once

#include "stdafx.h"

#include "grid_spec.cpp"
#include "sparse_grid.cpp"
#include "point_grid.cpp"


template<typename spec_t>
class ObjectGrid {
    /*
    abstract base class
    map extended objects, such as bounding boxes, to a sparse grid
    the defining distinction with the PointGrid class is that every object may occupy multiple cells
    */
public:
    typedef ObjectGrid<spec_t>				self_t;

    typedef typename spec_t::real_t         real_t;
	typedef typename spec_t::index_t		index_t;
	typedef typename spec_t::fixed_t		fixed_t;

	typedef typename spec_t::box_t			box_t;
	typedef typename spec_t::vector_t		vector_t;
	typedef typename spec_t::cell_t			cell_t;
	typedef erow<index_t, 2>                pair_t;

	const spec_t				            spec;

public:
    // allocate as single 2xn array?
	const ndarray<fixed_t>                  cell_id;     // the cell coordinates a box resides in
	ndarray<index_t>                        object_id;   // id of box generating this grid entry
    const index_t                           n_objects;
	const SparseGrid                        grid;

public:

    // given a cell hash key, return a range of the object indices located there
    auto objects_from_key(fixed_t key) const {
        return grid.indices_from_key | transformed([&](index_t i){return object_id[i];})
    }

    // ndarray of unique pairs from vector of non-unique pairs
    auto unique_pairs(std::vector<pair_t>& pairs) const {
        auto pair_order(const pair_t& i, const pair_t& j) {
            const pair_t d = i - j;
            return d(0) * n_objects + d(1) > 0;
        }
        auto pair_not_equal(const pair_t& i, const pair_t& j) {
            return (i != j).any();
        }
        boost::sort(pairs, pair_order);
        return ndarray_from_range(pairs | adjecent_filtered(pair_not_equal)).unview<index_t>();
    }

	// self-intersection; return [n, 2] of object indices
	ndarray<index_t, 2> intersect() const {
	    std::vector<pair_t> pairs;
	    // for each cell in grid
	    for (const fixed_t c : grid.unique_keys()) {
	        // generate each object pair in cell
	        const auto objects = objects_from_key(c);
	        for (index_t i : objects)
				for (index_t j : objects)
					if (i == j)
					    break;
					else
						pairs.push_back((i < j) ? pair_t(i, j) : pair_t(j, i));
        }
        return unique_pairs(pairs);
	}

	// other-intersection, where other is an object-grid
	ndarray<index_t, 2> intersect(const self_t& other) const {
	    const self_t& self = *this;

	    if (self.spec == other.spec)
	        throw python_exception('Grids to be intersected do not have identical specifications')

        // for each intersection of cell hashes
	    std::vector<fixed_t> intersection;
        boost::range::set_intersection(
            self.cell_id.range(),
            other.cell_id.range(),
            std::back_inserter(intersection)
        );

        // generate pairs
	    std::vector<pair_t> pairs;
	    // for each cell in grid
	    for (const fixed_t c : intersection) {
	        // generate each object pair in cell
	        for (index_t i : self.objects_from_key(c))
				for (index_t j : other.objects_from_key(c))
					pairs.push_back(pair_t(i, j));
        }
        return unique_pairs(pairs);
	}

	// other-intersection, where other is a pointgrid
	std::vector<fixed_t> intersect_base(const PointGrid<spec>& other) const {
	    const self_t& self = *this;

	    if (self.spec == other.spec)
	        throw python_exception('Grids to be intersected do not have identical specifications')

        // for each intersection of cell hashes
	    std::vector<fixed_t> intersection;
        boost::range::set_intersection(
            self.cell_id.range(),
            other.cell_id.range(),
            std::back_inserter(intersection)
        );

        return intersection;
	}

};
