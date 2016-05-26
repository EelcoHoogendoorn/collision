#pragma once

#include "stdafx.h"

#include "grid_spec.cpp"
#include "sparse_grid.cpp"


template<typename spec_t>
class BaseGrid {
    /*
    abstract base class
    map extended objects, such as bounding boxes, to a sparse grid
    the defining distinction with the PointGrid class is that every object may occupy multiple cells
    */
public:
    typedef BaseGrid<spec_t>				self_t;
    typedef typename spec_t::real_t         real_t;
	typedef typename spec_t::index_t		index_t;
	typedef typename spec_t::fixed_t		fixed_t;

	typedef typename spec_t::box_t			box_t;
	typedef typename spec_t::vector_t		vector_t;
	typedef typename spec_t::cell_t			cell_t;
	typedef erow<index_t, 2>                pair_t;

	const spec_t				            spec;

public:
	const ndarray<fixed_t>                  cell_id;
	const SparseGrid                        grid;

public:

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

	// intersect two sparse grids, to get shared occupied cells
	std::vector<fixed_t> intersect_cells(const self_t& other) const {
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
