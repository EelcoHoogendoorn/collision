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


#include "../typedefs.cpp"
#include "../numpy_eigen/array.cpp"
#include "../numpy_boost/ndarray.cpp"
#include "../numpy_boost/exception.cpp"

#include "grid_spec.cpp"
#include "sparse_grid.cpp"



using namespace boost;
using namespace boost::adaptors;


template<typename spec_t>
class ObjectGrid {
    /*
    extend pointgrid with object labels
    every box generates a set of cell-ids
    */
public:
    typedef ObjectGrid<spec_t>				self_t;
    typedef typename spec_t::real_t         real_t;
	typedef typename spec_t::index_t		index_t;
	typedef typename spec_t::fixed_t		fixed_t;

	typedef typename spec_t::box_t			box_t;
	typedef typename spec_t::vector_t		vector_t;
	typedef typename spec_t::cell_t			cell_t;

	const spec_t				            spec;

public:
    // allocate as single 2xn array?
	const ndarray<fixed_t>                  cell_id;     // the cell coordinates a box resides in
	ndarray<index_t>                        object_id;   // id of box generating this grid entry

	const SparseGrid                        grid;

public:


	// self-intersection
	ndarray<index_t, 2> intersect() const {
	    // for each cell in grid
	    // generate each object pair in cell
        // filter duplicates
	}

	// other-intersection
	ndarray<index_t, 2> intersect(self_t& other) const {
	    self_t& self = *this;
	    // assert that specs are identical
        // for each intersection of cell hashes
        // generate pairs
        // filter duplicates
	}


};


template<typename spec_t>
class BoxGrid : public ObjectGrid<spec_t> {

    const ndarray<box_t>                    boxes;
	const index_t                           n_boxes;    // number of boxes

	// constructor
	explicit BoxGrid(spec_t spec, ndarray<real_t, 3> boxes) :
		spec(spec),
		boxes(boxes.view<box_t>()),
		n_boxes(boxes.size()),
		cell_id(init_cells()),
        grid(cell_id)
	{
	}

	// determine grid cells and corresponding object ids
	auto init_cells() const {
		std::vector<fixed_t> cell_id(n_boxes);
		std::vector<index_t> object_id(n_boxes);
		for (index_t b : irange(0, n_boxes))
		    for (cell_t c: cells_from_box(boxes[b])) {
		        cell_id.push_back(spec.hash_from_cell(c));
		        object_id.push_back(b);
		    }
		object_id = ndarray_from_range(object_id);
		return ndarray_from_range(cell_id);
	}

    // range of all cells in a box in world space; ndim compatible, and minimal branching.
    auto cells_from_box(const box_t& box) const {
    	const cell_t lb = spec.cell_from_position(box.row(0)).max(cell_t(0, 0, 0));
		const cell_t ub = spec.cell_from_position(box.row(1)).min(spec.shape) + 1;

        const cell_t shape = ub - lb;
        const cell_t strides = spec.compute_strides(shape);
        const index_t size = strides(NDim - 1) * shape(NDim - 1);

		return irange(0, size)
		    | transformed([&](auto h){return lb + (h % strides);});
    }


	template <class F>
	void for_each_point_in_box(const box_t& box, const F& body) const
	{
		const auto in_box = [&](auto p)
		{
			const vector_t& vp = position[p];
			return !((vp < box.row(0)).any() || (vp > box.row(1)).any());
		};

	    for (cell_t c : cells_from_box(box))
            for (index_t p : grid.indices_from_key(spec.hash_from_cell(c)))
                if (in_box(p))
                    body(p);
	}
};