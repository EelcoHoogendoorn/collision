
#pragma once

#include "stdafx.h"

#include "object_grid.cpp"


template<typename spec_t>
class BoxGrid : public ObjectGrid<spec_t, BoxGrid<spec_t>> {

    /* axis-aligned bounding box grid */
public:
    typedef spec_t                          spec_t;
    typedef box_t                           object_t;

    const ndarray<box_t>                    objects;
    // allocate as single 2xn array?
	ndarray<index_t>                        object_id;   // id of box generating this grid entry
	const ndarray<fixed_t>                  cell_id;     // the cell coordinates a box resides in
	const SparseGrid<fixed_t, index_t>      grid;        // defines buckets

	auto get_object_id()  const { return self.object_id; }


	// constructor
	explicit BoxGrid(
	    const spec_t spec,
	    const ndarray<real_t, 2> boxes) :
	    ObjectGrid  (spec, boxes.size()),
		objects     (boxes.view<box_t>()),
		object_id	(ndarray<index_t>({0})),
		cell_id     (init_cells()),
        grid        (cell_id)
	{
	}

	// determine grid cells and corresponding object ids
	auto init_cells() const {
		std::vector<fixed_t> _cell_id(0);
		std::vector<index_t> _object_id(0);
		for (const index_t o : irange(0, self.n_objects))
		    for (const cell_t c: self.cells_from_box(self.objects[o])) {
		        _cell_id.push_back(self.spec.hash_from_cell(c));
		        _object_id.push_back(o);
		    }
		self.object_id = ndarray_from_range(_object_id);
		return ndarray_from_range(_cell_id);
	}

    // range of all cells in a box in world space; ndim compatible, and minimal branching.
    auto cells_from_box(const box_t& box) const {
    	const cell_t lb = self.spec.cell_from_position(box.row(0)).max(self.spec.zeros);
		const cell_t ub = self.spec.cell_from_position(box.row(1)).min(self.spec.shape) + 1;

        const cell_t shape(ub - lb);
        const cell_t strides(self.spec.compute_strides(shape));
        const cell_t prod(shape * strides); // can we replace double int division by this?
        const index_t size(prod(prod.size() - 1));

        auto mod = [](cell_t l, cell_t r) {return l - ((l / r) * r);};

		return irange(0, size)
		        | transformed([&](auto h){return lb + mod(h / strides, shape);});
    }

    inline static bool object_intersects_point(const box_t& box, const vector_t& point) {
        return !((point < box.row(0)).any() || (box.row(1) < point).any());
    }

    inline static bool object_intersects_object(const box_t& l, const box_t& r) {
        return !((l.row(1) < r.row(0)).any() || (r.row(1) < l.row(0)).any());
    }

};