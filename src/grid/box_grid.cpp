
#pragma once

#include "stdafx.h"

#include "grid_spec.cpp"
#include "object_grid.cpp"


template<typename spec_t>
class BoxGrid : public ObjectGrid<spec_t, BoxGrid> {

    typedef box_t                           object_t;
    const ndarray<box_t>                    objects;
	const index_t                           n_objects;    // number of boxes

	// constructor
	explicit BoxGrid(spec_t spec, ndarray<real_t, 3> boxes) :
		spec(spec),
		objects(boxes.view<box_t>()),
		n_objects(objects.size()),
		cell_id(init_cells()),
        grid(cell_id)
	{
	}

	// determine grid cells and corresponding object ids
	auto init_cells() const {
		std::vector<fixed_t> cell_id(n_objects);
		std::vector<index_t> object_id(n_objects);
		for (index_t b : irange(0, n_objects))
		    for (cell_t c: cells_from_box(objects[b])) {
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
        const index_t size = strides(n_dim - 1) * shape(n_dim - 1);

		return irange(0, size)
		        | transformed([&](auto h){return lb + ((h / strides) % shape);});
    }

    inline static bool object_intersects_point(const box_t& box, const vector_t& point) {
        return !((point < box.row(0)).any() || (box.row(1) < point).any());
    };

    inline static bool object_intersects_object(const box_t& l, const box_t& r) {
        return !(l.row(1) < r.row(0)).any() || (r.row(1) < l.row(0)).any());
    };


};