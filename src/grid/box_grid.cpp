
#pragma once

#include "stdafx.h"

#include "grid_spec.cpp"
#include "object_grid.cpp"


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

    inline static bool point_in_box(const vector_t& point, const box_t& box) {
        return !((point < box.row(0)).any() || (point > box.row(1)).any());
    };

    // generate box-point pairs, given a set of intersecting grid cells
    ndarray<index_t, 2> intersect(const PointGrid<spec>& other) const {
    	const self_t& self = *this;
        const auto intersection = self.intersect_base(other);

        // generate pairs
	    std::vector<pair_t> pairs;
	    // for each cell in grid
	    for (const fixed_t c : intersection) {
	        // generate each object pair in cell
	        for (index_t i : self.objects_from_key(c))
				for (index_t j : other.indices_from_existing_key(c))
				    if point_in_box(other.position[j], self.boxes[i])
    					pairs.push_back(pair_t(i, j));
        }
        return ndarray_from_range(pairs).unview<index_t>();

};