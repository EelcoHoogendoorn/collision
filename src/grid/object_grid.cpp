#pragma once

#include "stdafx.h"

#include "grid_spec.cpp"
#include "sparse_grid.cpp"


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
	    if ((self.box != other.box).any())
	        throw python_exception('Grids to be intersected do not have identical specifications')
        // for each intersection of cell hashes
        // generate pairs
        // filter duplicates
	}


};
