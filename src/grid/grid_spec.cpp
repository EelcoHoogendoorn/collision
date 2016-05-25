#pragma once

#include "stdafx.h"


template<typename real_t, typename fixed_t, typename index_t, int n_dim>
class GridSpec {
    /* helper class which stores the defining properties of a grid
    */

public:
	typedef GridSpec<real_t, fixed_t, index_t, n_dim> self_t;

	// expose as public
	typedef fixed_t                         fixed_t;
	typedef real_t                          real_t;     
	typedef index_t                         index_t;

	typedef earray  <real_t , 2, n_dim>	    box_t;
	typedef erow    <real_t ,    n_dim>	    vector_t;
	typedef erow    <fixed_t,    n_dim>	    cell_t;

	const real_t scale;    // size of a virtual voxel
	const box_t  box;      // maximum extent of pointcloud; used to map coordinates to positive integers
	const cell_t shape;    // number of virtual buckets in each direction; used to prevent out-of-bound lookup
	const cell_t strides;  // for lex-ranking cells

	GridSpec(ndarray<real_t, 2> position, real_t scale) :
		scale	(scale),
		box		(init_box(position)),
		shape	(init_shape()),
		strides	(init_strides())
	{
	}

private:
	//determine bounding box from point cloud positions
	auto init_box(ndarray<real_t, 2> position) const {
	    return compute_bounding_box(position.view<vector_t>());
	}
	// integer shape of the domain
	cell_t init_shape() const {      // interestingly, using auto as return type fails spectacularly
		return transform(box.row(1) - box.row(0)).cast<fixed_t>() + 1;	// use +0.5 before cast?
	}
	// find strides for efficient lexsort
	auto init_strides() const {
	    return compute_strides(shape);
	}


public:
    static auto compute_strides(const cell_t shape) {
		//		boost::partial_sum(shape.cast<int>(), begin(strides), std::multiplies<int>());   // doesnt work somehow
		cell_t strides;
		strides(0) = 1;
		for (auto i : irange(1, n_dim))
			strides(i) = strides(i - 1) * shape(i - 1);
		return strides;
    }
    static auto compute_bounding_box(ndarray<vector_t> points) {
		real_t inf = std::numeric_limits<real_t>::infinity();
		box_t box;
		box.row(0).fill(+inf);
		box.row(1).fill(-inf);
		for (vector_t p : points) {
			box.row(0) = box.row(0).min(p);
			box.row(1) = box.row(1).max(p);
		}
		return box;
	}

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

	inline bool operator==(const self_t& l, const self_t& r) {
	    return (l.scale == r.scale) &
	        (l.box.row(0) == r.box.row(0)).all() &
	        (l.strides ==  r.strides).all()
	}

	// initialize the stencil of hash offsets
	// boil this info further down to contiguous stretches?
	ndarray<fixed_t> compute_offsets(ndarray<fixed_t, 2> stencil) const {
        auto arr = ndarray_from_range(
            stencil.view<cell_t>()
                | transformed([&](cell_t c){return hash_from_cell(c);})
                | filtered([&](fixed_t h){return h > 0;})
                );
        boost::sort(arr);
        return arr;
	}

};
