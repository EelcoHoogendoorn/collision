#pragma once

#include "stdafx.h"

#include "base_grid.cpp"


template<typename spec_t>
class PointGrid : public BaseGrid<spec_t, PointGrid<spec_t>> {
    /*
    provide spatial lookup in O(1) time for n-dimensional point clouds
    */



public:
    typedef PointGrid<spec_t>				self_t;

    const ndarray<vector_t>                 position;    // positions

public:
    const ndarray<fixed_t>                  cell_id;     // the hash of a cell a point resides in
    const SparseGrid<fixed_t, index_t>      grid;        // defines buckets
    const ndarray<index_t>                  offsets;	 // determines stencil

public:
    //interface methods
    auto get_cells()        const { return self.cell_id; }


    // constructor
    explicit PointGrid(
        const spec_t spec,
        const ndarray<real_t, 2> position) :
        BaseGrid    (spec, position.size()),
        position    (position.view<vector_t>()),
        cell_id     (init_cells()),
        grid        (cell_id)
    {
    }
    // constructor with nonzero stencil, for self intersection
    explicit PointGrid(
        const spec_t spec,
        const ndarray<real_t, 2> position,
        const ndarray<index_t> offsets) :
        BaseGrid    (spec, position.size()),
        position    (position.view<vector_t>()),
        cell_id     (init_cells()),     // grid and cells can be composed
        grid        (cell_id),
        offsets     (offsets)
    {
    }

    // create a new pointgrid, using the permutation of existing pointgrid as initial guess
    self_t update(const ndarray<real_t, 2> position) const {
        return self_t(self.spec, position, self.grid.permutation, self.offsets);
    }
    explicit PointGrid(
        const spec_t spec,
        const ndarray<real_t, 2> position,
        const ndarray<index_t> permutation,
        const ndarray<index_t> offsets) :
        BaseGrid    (spec, position.size()),
        position	(position.view<vector_t>()),
        cell_id		(init_cells()),
        grid        (cell_id, permutation),
        offsets     (offsets)
    {
    }

private:
    // determine grid cells
    auto init_cells() const {
        auto cell_id = ndarray<fixed_t>({ self.n_objects });
        for (index_t v : irange(0, self.n_objects))
            cell_id[v] = self.spec.hash_from_cell(self.spec.cell_from_position(self.position[v]));
        return cell_id;
    }


public:
    // public traversal interface; what this class is all about
    template <class F>
    void for_each_point_in_cell(fixed_t cell, const F& body) const {
        for (index_t p : self.grid.indices_from_key(cell))
            body(p);
    }

    // symmetric iteration over all point pairs
    template <class F>
    void for_each_pair(const F& body) const {
        const real_t ls2 = self.spec.scale * self.spec.scale;

        const auto wrapper = [&](const index_t i, const index_t j){
            const vector_t rp = self.position[i] - self.position[j];
            const real_t d2 = (rp*rp).sum();
            if (d2 > ls2) return;
            body(i, j, d2);
        };

        //loop over all buckets
        for (const fixed_t ci : self.grid.unique_keys()) {
            const auto bi = self.grid.indices_from_existing_key(ci);
            //interaction within bucket
            for (const index_t pi : bi)
                for (index_t pj : bi)
                    if (pi == pj) break; else
                        wrapper(pi, pj);
            //loop over all neighboring buckets
            for (const fixed_t o : offsets) {
                const auto bj = self.grid.indices_from_key(ci + o);
                for (const index_t pj : bj)		//loop over other guy first; he might be empty, giving early exit
                    for (const index_t pi : bi)
                        wrapper(pi, pj);
            }
        }
    }
    // add optimized code path, which uses sorted nature of unique cell keys;
    // for 3d, we only need to incrementally update 5 contiguous ranges around the center cell

    // compute [n, 2] array of all pairs within length_scale distance
    ndarray<index_t, 2> get_pairs() const {
        std::vector<pair_t> pairs;
        for_each_pair([&](index_t i, index_t j, real_t d2) {
            pairs.push_back(pair_t(i, j));
        });
        return self.as_pair_array(pairs);
    }

};
