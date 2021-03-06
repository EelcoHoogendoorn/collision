#pragma once

#include "stdafx.h"

#include "base_grid.cpp"
#include "point_grid.cpp"


template<typename spec_t, typename sub_t>
class ObjectGrid : public BaseGrid<spec_t, sub_t> {
    /*
    map extended objects, such as bounding boxes, to a sparse grid
    the defining distinction with the PointGrid class is that every object may occupy multiple cells
    */
public:
    typedef ObjectGrid<spec_t, sub_t>		self_t;


public:
    ObjectGrid(
        const spec_t spec,
        const index_t n_objects) :
        BaseGrid(spec, n_objects)
    {
    }


    // given a cell hash key, return a range of the object indices located there
    auto objects_from_key(const fixed_t key) const {
        return self.grid.indices_from_key(key)
            | transformed([&](index_t i){return self.object_id[i];});
    }
    auto objects_from_existing_key(const fixed_t key) const {
        return self.grid.indices_from_existing_key(key)
            | transformed([&](index_t i){return self.object_id[i];});
    }


    // self-intersection; return [n, 2] of object indices
    ndarray<index_t, 2> intersect_self() const {
        std::vector<pair_t> pairs;
        for (const fixed_t c : self.grid.unique_keys()) {
            const auto obj = self.objects_from_key(c);
            for (const index_t i : obj)
                for (const index_t j : obj)
                    if (i == j)
                        break;
                    else
                        if (self.object_intersects_object(self.objects[i], self.objects[j]))
                            pairs.push_back((i < j) ? pair_t(i, j) : pair_t(j, i));
        }

        return self.unique_pairs(pairs);
    }

    // other-intersection, where other is some sub-type of object-grid
    template<typename other_t>
    ndarray<index_t, 2> intersect(const other_t& other) const {
        std::vector<pair_t> pairs;
        for (const fixed_t c : self.intersect_cells(other))
            for (index_t i : self.objects_from_existing_key(c))
                for (index_t j : other.objects_from_existing_key(c))
                    if (self.object_intersects_object(self.objects[i], other.objects[j]))
                        pairs.push_back(pair_t(i, j));
        return self.unique_pairs(pairs);
    }

    // other-intersection, where other is a point-grid
    ndarray<index_t, 2> intersect_points(const PointGrid<spec_t>& other) const {
        std::vector<pair_t> pairs;
        for (const fixed_t c : self.intersect_cells(other))
            for (const index_t i : self.objects_from_existing_key(c))
                for (const index_t j : other.grid.indices_from_existing_key(c))
                    if (self.object_intersects_point(self.objects[i], other.position[j]))
                        pairs.push_back(pair_t(i, j));
        return self.as_pair_array(pairs);
    }

};
