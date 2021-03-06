
#pragma once

#include "stdafx.h"


template<typename key_t, typename value_t, typename index_t>
class HashMap {
    /* key-value map based on hashing of keys */

    const index_t    n_items;   // number of items
    const index_t    n_entries; // number of entries
    const index_t    mask;      // bitmask for valid range
    const index_t    prime;

    ndarray<key_t>   keys;      // voxel coordinates uniquely identifying a bucket
    ndarray<value_t> values;    // bucket description, or where to look in pivot array

public:
    // construct by zipping keys and values range
    template<typename items_range>
    HashMap(const items_range& items) :
        prime       (73856093), // other options: { 73856093, 19349663, 83492791 }
        n_items     (boost::distance(items)),
        n_entries   (init_entries()),
        mask        (n_entries - 1),
        keys        (init_keys()),
        values      (init_values())
    {
        for (auto item : items)
            write(boost::get<0>(item), boost::get<1>(item));
    }

    inline value_t operator[](const key_t key) const {
        for (index_t entry : circular_view(key)) {
            if ((keys[entry] == key))
                return values[entry];	// we found it; this should be the most common code path
            if (values[entry] == -1)
                break;				// if we didnt find it yet by now, we never will
        }
        return -1;
    }

private:
    inline void write(const key_t key, const value_t value) {
        for (index_t entry : circular_view(key))
            if (values[entry] == -1) {
                values[entry] = value;
                keys[entry] = key;
                return;
            }
    }

    auto circular_view(const key_t key) const {
        index_t entry = get_hash(key);		    // get entry initial guess
        return irange(entry, entry + n_entries)
                    | transformed([&](index_t i){return i & mask;});
    }

    inline auto get_hash(const key_t key) const {
        return key * prime;
    }

    auto init_entries() const{
        //calc number of entries in hashmap. hashmap should have twice the number of items, at mimimum.
        int entries = 64;
        while (entries < n_items * 2) entries <<= 1;
        return entries;
    }

    auto init_keys() const {
        return ndarray<key_t>( list_of(n_entries) );
    }

    auto init_values() const {
        ndarray<value_t> values( list_of(n_entries) );
        fill(values, -1); 		//mark grid as unoccupied
        return values;
    }

};


template<typename key_t, typename value_t, typename index_t>
class BaseMap {
    /* interface of our custom map types; constructable from a range of key-value tuples, and read-only value lookup by key */

//	template<typename items_range>
//	BaseMap(const items_range& items) {};
//
//	inline value_t operator[](const key_t key) const {};
};


class SortMap {
    /* key-value map based on sorting of keys */

};

class StdMap {
    /* key-value map based on std implementation */

};

class DenseMap {
    /*
    allocate ndarray<index_t>(prod(shape)), and write hashes directly to that
    more cache friendly than hasing based map, and for a fluid-in-square-box type simulation,
    more memory efficient too
    */
};