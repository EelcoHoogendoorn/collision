#pragma once

#include "stdafx.h"
#include "maps.cpp"


template<typename key_t, typename index_t>
class SparseGrid {
    /*
    this class groups a set of keys using an indirect sort
    it allows for querying of unique keys,
    and getting all permutation indices for a given key
    */

public:
	typedef key_t	    key_t;
	typedef index_t		index_t;

public:
    const ndarray<key_t>         keys;
    const index_t                n_keys;

	const ndarray<index_t>       permutation; // index array mapping the keys to lexographically sorted order
	const ndarray<index_t>       pivots;	  // boundaries between buckets of keys as viewed under permutation
	const index_t                n_groups;    // number of unique keys

	const HashMap<key_t, index_t, index_t> group_from_key; // maps key to group index

	auto get_permutation()  const { return permutation; }
	auto get_pivots()       const { return pivots; }

public:
	// constructor
	explicit SparseGrid(ndarray<key_t> keys) :
        keys        (keys),
        n_keys      (keys.size()),
		permutation (init_permutation()),
		pivots      (init_pivots()),
		n_groups    (pivots.size() - 1),
		group_from_key(       // create a map to invert the key_from_group function
			boost::combine(
				unique_keys(),
				irange(0, n_groups)
			)
		)
	{
	}

    // construct using permutation initial guess
	explicit SparseGrid(ndarray<key_t> keys, ndarray<index_t> permutation) :
        keys        (keys),
        n_keys      (keys.size()),
		permutation (init_permutation(permutation)),
		pivots      (init_pivots()),
		n_groups    (pivots.size() - 1),
		group_from_key(       // create a map to invert the key_from_group function
			boost::combine(
				unique_keys(),
				irange(0, n_groups)
			)
		)
	{
	}

private:
	// finds the index vector that puts the vertices in a lexographically sorted order
	auto init_permutation() const {
	    return init_permutation(irange(0, n_keys));
	}
	template<typename range_t>
	auto init_permutation(const range_t initial_permutation) const {
		ndarray<index_t> _permutation({ n_keys });
		// init with initial order; 0 to n
		boost::copy(initial_permutation, _permutation.begin());
		// branching-free lex sorting ftw
		auto _keys = keys.range();
		auto lex = [&](index_t l, index_t r) {return _keys[r] > _keys[l];};
        // wow, casting permutation to raw range yield factor 3 performance in sorted case
		boost::sort(_permutation.range(), lex);
		return _permutation;
	}
	//divide the sorted vertices into buckets, containing vertices in the same virtual voxel
	auto init_pivots() const {
		// allocate array of size n_keys, becuase it plays nicely with the rest of our numpy mempool
		// changes this into push-back into growing ndarray instead
		ndarray<index_t> pivots({ n_keys+1 });

		index_t n_pivots = 0;		//number of pivots
		auto add_pivot = [&](index_t p) {pivots[n_pivots++] = p;};

		auto res = permutation
			| transformed([&](auto i) {return keys[i];})
			| indexed(0)
			| adjacent_filtered([](auto a, auto b) {return a.value() != b.value();})
			| transformed([](auto i) {return i.index();});

		for (index_t p : res)
			add_pivot(p);
		add_pivot(n_keys);

		return pivots.resize(n_pivots);
	}


	inline auto indices_from_group(index_t g) const {
		return (g == -1) ? irange(0, 0) : irange(pivots[g], pivots[g + 1]);
	}
	// get n-th unique key
	inline key_t key_from_group(index_t g) const {
		return keys[permutation[pivots[g]]];
	}
public:
	// range over each unique key in the grid, in sorted order
	auto unique_keys() const {
		return irange(0, n_groups)
		    | transformed([&](auto g){return key_from_group(g);});
	}
    // return a range of the permutation indices within a key-group
	inline auto indices_from_key(key_t key) const {
		return indices_from_group(group_from_key[key])
			| transformed([&](index_t i) {return permutation[i];});
	}

};
