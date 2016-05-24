#pragma once

#include <functional>
#include <algorithm>

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/algorithm.hpp>

#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/adjacent_filtered.hpp>
//#include <boost/range/adaptors.hpp>       // somehow gives a link error?

#include "typedefs.cpp"
#include "numpy_eigen/array.cpp"
#include "numpy_boost/ndarray.cpp"
#include "numpy_boost/exception.cpp"
#include "maps.cpp"


using namespace boost;
using namespace boost::adaptors;


/*

*/

// no need for spec really...
template<typename spec_t>
class SparseGrid {
    /* this is some kind of multimap. basically just grouping some sortable thing by permutation order */

public:
	typedef typename spec_t::index_t		index_t;
	typedef typename spec_t::key_t			key_t;

public:
	const spec_t				 spec;
    const ndarray<key_t>         keys;
    const index_t                n_keys;

	const ndarray<index_t>       permutation; // index array mapping the keys to lexographically sorted order
	const ndarray<index_t>       pivots;	  // boundaries between buckets of keys as viewed under permutation
	const index_t                n_buckets;   // number of unique keys
	const HashMap<key_t, index_t, index_t> bucket_from_key; // maps key to bucket index

	auto get_permutation()  const { return permutation; }
	auto get_pivots()       const { return pivots; }
	void set_permutation    (ndarray<index_t> permutation)  {}
	void set_pivots         (ndarray<index_t> pivots)       {}

public:
	// constructor
	explicit SparseGrid(spec_t spec, ndarray<key_t> keys) :
		spec(spec),
        keys(keys),
        n_keys(keys.size()),
		permutation(init_permutation()),
		pivots(init_pivots()),
		n_buckets(pivots.size() - 1),
		bucket_from_key(       // create a map to invert the key_from_bucket function
			boost::combine(
				irange(0, n_buckets) | transformed([&](index_t b) {return key_from_bucket(b);}),
				irange(0, n_buckets)
			)
		)
	{
	}

    // construct using permutation initial guess
	explicit SparseGrid(spec_t spec, ndarray<key_t> keys, ndarray<index_t> permutation) :
		spec(spec),
        keys(keys),
        n_keys(keys.size()),
		permutation(init_permutation(permutation)),
		pivots(init_pivots()),
		n_buckets(pivots.size() - 1),
		bucket_from_key(       // create a map to invert the key_from_bucket function
			boost::combine(
				irange(0, n_buckets) | transformed([&](index_t b) {return key_from_bucket(b);}),
				irange(0, n_buckets)
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
		add_pivot(n_points);

		return pivots.resize(n_pivots);
	}


public:
	// get n-th unique key
	inline key_t key_from_bucket(index_t b) const {
		return keys[permutation[pivots[b]]];
	}
	inline auto indices_from_bucket(index_t b) const {
		return (b == -1) ? irange(0, 0) : irange(pivots[b], pivots[b + 1]);
	}
	inline auto indices_from_key(key_t key) const {
		return indices_from_bucket(bucket_from_key[key])
			| transformed([&](index_t i) {return permutation[i];});
	}

};
