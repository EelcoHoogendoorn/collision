
#pragma once

#include <functional>

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/combine.hpp>

#include "ndarray.cpp"
#include "linalg.cpp"

using namespace boost;
using namespace boost::adaptors;


//used for hashing calcs
const int3 primes(73856093, 19349663, 83492791);


template<typename key_type, typename value_type>
class HashMap {

	const int n_entries;                    //number of entries in hashmap
	ndarray<2, int>          keys;      //voxel coordinates uniquely identifying a bucket; why does key_type::scalar not work?
	ndarray<1, value_type>   values;    //bucket description, or where to look in pivot array
//    typedef std::tuple<int3, int> entry_type;
//    const numpy_boost<entry_type, 1> entries;

public:
    // construct by zipping keys and values range
    template<class K, class V>
	HashMap(const K& ikeys, const V& ivalues):
	    n_entries(calc_entries(boost::distance(ivalues))),
		keys({n_entries, 3}),
		values({n_entries})
	{
		//mark grid as unoccupied
		fill(values, -1);
        for (const auto pair : boost::combine(ikeys, ivalues))
            write(boost::get<0>(pair), boost::get<1>(pair));
	}

	inline int read(const key_type& key) const
	{
	    const auto _keys = keys.view<const key_type>();
		int entry = get_hash(key);			//hash guess
		while (true)						//find the right entry
		{
			if ((_keys[entry] == key).all())
			    return values[entry];	                // we found it; this should be the most common code path
			if (values[entry] == -1) return -1;	    	// if we didnt find it yet by now, we never will
			entry = (entry + 1) & (n_entries - 1);		// circular increment
		}
	}

private:
	inline void write(const key_type& key, const value_type value)
	{
	    auto _keys = keys.view<key_type>();
		int entry = get_hash(key);				// get entry initial guess
		while (true)                            // find an empty entry
		{
			if (values[entry] == -1) break;         // found an empty entry
			entry = (entry + 1) & (n_entries - 1);	// circular increment
		}
		values[entry] = value;
		_keys[entry] = key;
	}

	inline int get_hash(const key_type& key) const
	{
		const key_type c = key * primes;
		return (c[0]^c[1]^c[2]) & (n_entries - 1);
	}

	static int calc_entries(const int size)
	{
		//calc number of entries in hashmap. hashmap should have twice the number of buckets, at mimimum.
		//absolute size is proportional to the number of vertices in the dataset, not to the space they occupy
		int entries = 64;
		while (entries < size * 2) entries <<= 1;
		return entries;
	}
};
