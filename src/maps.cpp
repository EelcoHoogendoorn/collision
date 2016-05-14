
#pragma once

#include <functional>
#include <array>

#include <boost/range/irange.hpp>
#include <boost/range/difference_type.hpp>
#include <boost/tuple/tuple.hpp>

#include "ndarray.cpp"
#include "linalg.cpp"


//used for hashing calcs
const std::array<int, 3> PRIMES = { 73856093, 19349663, 83492791 };


template<class key_type, class value_type, int NDim>
class HashMap {
	typename typedef key_type::Scalar key_type_scalar;
	typedef RowArray<int, NDim> primes_type;

	const primes_type primes;       // for hashing
	const int n_items;              // number of items
	const int n_entries;            // number of entries
	const int mask;                 // bitmask for valid range
	ndarray<1, key_type>    keys;   // voxel coordinates uniquely identifying a bucket
	ndarray<1, value_type>  values; // bucket description, or where to look in pivot array

public:
	// construct by zipping keys and values range
	template<typename items_range>
	HashMap(const items_range& items) :
		primes(init_primes()),
		n_items(boost::distance(items)),
		n_entries(init_entries()),
		mask(n_entries - 1),
		keys(init_keys()),
		values(init_values())
	{
		//mark grid as unoccupied
		for (auto item : items)
			write(boost::get<0>(item), boost::get<1>(item));
	}

	inline const value_type operator[](const key_type& key) const
	{
		int entry = get_hash(key);			//hash guess
		while (true)						//find the right entry
		{
			if ((keys[entry] == key).all())
				return values[entry];	                // we found it; this should be the most common code path
			if (values[entry] == -1) return -1;	    	// if we didnt find it yet by now, we never will
			entry = (entry + 1) & mask;		            // circular increment
		}
	}

private:
	// copy required number of primes into constant array
	primes_type init_primes() const {
		primes_type primes;
		for (auto i : boost::irange(0, NDim))
			primes(i) = PRIMES[i];
		return primes;
	}

	inline void write(const key_type& key, const value_type value)
	{
		int entry = get_hash(key);				// get entry initial guess
		while (true)                            // find an empty entry
		{
			if (values[entry] == -1) break;     // found an empty entry
			entry = (entry + 1) & mask;	        // circular increment
		}
		values[entry] = value;
		keys[entry] = key;
	}

	inline int get_hash(const key_type& key) const {
		return (key.cast<int>() * primes).redux(std::bit_xor<int>()) & mask;
	}

	int init_entries() const
	{
		//calc number of entries in hashmap. hashmap should have twice the number of items, at mimimum.
		int entries = 64;
		while (entries < n_items * 2) entries <<= 1;
		return entries;
	}

	ndarray<1, key_type> init_keys() const {
		return ndarray<2, key_type_scalar>({n_entries, NDim}).view<key_type>();
	}

	ndarray<1, value_type> init_values() const {
		ndarray<1, value_type> values({n_entries});
		fill(values, -1);
		return values;
	}

};

//// for checking if it matters anything in terms of performance
//class StdMap : public std::map
//{
//
//};
