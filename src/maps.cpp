
#pragma once

#include <functional>
#include <array>

#include <boost/range/irange.hpp>
#include <boost/range/difference_type.hpp>
#include <boost/tuple/tuple.hpp>

#include "typedefs.cpp"
#include "numpy_eigen/array.cpp"
#include "numpy_boost/ndarray.cpp"

//used for hashing calcs
const std::array<int, 3> PRIMES = { 73856093, 19349663, 83492791 };


template<typename key_type, typename value_type, int NDim>
class HashMap {

	typedef int64                                       primes_type_scalar;
	typedef Eigen::Array<primes_type_scalar, 1, NDim>   primes_type;
	typename typedef key_type::Scalar                   key_type_scalar;

	const primes_type   primes;    // for hashing
	const int           n_items;   // number of items
	const int           n_entries; // number of entries
	const int           mask;      // bitmask for valid range
	ndarray<key_type>   keys;      // voxel coordinates uniquely identifying a bucket
	ndarray<value_type> values;    // bucket description, or where to look in pivot array

public:
	// construct by zipping keys and values range
	template<typename items_range>
	HashMap(const items_range& items) :
		primes      (init_primes()),
		n_items     (boost::distance(items)),
		n_entries   (init_entries()),
		mask        (n_entries - 1),
		keys        (init_keys()),
		values      (init_values())
	{
		for (auto item : items)
			write(boost::get<0>(item), boost::get<1>(item));
	}

	inline value_type operator[](const key_type& key) const {
		int entry = get_hash(key);			//hash guess
		while (true) {						//find the right entry
			if ((keys[entry] == key).all())
				return values[entry];	                // we found it; this should be the most common code path
			if (values[entry] == -1) return -1;	    	// if we didnt find it yet by now, we never will
			entry = (entry + 1) & mask;		            // circular increment
		}
	}

private:
	// copy required number of primes into constant array
	auto init_primes() const {
		primes_type primes;
		for (auto i : boost::irange(0, NDim))
			primes(i) = PRIMES[i];
		return primes;
	}

	inline void write(const key_type& key, value_type value) {
		auto entry = get_hash(key);				// get entry initial guess
		while (true) {                          // find an empty entry
			if (values[entry] == -1) break;     // found an empty entry
			entry = (entry + 1) & mask;	        // circular increment
		}
		values[entry] = value;
		keys[entry] = key;
	}

	inline auto get_hash(const key_type& key) const {
		return (key.cast<primes_type_scalar>() * primes).redux(std::bit_xor<primes_type_scalar>()) & mask;
	}

	auto init_entries() const{
		//calc number of entries in hashmap. hashmap should have twice the number of items, at mimimum.
		int entries = 64;
		while (entries < n_items * 2) entries <<= 1;
		return entries;
	}

	auto init_keys() const {
		return ndarray<key_type_scalar, 2>({n_entries, NDim}).view<key_type>();
	}

	auto init_values() const {
		ndarray<value_type> values({n_entries});
		fill(values, -1); 		//mark grid as unoccupied
		return values;
	}

};

//// for checking if it matters anything in terms of performance
//class StdMap : public std::map
//{
//
//};
