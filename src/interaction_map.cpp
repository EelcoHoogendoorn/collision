#pragma once

#include <limits>

#include <boost/shared_ptr.hpp>

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/any_range.hpp>

#include <boost/iterator.hpp>
#include <boost/iterator/permutation_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "linalg.cpp"
#include "ndarray.cpp"



/*
vertex grid spatial collision map
this provides O(1) spatial lookups with no requirements at all on the point layout, 
the datastructures are easy to construct and avoid dynamic allocation or cache-unfriendly structures

since we use this for surfaces embedded in 3d space, the virtual voxel grid is quite sparse
this makes a dense grid less attractive; plus the part of the datastructure that will require random access
fits snugly in L1 cache

the only way to make this faster would be to actively reorder the input points, to exploit temporal coherency in the lexsort
*/


//used for hashing calcs
const int3 primes(73856093, 19349663, 83492791);

class HashMap {
	const int entries;	//number of entries in hashmap
	int_2 keys;	        //voxel coordinates uniquely identifying a bucket
	int_1 values;	    //bucket description, or where to look in pivot array

	inline int get_hash(const int3& p) const
	{
		const int3 c = p*primes;
		return (c[0]^c[1]^c[2]) & (entries-1);
	}

public:
	HashMap(const int entries):
	    entries(entries),
		keys({entries, 3}),
		values({entries})
	{
		//mark grid as unoccupied
		fill(values, -1);
	}

	inline void write(const int3& key, const int value)
	{
		const auto _keys = keys.range<int3>();
		int entry = get_hash(key);				//get entry initial guess
		// find an empty entry
		while (true)
		{
			if (values[entry] == -1) break;
			entry = (entry+1) & (entries-1);					//circular increment
		}
		values[entry] = value;
		_keys[entry] = key;
	}
	inline int read(const int3& key) const
	{
		const auto _keys  = keys.range<const int3>();
		int entry = get_hash(key);			//hash guess

		while (true)						//find the right entry
		{
			if ((_keys[entry] == key).all()) return values[entry];	//we found it; this should be the most common code path
			if (values[entry] == -1) return -1;	    				//if we didnt find it yet, we never will
			entry = (entry+1) & (entries-1);						//circular increment
		}
	}
};


class VertexGridHash {
    /*
    this datastructure allows for coarse/braod collision queries
    it is one of the most simple datastructures to implement and debug,
    and is chosen here with future GPU-porting in mind, as it maps well to parallel architectures
    */

public:
	float_2 position; //vertex positions
	const int vertices;     //number of vertices

	const float lengthscale;//size of a virtual voxel, which equals the maximum interaction range of a vertex-vertex pair query

private:
	float3 pmin, pmax;	    //maximum extent of pointcloud; used to map coordinates to positive integers
	int3 size;			    //number of virtual buckets in each direction; used to prevent out-of-bound lookup

	int_2 cell_id;	//the cell coordinates a vertex resides in
	int_1 indices;	//index array mapping the vertices to lexographically sorted order
	int_1 pivots;		//boundaries between buckets in vertices as viewed under indices
	int buckets;		    //number of buckets

	boost::shared_ptr<HashMap> hashmap;
	
public:
	//interface methods
	float_2 get_position() const {return this->position;}
	void set_position(float_2 position){this->position = position;}
	int_2 get_cell_id() const {return this->cell_id;}
	void set_cell_id(int_2 cell_id){this->cell_id = cell_id;}
	int_1 get_indices() const {return this->indices;}
	void set_indices(int_1 indices){this->indices = indices;}
	int_1 get_pivots() const {return this->pivots;}
	void set_pivots(int_1 pivots){this->pivots = pivots;}


	VertexGridHash(float_2 position, float lengthscale):
		position(position),
		vertices(position.size()),
		lengthscale(lengthscale),
		cell_id({vertices, 3}),
		indices({vertices}),
		pivots({vertices+1})
	{
		//build the datastructure
		sizing();
		indexing();
		pivoting();
		writegrid();
	}

private:
	//map a global coord into the grid local coords
	inline float3 transform(const float3& v) const
	{
		return (v-pmin) / lengthscale;
	}
	inline int3 cell_from_local_position(const float3& v) const
	{
		return (v-0.5).cast<int>();	//defacto floor
	}
	inline int3 cell_from_position(const float3& v) const
	{
		return cell_from_local_position(transform(v));
	}
	//convert bucket index into cell coords
	inline int3 cell_from_bucket(const int b) const
	{
		const auto _cell_id  = cell_id.range<const int3>();
		return _cell_id[indices[pivots[b]]];
	}


	//determine extents of data
	void sizing()
	{
		pmin.fill(+std::numeric_limits<float>::infinity()); 
		pmax.fill(-std::numeric_limits<float>::infinity());
		for (const float3 p: position.range<const float3>())
		{
			pmin = pmin.min(p);
			pmax = pmax.max(p);
		}

		size = transform(pmax).cast<int>()+1;
	}
	//finds the index vector that puts the vertices in a lexographically sorted order
	void indexing()
	{
		const auto _position = position.range<const float3>();
		const auto _cell_id  = cell_id .range<int3>();

		//determine grid cells
		for (const int v: boost::irange(0, vertices))
			_cell_id[v] = cell_from_position(_position[v]);

		//create index array, based on lexographical ordering
		boost::copy(boost::irange(0, vertices), indices.begin());
		boost::sort(
	        indices.range<int>(),
			[&](const int l, const int r) {
				const int3 cl = _cell_id[l];
				const int3 cr = _cell_id[r];
				return 
					cl[0]!=cr[0] ?
						cl[0]<cr[0]:
						cl[1]!=cr[1] ?
							cl[1]<cr[1]:
							cl[2]<cr[2];
        });
	}
	//divide the sorted vertices into buckets, containing vertices in the same virtual voxel
	void pivoting()
	{
		int np = 0;		//number of pivots
		const auto add_pivot = [&](const int b) {pivots[np] = b; np += 1;};
		const auto _cell_id  = cell_id .range<const int3>();

		add_pivot(0);
		for (const int i: boost::irange(1, vertices))
			if ((_cell_id[indices[i]] != _cell_id[indices[i-1]]).any())		//if different from the previous one in sorted space
				add_pivot(i);
		if (np==vertices) throw my_exception("every vertex is in its own cell; that cant be right, can it? lengthscale probably needs to go way up");
		add_pivot(vertices);
		buckets = np - 1;		//number of buckets is pivots-1
	}
	//build the hashmap to quickly access buckets of vertices at arbitrary locations
	void writegrid()
	{
		//calc number of entries in hashmap. hashmap should have twice the number of buckets, at mimimum.
		//absolute size is proportional to the number of vertices in the dataset, not to the space they occupy
		int entries = 64;
		while (entries < buckets*2) entries <<= 1;

        hashmap = boost::shared_ptr<HashMap>(new HashMap(entries));

		//write buckets to the hash map
		for (const int b: boost::irange(0, buckets))
			hashmap->write(cell_from_bucket(b), b);
	}


public:
    // public traversal interface; what this class is all about
	boost::integer_range<int> indices_from_bucket(const int b) const
	{
		if (b==-1){
			return boost::irange(0,0);			//create empty range
		}else{
			return boost::irange(pivots[b], pivots[b+1]);
		}
	}
	//abstract away the iteration over a gridcell
	auto vertices_from_cell(const int3& cell) const
	{
		const auto piv = indices_from_bucket(hashmap->read(cell));
		const auto ind = indices.range<int>();			
		return boost::make_iterator_range(
			boost::make_permutation_iterator( ind.begin(), piv.begin() ), 
			boost::make_permutation_iterator( ind.end(),   piv.end()   ));
	}
	template <class F>
	void for_each_vertex_in_cell(const int3& cell, const F& body) const
	{
		//would it help here if we permuted the bucket indices to achieve lexographic iteration?
		for (const int i: indices_from_bucket(hashmap->read(cell)))
			body(indices[i]);
	}

	//loop over each occupied cell in the grid
	template <class F>
	void for_each_cell(const F& body) const
	{
		//iterate over entries instead? only helps if entry stores full cell info
		for (const int b: boost::irange(0,buckets))
			body(cell_from_bucket(b));
	}
	//loop over bounding box and apply body
	template <class F>
	void for_each_vertex_in_bounding_box(const float3& gmin, const float3& gmax, const F& body) const
	{
		const auto _position = position.range<const float3>();
		const auto in_box = [&](const int v)
		{
			const float3 vp = _position[v];
			return !((vp<gmin).any() || (vp>gmax).any());
		};

		const float3 lmin = transform(gmin);
		const float3 lmax = transform(gmax);

		//intersected volume is not positive; bail
		if ((lmin.max(float3(0,0,0)) > lmax.min(size.cast<float>())).any()) return;				

		//compute local cell coords; constrain to [0-size)
		const int3 lb =  cell_from_local_position(lmin)   .max(int3(0,0,0));
		const int3 ub = (cell_from_local_position(lmax)+1).min(size       );

		//loop over all cells that intersect with bb
		for (const int x: boost::irange(lb[0],ub[0]))
			for (const int y: boost::irange(lb[1],ub[1]))
				for (const int z: boost::irange(lb[2],ub[2]))
					for (const int v: vertices_from_cell(int3(x,y,z)))
						if (in_box(v)) body(v);
	}
	//for unit testing purposes
	template <class F>
	void for_each_vertex_in_bounding_box_naive(const float3& gmin, const float3& gmax, const F& body) const
	{
		const auto _position = position.range<const float3>();
		const auto in_box = [&](const int v)
		{
			const float3 vp = _position[v];
			return !((vp<gmin).any() || (vp>gmax).any());
		};

		if ((gmin>pmax).any() || (gmax<pmin).any()) return;

		for (const int v: boost::irange(0, vertices))
			if (in_box(v)) body(v);
	}

};

