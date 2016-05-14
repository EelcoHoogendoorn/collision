#pragma once

#include <limits>
#include <iostream>
#include <functional>
#include <algorithm>


#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>

#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/adjacent_filtered.hpp>
//#include <boost/range/adaptors.hpp>       // somehow gives a link error?

#include "linalg.cpp"
#include "ndarray.cpp"
#include "maps.cpp"


using namespace boost;
using namespace boost::adaptors;


/*
vertex grid spatial collision map
this provides O(1) spatial lookups with no requirements at all on the point layout, 
the datastructures are easy to construct and avoid dynamic allocation or cache-unfriendly structures

since we use this for surfaces embedded in 3d space, the virtual voxel grid is quite sparse
this makes a dense grid less attractive; plus the part of the datastructure that will require random access
fits snugly in L1 cache

the only way to make this faster would be to actively reorder the input points, to exploit temporal coherency in the lexsort
*/

template<typename cell_type_scalar, typename vector_type_scalar, int NDim>
class VertexGridHash {
    /*
    this datastructure allows for coarse/braod collision queries
    it is one of the most simple datastructures to implement and debug,
    and is chosen here with future GPU-porting in mind, as it maps well to parallel architectures
    */

public:
    typedef int32 index_type;       // 32 bit int is a fine size type
    typedef Eigen::Array<vector_type_scalar, 2, NDim> extents_type;
    typedef Eigen::Array<vector_type_scalar, 1, NDim> vector_type;
    typedef Eigen::Array<cell_type_scalar,   1, NDim> cell_type;
    typedef Eigen::Array<int,                1, NDim> strides_type;

	const ndarray<1, const vector_type>  position;    //  positions
	const index_type                n_points;
	const float                     lengthscale; // size of a virtual voxel

    const extents_type              extents;     // maximum extent of pointcloud; used to map coordinates to positive integers
	const cell_type                 size;        // number of virtual buckets in each direction; used to prevent out-of-bound lookup
	const strides_type              strides;     // for lex-ranking cells
public:
	ndarray<2, cell_type_scalar>   _cell_id; // the cell coordinates a vertex resides in
	const ndarray<1, const cell_type>   cell_id; // the cell coordinates a vertex resides in
	const ndarray<1,       index_type>  indices; // index array mapping the vertices to lexographically sorted order
	      ndarray<1,       index_type>  pivots;	 // boundaries between buckets in vertices as viewed under indices
	const index_type                    n_buckets;// number of cells
	const HashMap<cell_type, index_type, NDim> bucket_from_cell; // maps cell coordinates to bucket index

public:
	//interface methods
	ndarray<2, cell_type_scalar> get_cell_ids() {
		//std::cout << "converted array" << std::endl;
	 //   auto arr = ndarray<2, cell_type_scalar>((PyObject*)&(this->cell_id.array));		// crash at this line
	 //   std::cout << "converted array" << std::endl;
	    return this->_cell_id;
	}
	void set_cell_ids(ndarray<2, cell_type_scalar> cells) {int a=3;}
	ndarray<1, index_type> get_indices() const {return this->indices;}
	void set_indices(ndarray<1, index_type> indices) {int a=3;}
	ndarray<1, index_type> get_pivots() const {return this->pivots;}
	void set_pivots(ndarray<1, index_type> pivots){this->pivots = pivots;}

    // grid constructor
	explicit VertexGridHash(const ndarray<2, vector_type_scalar> position, const vector_type_scalar lengthscale):
		position    (position.view<const vector_type>()),
		n_points    (position.size()),
		lengthscale (lengthscale),
		extents     (init_extents()),
		size        (init_size()),
		strides     (init_strides()),
		_cell_id    ({ n_points, NDim }),
		cell_id     (init_cells()),
		indices     (init_indices()),
		pivots      ({n_points}),
		n_buckets   (init_pivots()),
		bucket_from_cell(
		    n_buckets,
		    boost::combine(       // create a map to invert the cell_from_bucket function
                irange(0, n_buckets) | transformed([&](auto b){return cell_from_bucket(b);}),
                irange(0, n_buckets)
		    )
		)
	{
	    // empty constructor; noice
	}

private:
	//determine extents of data
	extents_type init_extents() const
	{
	    extents_type ext;
		ext.row(0).fill(+std::numeric_limits<vector_type_scalar>::infinity());
		ext.row(1).fill(-std::numeric_limits<vector_type_scalar>::infinity());
		for (auto p: position)
		{
			ext.row(0) = ext.row(0).min(p);
			ext.row(1) = ext.row(1).max(p);
		}
		return ext;
	}
    // integer size of the domain
	cell_type init_size() const {
	    return transform(extents.row(1) - extents.row(0)).cast<cell_type_scalar>() + 1;	// use +0.5 before cast?
	}

	strides_type init_strides() const {
//		boost::partial_sum(size.cast<int>(), begin(strides), std::multiplies<int>());   // doesnt work somehow
        strides_type strides;
        strides(0) = 1;
		for (auto i : irange(1, NDim))
			strides(i) = strides(i-1) * size(i-1);
		return strides;
	}

	// determine grid cells
	ndarray<1, const cell_type> init_cells() const
	{
		// silly indirection, because we cannot yet allocate custom type
	    //_cell_id = ndarray<2, cell_type_scalar>({n_points, NDim});
	    ndarray<1, cell_type> cidx = _cell_id.view<cell_type>();
		for (auto v: irange(0, n_points))
			cidx[v] = cell_from_position(position[v]);
	    return _cell_id.view<const cell_type>();
	}
	// finds the index vector that puts the vertices in a lexographically sorted order
	ndarray<1, index_type> init_indices() const
	{
		// create index array, based on lexographical ordering
	    ndarray<1, index_type> idx({n_points});
	    // init with initial order; 0 to n
		boost::copy(irange(0, n_points), idx.begin());
        // branching-free lex ordering
		auto lex = [&](auto l, auto r) {return ((cell_id[l] - cell_id[r]).cast<int>() * strides).sum() > 0;};
        // do the sort; boost::integer_sort might be preferable?
		boost::sort(idx, lex);
		return idx;//.view<index_type>();
	}
	//divide the sorted vertices into buckets, containing vertices in the same virtual voxel
	index_type init_pivots()
	{
		index_type np = 0;		//number of pivots
		auto add_pivot = [&](auto b) {pivots[np] = b; np += 1;};

        auto res = indices
                    | transformed([&](auto i){return cell_id[i];})
                    | indexed(0)
                    | adjacent_filtered([](auto a, auto b){return (a.value() != b.value()).any();})
                    | transformed([](auto i){return i.index();});

        for (auto i : res)
			add_pivot(i);
		add_pivot(n_points);

		if (np >= n_points)
		    throw python_exception("every vertex is in its own cell; lengthscale probably needs to go way up");

//        int_1::extent_gen extents;
//        pivots.view(boost::extents[np - 1]);;;

		return np - 1;
	}


protected:
	//map a global coord into the grid local coords
	inline const vector_type transform(const vector_type& v) const {
		return (v - extents.row(0)) / lengthscale;
	}
	inline const cell_type cell_from_local_position(const vector_type& v) const {
		return (v).cast<cell_type_scalar>();	// defacto floor
	}
	inline const cell_type cell_from_position(const vector_type& v) const {
		return cell_from_local_position(transform(v));
	}
	//convert bucket index into cell coords
	inline const cell_type cell_from_bucket(index_type b) const {
		return cell_id[indices[pivots[b]]];
	}

protected:
	auto indices_from_bucket(index_type b) const {
		return (b == -1) ? irange(0, 0) : irange(pivots[b], pivots[b+1]);
    }
    auto indices_from_cell(const cell_type& cell) const {
        return indices_from_bucket(bucket_from_cell[cell]);
    }
	auto vertices_from_cell(const cell_type& cell) const {
		return indices_from_cell(cell)
		        | transformed([&](index_type i){return indices[i];});
    }

public:
    // public traversal interface; what this class is all about
	template <class F>
	void for_each_vertex_in_cell(const cell_type& cell, const F& body) const {
		for (const int i: indices_from_cell(cell))
			body(indices[i]);
	}

	//loop over each occupied cell in the grid
	template <class F>
	void for_each_cell(const F& body) const {
		for (index_type b: irange(0, n_buckets))
			body(cell_from_bucket(b));
	}

	//loop over bounding box and apply body
	template <class F>
	void for_each_vertex_in_bounding_box(const vector_type& gmin, const vector_type& gmax, const F& body) const
	{
		const auto in_box = [&](auto v)
		{
			const vector_type& vp = position[v];
			return !((vp < gmin).any() || (vp > gmax).any());
		};

		const vector_type lmin = transform(gmin);
		const vector_type lmax = transform(gmax);

		//intersected volume is not positive; bail
		if ((lmin.max(vector_type(0,0,0)) > lmax.min(size.cast<vector_type_scalar>())).any()) return;

		//compute local cell coords; constrain to [0-size)
		const cell_type lb =  cell_from_local_position(lmin)   .max(cell_type(0, 0, 0));
		const cell_type ub = (cell_from_local_position(lmax)+1).min(size              );

		//loop over all cells that intersect with bb
		for (auto x: irange(lb[0], ub[0]))
			for (auto y: irange(lb[1], ub[1]))
				for (auto z: irange(lb[2], ub[2]))
					for (auto v: vertices_from_cell(cell_type(x, y, z)))
						if (in_box(v))
						    body(v);
	}
	//for unit testing purposes
	template <class F>
	void for_each_vertex_in_bounding_box_naive(const vector_type& gmin, const vector_type& gmax, const F& body) const
	{
		const auto in_box = [&](auto v)
		{
			const vector_type& vp = position[v];
			return !((vp < gmin).any() || (vp > gmax).any());
		};

		if ((gmin > extents.row(1)).any() || (gmax < extents.row(0)).any()) return;

		for (auto v: irange(0, n_points))
			if (in_box(v))
			    body(v);
	}

};