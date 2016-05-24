template<typename spec_t>
class PointGrid {
	/*
	provide spatial lookup in O(1) time for n-dimensional point clouds
	*/

public:
    typedef PointGrid<spec_t>				self_t;
    typedef typename spec_t::real_t         real_t;
	typedef typename spec_t::index_t		index_t;
	typedef typename spec_t::hash_t			hash_t;
	typedef typename spec_t::fixed_t		fixed_t;

	typedef typename spec_t::box_t			box_t;
	typedef typename spec_t::vector_t		vector_t;
	typedef typename spec_t::cell_t			cell_t;

	const spec_t				 spec;
	const ndarray<vector_t>      position;    // positions
	const index_t                n_points;    // number of points

public:
	const ndarray<fixed_t>       cell_id;     // the cell coordinates a vertex resides in
	const ndarray<index_t>       offsets;	  // determines stencil

public:
	//interface methods
	auto get_cells()        const { return cell_id; }
	void set_cells          (ndarray<fixed_t> cells)        {}


	// constructor
	explicit PointGrid(spec_t spec, ndarray<real_t, 2> position) :
		spec(spec),
		position(position.view<vector_t>()),
		n_points(position.size()),
		cell_id(init_cells()),
	{
	}
	// constructor
	explicit PointGrid(spec_t spec, ndarray<real_t, 2> position, ndarray<index_t> offsets) :
		spec(spec),
		position(position.view<vector_t>()),
		n_points(position.size()),
		cell_id(init_cells()),
	{
	}

	// create a new pointgrid, using the permutation of existing pointgrid as initial guess
	self_t update(const ndarray<real_t, 2> position) {
	    return self_t(spec, position, permutation, offsets);
	}
	explicit PointGrid(spec_t spec, ndarray<real_t, 2> position, ndarray<index_t> permutation, ndarray<index_t> offsets) :
		spec		(spec),
		position	(position.view<vector_t>()),
		n_points	(position.size()),
		cell_id		(init_cells()),
	{
	}

private:
	// determine grid cells
	auto init_cells() const {
		auto cell_id = ndarray<fixed_t>({ n_points });
		for (index_t v : irange(0, n_points))
			cell_id[v] = spec.hash_from_cell(spec.cell_from_position(position[v]));
		return cell_id;
	}


protected:
	//convert bucket index into cell coords
	inline fixed_t cell_from_bucket(index_t b) const {
		return cell_id[permutation[pivots[b]]];
	}
	auto indices_from_bucket(index_t b) const {
		return (b == -1) ? irange(0, 0) : irange(pivots[b], pivots[b + 1]);
	}
	auto indices_from_cell(fixed_t cell) const {
		return indices_from_bucket(bucket_from_cell[cell]);
	}
	auto vertices_from_cell(fixed_t cell) const {
		return indices_from_cell(cell)
			| transformed([&](index_t i) {return permutation[i];});
	}

public:
	// public traversal interface; what this class is all about
	template <class F>
	void for_each_vertex_in_cell(fixed_t cell, const F& body) const {
		for (index_t v : vertices_from_cell(cell))
			body(v);
	}

	//loop over each occupied cell in the grid
	template <class F>
	void for_each_cell(const F& body) const {
		for (index_t b : irange(0, n_buckets))
			body(cell_from_bucket(b));
	}


	//symmetric iteration over all vertex pairs, summing reaction forces
	template <class F>
	void for_each_vertex_pair(const F& body) const
	{
		const real_t ls2 = spec.scale*spec.scale;

		//this function wraps sign conventions regarding relative position and force
		const auto wrapper = [&](const index_t vi, const index_t vj){
			const vector_t rp = position[vi] - position[vj];
			const real_t d2 = (rp*rp).sum();
			if (d2 > ls2) return;
			body(vi, vj, d2);  // store index pair
		};

		//loop over all buckets
		for_each_cell([&](const fixed_t ci) {
            const auto bi = vertices_from_cell(ci);
			//interaction within bucket
			for (index_t vi : bi)
				for (index_t vj : bi)
					if (vi==vj) break; else
						wrapper(vi, vj);
			//loop over all neighboring buckets
			for (index_t o : offsets) {
				const auto bj = vertices_from_cell(ci + o);
				for (const index_t vj : bj)		//loop over other guy first; he might be empty, giving early exit
					for (const index_t vi : bi)
						wrapper(vi, vj);
			}
		});
	}
    // compute [n, 2] array of all paris within length_scale distance
	ndarray<index_t, 2> get_pairs() const
	{
	    typedef Eigen::Array<index_t, 1, 2> pair_t;

	    std::vector<pair_t> pairs;
	    for_each_vertex_pair([&](index_t i, index_t j, real_t d2) {
	        pairs.push_back(pair_t(i, j));
	    });
	    index_t n_pairs(pairs.size());
//	    ndarray<pair_t> _pairs({n_pairs});
//        boost::copy(pairs, _pairs.begin());
//        ndarray<index_t, 2> output = _pairs.unview<index_t>();
        ndarray<index_t, 2> output({ n_pairs, 2});
        for (index_t i : irange(0, n_pairs))
        {
            output[i][0] = pairs[i][0];
            output[i][1] = pairs[i][1];
        }
        return output;
	}
};


