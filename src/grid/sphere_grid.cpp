
#pragma once

#include "stdafx.h"

#include "grid_spec.cpp"
#include "object_grid.cpp"


template<typename spec_t>
class SphereGrid : public ObjectGrid<spec_t> {

    /*
    closest analogue to pointgrid, but with higher construction cost
    allows for variable radius though
    place a marker on any cell that intersects with the sphere
    */

};