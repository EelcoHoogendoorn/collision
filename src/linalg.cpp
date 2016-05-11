/*
defines vector types
*/

#pragma once
#include "typedefs.cpp"
#include <Eigen/Dense>


typedef Eigen::Matrix<float,1,3> _float3;	//workhorse
typedef Eigen::Matrix<int,1,3>   _int3;		//used as intvec and index tuple
typedef Eigen::Matrix<float,3,3> _float33;	//used for triangle geometry computations

typedef Eigen::Array<float,1,3> float3;		//workhorse
typedef Eigen::Array<int,1,3>   int3;		//used as intvec and index tuple
typedef Eigen::Array<float,3,3> float33;	//used for triangle geometry computations
typedef Eigen::Array<float,2,3> float23;	//min/max extents
//typedef Eigen::Array<float,3,2> float32;	//min/max extents