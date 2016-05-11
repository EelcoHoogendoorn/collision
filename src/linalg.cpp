/*
defines vector types
*/

#pragma once
#include "typedefs.cpp"
#include <Eigen/Dense>


typedef Eigen::Matrix<float,3,1> _float3;	//workhorse
typedef Eigen::Matrix<int,3,1>   _int3;		//used as intvec and index tuple
typedef Eigen::Matrix<float,3,3> _float33;	//used for triangle geometry computations

typedef Eigen::Array<float,3,1> float3;		//workhorse
typedef Eigen::Array<int,3,1>   int3;		//used as intvec and index tuple
typedef Eigen::Array<float,3,3> float33;	//used for triangle geometry computations
