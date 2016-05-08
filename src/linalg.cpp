/*
defines vector types
*/

#pragma once
#include "typedefs.cpp"
#include "Eigen/Dense"


typedef Eigen::Matrix<float,3,1> _float3;	//workhorse
typedef Eigen::Matrix<int,3,1>   _int3;		//used as intvec and index tuple
typedef Eigen::Matrix<float,3,3> _float33;	//used for triangle geometry computations

typedef Eigen::Array<float,3,1> float3;		//workhorse
typedef Eigen::Array<int,3,1>   int3;		//used as intvec and index tuple
typedef Eigen::Array<float,3,3> float33;	//used for triangle geometry computations


//add iteration support for int vector types
/*
namespace std
{
	template <>
	const int* begin(const int3& v)
	{
		return &(v[0]);
	}
	template <>
	const int* end(const int3& v)
	{
		return &(v[3]);
	}
	template <>
	int* begin(int3& v)
	{
		return &(v[0]);
	}
	template <>
	int* end(int3& v)
	{
		return &(v[3]);
	}

}*/

/*
class X
{
    // ...
	friend const int* begin(const int3& v)
	{
		return &(v[0]);
	}
	friend const int* end(const int3& v)
	{
		return &(v[3]);
	}
	friend int* begin(int3& v)
	{
		return &(v[0]);
	}
	friend int* end(int3& v)
	{
		return &(v[3]);
	}
};
*/