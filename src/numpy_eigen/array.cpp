/*
defines vector types
*/

#pragma once

#include <Eigen/Dense>
#include <boost/range.hpp>


template <typename T, int C>
using erow = Eigen::Array<T, 1, C, Eigen::RowMajor>;

template <typename T, int R>
using ecol = Eigen::Array<T, R, 1, Eigen::RowMajor>;

template <typename T, int R, int C>
using earray = Eigen::Array<T, R, C, Eigen::RowMajor>;



// add range support to arrays: https://forum.kde.org/viewtopic.php?f=74&t=111602
namespace std {
	template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
	typename Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>::Scalar*
		begin(Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& v) {
		return v.data();
	}
	template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
	typename Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>::Scalar*
		end(Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& v) {
		return v.data() + v.size();
	}
} // std


namespace boost {
	template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
	inline BOOST_DEDUCED_TYPENAME range_difference<
		Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>>::type
		range_calculate_size(const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& rng)
	{
		return rng.size();
	}

	template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
	struct range_mutable_iterator< Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
	{
		typedef _Scalar* type;
	};

	template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
	struct range_const_iterator< Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
	{
		typedef _Scalar* type;
	};
} // boost
