/**
 * @file VecNd.h
 * @brief Type definition for multi-D vectors
 *
 * Type definitions of multi dimensional vectors prdouce difficult to read
 * code. Instead, this VecNd serves a short hand for them.
 */

#pragma once

#include <vector>
namespace {

template<typename T, std::size_t dims>
struct MultiDVector{
	using type = std::vector<typename MultiDVector<T,dims-1>::type>;
};

template<typename T>
struct MultiDVector<T,0>{
	using type = T;
};

}

template<typename T,std::size_t dim>
using VectorNd = typename MultiDVector<T,dim>::type;
