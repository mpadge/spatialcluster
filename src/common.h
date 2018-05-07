#pragma once

#include <limits>
#include <map>
#include <algorithm> // std::find
#include <vector>
#include <limits>
#include <random>
#include <string> // stoi
#include <cmath> // round
#include <unordered_set>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


constexpr float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
constexpr double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
constexpr int INFINITE_INT =  std::numeric_limits<int>::max ();

typedef int node_id_t;
typedef int cluster_id_t;
typedef size_t index_t;

typedef std::unordered_map <index_t, int> indx2int_map_t;
typedef std::unordered_map <int, index_t> int2indx_map_t;
typedef std::unordered_set <int> intset_t;
typedef std::unordered_map <int, intset_t> int2intset_map_t;

typedef std::unordered_map <unsigned int, unsigned int> uint_map_t;
typedef std::unordered_map <unsigned int,
        std::unordered_set <unsigned int> > uint_set_map_t;
