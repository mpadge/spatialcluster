#pragma once

#include <limits>
#include <map>
#include <algorithm> // std::find
#include <vector>
#include <limits>
#include <random>
#include <string> // stoi
#include <cmath> // round

#include <Rcpp.h>

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
const int INFINITE_INT =  std::numeric_limits<int>::max ();

typedef std::string vertex_id_t, edge_id_t;
typedef std::unordered_map <unsigned int,
    std::unordered_set <unsigned int> > int2ints_map_t;

