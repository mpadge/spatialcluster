#pragma once

#include <limits>
#include <map>
#include <algorithm> // std::find
#include <vector>
#include <limits>
#include <random>
#include <string> // stoi
#include <cmath> // round

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
const int INFINITE_INT =  std::numeric_limits<int>::max ();

typedef std::string vertex_id_t, edge_id_t;
typedef std::unordered_map <unsigned int,
    std::unordered_set <unsigned int> > int2ints_map_t;

typedef arma::Mat <unsigned short> sint_mat_t;

typedef std::unordered_map <unsigned int, unsigned int> uint_map_t;
typedef std::unordered_map <unsigned int,
        std::unordered_set <unsigned int> > uint_set_map_t;

typedef std::unordered_map <unsigned int, double> intd_map_t;
typedef std::unordered_map <double, unsigned int> dint_map_t;

/* The main matrices (contig, num_edges, dmat, avg_dist) are all referenced by
 * direct indices throughout, not by vertex numbers. The latter are mapped to
 * the former by vert2index_map.  The index2cl and cl2index then associate those
 * indices with clusters which are themselves also direct indices into the
 * matrices. Cluster merging simply re-directs multiple indices onto the same
 * cluster (index) numbers.
 */
struct CLDAT
{
    unsigned int n;

    std::unordered_map <double,
        std::pair <unsigned int, unsigned int> > edgewt2idx_pair_map;
    std::unordered_map <unsigned int, std::unordered_set <double> >
        idx2edgewt_map; // all wts associated with that cluster

    arma::Mat <unsigned short> contig_mat, num_edges;
    arma::Mat <double> dmat, avg_dist;

    uint_map_t index2cl_map, vert2index_map, index2vert_map;
    uint_set_map_t cl2index_map;
};

