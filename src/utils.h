#pragma once

unsigned int sets_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        uint_map_t &vert2index_map,
        uint_map_t &index2vert_map,
        uint_map_t &index2cl_map,
        uint_set_map_t &cl2index_map);

void mats_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const uint_map_t &vert2index_map,
        arma::Mat <unsigned short> &contig_mat,
        arma::Mat <double> &d_mat);

void dmat_full_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const uint_map_t &vert2index_map,
        arma::Mat <double> &d_mat);

int find_shortest_connection (
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        Rcpp::NumericVector &d,
        uint_map_t &vert2index_map,
        arma::Mat <double> &d_mat,
        uint_set_map_t &cl2index_map,
        const unsigned int cfrom,
        const unsigned int cto);

void merge_clusters (
        arma::Mat <unsigned short> &contig_mat,
        uint_map_t &index2cl_map,
        uint_set_map_t &cl2index_map,
        const unsigned int merge_from,
        const unsigned int merge_to);

struct oneEdge
{
    unsigned int from, to;
    double dist;
};

bool edge_sorter (oneEdge const & lhs, oneEdge const & rhs);
