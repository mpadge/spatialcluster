#pragma once

bool strfound (const std::string str, const std::string target);

int sets_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        int2indx_map_t &vert2index_map,
        indx2int_map_t &index2vert_map,
        indx2int_map_t &index2cl_map,
        int2indxset_map_t &cl2index_map);

void mats_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const int2indx_map_t &vert2index_map,
        arma::Mat <int> &contig_mat,
        arma::Mat <double> &d_mat);

void dmat_full_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const int2indx_map_t &vert2index_map,
        arma::Mat <double> &d_mat);

size_t find_shortest_connection (
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        Rcpp::NumericVector &d,
        int2indx_map_t &vert2index_map,
        arma::Mat <double> &d_mat,
        int2indxset_map_t &cl2index_map,
        const int cfrom,
        const int cto);

void merge_clusters (
        arma::Mat <int> &contig_mat,
        indx2int_map_t &index2cl_map,
        int2indxset_map_t &cl2index_map,
        int merge_from,
        int merge_to);

struct OneEdge
{
    int from, to;
    double dist;
};
