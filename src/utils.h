#pragma once

namespace utils {

bool strfound (const std::string str, const std::string target);

struct OneEdge
{
    int from, to;
    double dist;
};

size_t sets_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        int2indx_map_t &vert2index_map,
        indx2int_map_t &index2vert_map,
        indx2int_map_t &index2cl_map,
        int2indxset_map_t &cl2index_map);

// TODO: Next 2 are only used in slk, so move em there
void mats_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const int2indx_map_t &vert2index_map,
        arma::Mat <int> &contig_mat,
        arma::Mat <double> &d_mat,
        bool shortest);

void dmat_full_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const int2indx_map_t &vert2index_map,
        arma::Mat <double> &d_mat,
        bool shortest);

size_t find_shortest_connection (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const int2indx_map_t &vert2index_map,
        const arma::Mat <double> &d_mat,
        const int2indxset_map_t &cl2index_map,
        const int cfrom,
        const int cto,
        const bool shortest);

void merge_clusters (
        arma::Mat <int> &contig_mat,
        indx2int_map_t &index2cl_map,
        int2indxset_map_t &cl2index_map,
        int merge_from,
        int merge_to);

void reconnect_cluster (
        arma::Mat <int> &contig_mat,
        const arma::Mat <double> &d_mat,
        const indx2int_map_t &index2cl_map,
        const int2indxset_map_t &cl2index_map,
        const int clnum);

} // end namespace utils
