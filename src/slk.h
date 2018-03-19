#pragma once

typedef arma::Mat <unsigned short> sint_mat_t;

typedef std::unordered_map <unsigned int, unsigned int> uint_map_t;
typedef std::unordered_map <unsigned int,
        std::set <unsigned int> > uint_set_map_t;

unsigned int get_n (
    const Rcpp::IntegerVector &from,
    const Rcpp::IntegerVector &to);

void mats_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        arma::Mat <unsigned short> &contig_mat,
        arma::Mat <double> &d_mat,
        const unsigned int n);

void dmat_full_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        arma::Mat <double> &d_mat,
        const unsigned int n);

void sets_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        uint_map_t &edge2cl_map,
        uint_set_map_t &cl2edge_map);

int find_shortest_connection (
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        Rcpp::NumericVector &d,
        arma::Mat <double> &d_mat,
        uint_set_map_t &cl2edge_map,
        int cfrom,
        int cto);

void merge_clusters (
        arma::Mat <unsigned short> &contig_mat,
        uint_map_t &edge2cl_map,
        uint_map_t &cl2edge_map,
        int i,
        int merge_from,
        int merge_to);

Rcpp::IntegerVector rcpp_slk (
        const Rcpp::DataFrame &gr_full,
        Rcpp::DataFrame &gr);
