#pragma once

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
