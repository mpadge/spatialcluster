#pragma once

unsigned int sets_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        uint_map_t &vert2index_map,
        uint_map_t &index2vert_map,
        uint_map_t &edge2cl_map,
        uint_set_map_t &cl2edge_map);

void mats_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        uint_map_t &vert2index_map,
        arma::Mat <unsigned short> &contig_mat,
        arma::Mat <double> &d_mat);

void dmat_full_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        uint_map_t &vert2index_map,
        arma::Mat <double> &d_mat);

int find_shortest_connection (
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        Rcpp::NumericVector &d,
        uint_map_t &vert2index_map,
        arma::Mat <double> &d_mat,
        uint_set_map_t &cl2edge_map,
        int cfrom,
        int cto);

void merge_clusters (
        arma::Mat <unsigned short> &contig_mat,
        uint_map_t &vert2cl_map,
        uint_set_map_t &cl2vert_map,
        int merge_from,
        int merge_to);

