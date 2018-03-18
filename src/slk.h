#pragma once

typedef arma::Mat <unsigned short> sint_mat_t;

typedef std::unordered_map <unsigned int, unsigned int> uint_map_t;;
typedef std::unordered_map <unsigned int,
        std::set <unsigned int> > uint_set_map_t;;

void mats_init (const Rcpp::DataFrame &gr,
        arma::Mat <unsigned short> &contig_mat,
        arma::Mat <double> &d_mat);

bool merge_clusters (arma::Mat <unsigned short> &contig_mat,
        uint_map_t &edge2cl_map, uint_map_t &cl2edge_map,
        int i, int merge_from, int merge_to);

bool does_edge_connect (arma::Mat <unsigned short> &contig_mat,
        uint_map_t edge2cl_map,
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        int ei);

int find_shortest_connection (
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        Rcpp::NumericVector &d,
        uint_map_t edge2cl_map,
        uint_set_map_t cl2edge_map,
        int merge_from,
        int merge_to);

void rcpp_slk (const Rcpp::DataFrame &grfull,
        Rcpp::DataFrame &gr);
