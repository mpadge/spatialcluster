#pragma once

typedef arma::Mat <unsigned short> sint_mat_t;

void contig_mat_init (const Rcpp::DataFrame &gr,
        arma::Mat <unsigned short> &contig_mat);

void rcpp_slk (const Rcpp::DataFrame &grfull,
        Rcpp::DataFrame &gr);
