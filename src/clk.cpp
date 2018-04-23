#include "common.h"
#include "utils.h"
#include "clk.h"

// --------- COMPLETE LINKAGE CLUSTER ----------------


//' rcpp_clk
//'
//' Full-order complete linkage cluster redcap algorithm
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_clk (
        const Rcpp::DataFrame gr_full,
        const Rcpp::DataFrame gr)
{
    Rcpp::IntegerVector from_full_ref = gr_full ["from"];
    Rcpp::IntegerVector to_full_ref = gr_full ["to"];
    Rcpp::NumericVector d_full_ref = gr_full ["d"];
    Rcpp::IntegerVector from_ref = gr ["from"];
    Rcpp::IntegerVector to_ref = gr ["to"];
    Rcpp::NumericVector d_ref = gr ["d"];

    // Rcpp classes are always passed by reference, so cloning is necessary to
    // avoid modifying the original data.frames.
    Rcpp::IntegerVector from_full = Rcpp::clone (from_full_ref);
    Rcpp::IntegerVector to_full = Rcpp::clone (to_full_ref);
    Rcpp::NumericVector d_full = Rcpp::clone (d_full_ref);
    Rcpp::IntegerVector from = Rcpp::clone (from_ref);
    Rcpp::IntegerVector to = Rcpp::clone (to_ref);
    Rcpp::NumericVector d = Rcpp::clone (d_ref);

    // Index vectors are 1-indexed, so
    from_full = from_full - 1;
    to_full = to_full - 1;
    from = from - 1;
    to = to - 1;

    arma::Mat <unsigned short> contig_mat;
    arma::Mat <double> d_mat, d_mat_full;

    std::vector <int> treevec;

    return Rcpp::wrap (treevec);
}
