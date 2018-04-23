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
        const Rcpp::DataFrame gr)
{
    Rcpp::IntegerVector from_ref = gr ["from"];
    Rcpp::IntegerVector to_ref = gr ["to"];
    Rcpp::NumericVector d = gr ["d"];
    // Rcpp classes are always passed by reference, so cloning is necessary to
    // avoid modifying the original data.frames.
    Rcpp::IntegerVector from = Rcpp::clone (from_ref);
    Rcpp::IntegerVector to = Rcpp::clone (to_ref);
    // Index vectors are 1-indexed, so
    from = from - 1;
    to = to - 1;

    CLDAT clk_dat;

    std::vector <int> treevec;

    return Rcpp::wrap (treevec);
}
