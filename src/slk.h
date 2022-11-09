#pragma once

// --------- SINGLE LINKAGE CLUSTER ----------------

Rcpp::IntegerVector rcpp_slk (
        const Rcpp::DataFrame gr_full,
        const Rcpp::DataFrame gr,
        const bool shortest,
        const bool quiet);
