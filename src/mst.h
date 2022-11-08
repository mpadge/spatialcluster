#pragma once

#include <Rcpp.h>

struct MSTEdge {
    int from, to;
    double dist;
    bool operator < (MSTEdge const& other) {
        return dist < other.dist;
    };
};

std::vector <MSTEdge> mst (Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);
