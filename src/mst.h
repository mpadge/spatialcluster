#pragma once

#include <algorithm>
#include <Rcpp.h>

struct MSTEdge {
    int from, to;
    double dist;
    bool operator<(const MSTEdge& rhs) const { dist < rhs.dist; }
};

std::vector <MSTEdge> mst (Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);
