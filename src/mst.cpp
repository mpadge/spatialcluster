#include "mst.h"

std::vector <MSTEdge> mst (Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d) {
    const size_t n = static_cast <size_t> (from.size ());

    std::vector <MSTEdge> edges (n);
    for (size_t i = 0; i < n; i++) {
        MSTEdge ei;
        ei.from = from (i);
        ei.to = to (i);
        ei.dist = d (i);
        edges [i] = ei;
    }

    std::vector <size_t> cl_id (n);
    for (size_t i = 0; i < n; i++) {
        cl_id [i] = i;
    }

    std::sort (edges.begin (), edges.end ());

    std::vector <MSTEdge> result;

    for (MSTEdge e : edges) {
        const size_t cl_from = cl_id [static_cast <size_t> (e.from)],
            cl_to = cl_id [static_cast <size_t> (e.to)];

        if (cl_from != cl_to) {
            result.push_back (e);

            const size_t cl_min = std::min (cl_from, cl_to),
                  cl_max = std::max (cl_from, cl_to);

            for (size_t i = 0; i < n; i++) {
                if (cl_id [i] == cl_max) {
                    cl_id [i] = cl_min;
                }
            }
        }
    }

    return result;
}


//' rcpp_mst
//'
//' Minimum spanning tree
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::DataFrame rcpp_mst (Rcpp::DataFrame input) {
    Rcpp::IntegerVector from = input ["from"];
    Rcpp::IntegerVector to = input ["to"];
    Rcpp::NumericVector d = input ["d"];

    std::vector <MSTEdge> tree = mst (from, to, d);

    Rcpp::IntegerVector from_out (tree.size ());
    Rcpp::IntegerVector to_out (tree.size ());
    Rcpp::NumericVector d_out (tree.size ());
    size_t i = 0;
    for (auto t: tree) {
        from_out (i) = t.from;
        to_out (i) = t.to;
        d_out (i) = t.dist;
        i++;
    }

    Rcpp::DataFrame res = Rcpp::DataFrame::create (
        Rcpp::Named ("from") = from_out,
        Rcpp::Named ("to") = to_out,
        Rcpp::Named ("d") = d_out,
        Rcpp::_["stringsAsFactors"] = false);

    return res;
};
