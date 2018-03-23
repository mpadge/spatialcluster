#include "common.h"
#include "utils.h"
#include "bst.h"
#include "alk.h"

// --------- AVERAGE LINKAGE CLUSTER ----------------

void edge_tree_init (Edge_tree * edge_tree,
        const Rcpp::DataFrame gr)
{
    Rcpp::IntegerVector from = gr ["from"];
    Rcpp::IntegerVector to = gr ["to"];
    Rcpp::NumericVector d = gr ["d"];

    for (int i = 0; i < from.size (); i++)
    {
        if (i == 0)
            edge_tree->tree = treeNewNode (d [0]);
        else
            treeInsertNode (edge_tree->tree, d [i]);
        edge_tree->edgewt2id_map.emplace (d [i], i);
        edge_tree->id2edgewt_map.emplace (i, d [i]);
    }
}

//' rcpp_alk
//'
//' Full-order average linkage cluster redcap algorithm
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_alk (
        const Rcpp::DataFrame gr)
{
    Edge_tree edge_tree;
    edge_tree_init (&edge_tree, gr);

    Rcpp::IntegerVector res;
    return res;
}
