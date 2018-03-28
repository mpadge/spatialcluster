// --------- AVERAGE LINKAGE CLUSTER ----------------

#include <unordered_set>

struct Edge_tree
{
    unsigned int n;
	Tree <double> * tree;

    std::unordered_map <double,
        std::pair <unsigned int, unsigned int> > edgewt2clpair_map;

    arma::Mat <unsigned short> contig_mat, num_edges;
    arma::Mat <double> dmat, avg_dist;

    uint_map_t vert2cl_map, vert2index_map, index2vert_map;
    uint_set_map_t cl2vert_map;
};

void edge_tree_init (Edge_tree * edge_tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

int edge_tree_step (Edge_tree * edge_tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

Rcpp::IntegerVector rcpp_alk (
        const Rcpp::DataFrame gr);
