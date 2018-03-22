// --------- AVERAGE LINKAGE CLUSTER ----------------

struct Edge_tree
{
	Tree <double> *tree;
    dint_map_t edgewt2id_map;
    intd_map_t id2edgewt_map;
};

void edge_tree_init (Edge_tree *edge_tree,
        const Rcpp::DataFrame gr);

Rcpp::IntegerVector rcpp_alk (
        const Rcpp::DataFrame gr);
