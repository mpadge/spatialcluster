// --------- AVERAGE LINKAGE CLUSTER ----------------

struct Edge_tree
{
    unsigned int n;
	Tree <double> * tree;
    dint_map_t edgewt2id_map;
    intd_map_t id2edgewt_map;

    arma::Mat <unsigned short> contig_mat, num_edges;
    arma::Mat <double> dmat, avg_dist;

    uint_map_t vert2cl_map;
    uint_set_map_t cl2vert_map;
};

void edge_tree_init (Edge_tree * edge_tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

void edge_tree_step (Edge_tree * edge_tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d,
        std::unordered_set <unsigned int> &the_tree);

Rcpp::IntegerVector rcpp_alk (
        const Rcpp::DataFrame gr);
