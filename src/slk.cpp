#include "common.h"
#include "slk.h"

//' initial contiguity and distance matrices. The contiguity matrix is between
//' clusters, so is constantly modified, whereas the distance matrix is between
//' edges, so is fixed at load time.
void mats_init (const Rcpp::DataFrame &gr,
        arma::Mat <unsigned short> &contig_mat,
        arma::Mat <double> &d_mat)
{
    Rcpp::IntegerVector from = gr ["from"];
    Rcpp::IntegerVector to = gr ["to"];
    Rcpp::NumericVector d = gr ["d"];

    const int n = to.length ();

    contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    d_mat = arma::zeros <arma::Mat <double> > (n, n);

    for (int i = 0; i < n; i++)
    {
        contig_mat [from [i], to [i]] = 1;
        contig_mat [to [i], from [i]] = 1;
        d_mat [from [i], to [i]] = d [i];
        d_mat [to [i], from [i]] = d [i];
    }
}

//' merge two clusters in the contiguity matrix, reducing the size of the matrix
//' by one row and column.
bool merge_clusters (arma::Mat <unsigned short> &contig_mat,
        uint_map_t &edge2cl_map, uint_set_map_t &cl2edge_map,
        int i, int merge_from, int merge_to)
{
    bool merged = false;
    if (contig_mat [i, merge_from] == 1 && // TODO: <- or!
            contig_mat [merge_from, i] == 1 &&
            contig_mat [merge_to, i] == 1 &&
            contig_mat [merge_from, i] == 1)
    {
        contig_mat [i, merge_from] = 1;
        contig_mat [merge_from, i] = 1;
        contig_mat [merge_to, i] = 1;
        contig_mat [merge_from, i] = 1;
        contig_mat.shed_row (merge_from);
        contig_mat.shed_col (merge_from);

        std::set <unsigned int> edges_i = cl2edge_map.at (merge_from);
        cl2edge_map.erase (merge_from);
        std::set <unsigned int> edges_j = cl2edge_map.at (merge_to);
        for (auto e: edges_i)
            edges_j.insert (e);
        cl2edge_map [merge_to] = edges_j;

        merged = true;
    }

    return merged;
}

//' does the edge ei from graph_full connect two contiguous clusters?
bool does_edge_connect (arma::Mat <unsigned short> &contig_mat,
        uint_map_t edge2cl_map,
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        int ei)
{
    //int cl_fr = edge2cl_map.at [from [ei]],
    //    cl_to = edge2cl_map.at [to [ei]];

    //return (contig_mat [cl_fr, cl_to] == 1);
    return false;
}

int find_shortest_connection (
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        Rcpp::NumericVector &d,
        uint_map_t edge2cl_map,
        uint_set_map_t cl2edge_map,
        int merge_from,
        int merge_to)
{
    std::set <unsigned int> edges_i = cl2edge_map.at (merge_from);
    std::set <unsigned int> edges_j = cl2edge_map.at (merge_to);

    double dmin = INFINITE_DOUBLE;
    int shortest = INFINITE_INT;

    for (int i = 0; i < from.length (); i++)
        if (d (from [i], to [i]) < dmin)
        {
            dmin = d (from [i], to [i]);
            shortest = i;
        }
       

    return shortest;
}


//' rcpp_slk
//'
//' Full-order single linkage cluster redcap algorithm
//'
//' @noRd
void rcpp_slk (const Rcpp::DataFrame &grfull,
        Rcpp::DataFrame &gr)
{
    Rcpp::IntegerVector from = gr ["from"];
    Rcpp::IntegerVector to = gr ["to"];
    Rcpp::NumericVector d = gr ["d"];

    arma::Mat <unsigned short> contig_mat;
    arma::Mat <double> d_mat;
    mats_init (gr, contig_mat, d_mat);
}
