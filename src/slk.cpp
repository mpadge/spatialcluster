#include "common.h"
#include "slk.h"

void contig_mat_init (const Rcpp::DataFrame &gr,
        arma::Mat <unsigned short> &contig_mat)
{
    Rcpp::IntegerVector from = gr ["from"];
    Rcpp::IntegerVector to = gr ["to"];

    const int n = to.length ();

    contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);

    for (int i = 0; i < n; i++)
    {
        contig_mat [from [i], to [i]] = 1;
        contig_mat [to [i], from [i]] = 1;
    }
}

//' merge two clusters in the contiguity matrix, reducing the size of the matrix
//' by one row and column.
bool contig_mat_merge (arma::Mat <unsigned short> &contig_mat,
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

int find_shortest_connection (uint_map_t edge2cl_map,
        uint_map_t cl2edge_map)
{
    int shortest;

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
    contig_mat_init (gr, contig_mat);
}
