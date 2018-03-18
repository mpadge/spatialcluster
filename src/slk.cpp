#include "common.h"
#include "slk.h"

unsigned int get_n (
    Rcpp::IntegerVector &from,
    Rcpp::IntegerVector &to)
{
    std::set <unsigned int> nodes;
    for (auto i: from)
        nodes.insert (i);
    for (auto i: to)
        nodes.insert (i);
    return nodes.size ();
}

//' initial contiguity and distance matrices. The contiguity matrix is between
//' clusters, so is constantly modified, whereas the distance matrix is between
//' edges, so is fixed at load time.
void mats_init (
        const Rcpp::DataFrame &gr,
        arma::Mat <unsigned short> &contig_mat,
        arma::Mat <double> &d_mat,
        const unsigned int n)
{
    Rcpp::IntegerVector from = gr ["from"];
    Rcpp::IntegerVector to = gr ["to"];
    Rcpp::NumericVector d = gr ["d"];

    contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    d_mat = arma::zeros <arma::Mat <double> > (n, n);

    for (int i = 0; i < from.length (); i++)
    {
        contig_mat [from [i], to [i]] = 1;
        contig_mat [to [i], from [i]] = 1;
        d_mat [from [i], to [i]] = d [i];
        d_mat [to [i], from [i]] = d [i];
    }
}

void sets_init (
        const Rcpp::DataFrame &gr,
        uint_map_t &edge2cl_map,
        uint_set_map_t &cl2edge_map)
{
    Rcpp::IntegerVector from = gr ["from"];
    Rcpp::IntegerVector to = gr ["to"];

    edge2cl_map.clear ();
    cl2edge_map.clear ();

    unsigned int clnum = 0;

    for (int i = 0; i < from.length (); i++)
    {
        std::set <unsigned int> eset;
        eset.insert (from [i]);
        cl2edge_map.emplace (from [i], eset);
        edge2cl_map.emplace (i, clnum++);
    }
}

//' merge two clusters in the contiguity matrix, reducing the size of the matrix
//' by one row and column.
bool merge_clusters (
        arma::Mat <unsigned short> &contig_mat,
        uint_map_t &edge2cl_map,
        uint_set_map_t &cl2edge_map,
        int i,
        int merge_from,
        int merge_to)
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
bool does_edge_connect (
        arma::Mat <unsigned short> &contig_mat,
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
void rcpp_slk (
        const Rcpp::DataFrame &grfull,
        Rcpp::DataFrame &gr)
{
    Rcpp::IntegerVector from = gr ["from"];
    Rcpp::IntegerVector to = gr ["to"];
    Rcpp::NumericVector d = gr ["d"];

    const unsigned int n = get_n (from, to);

    arma::Mat <unsigned short> contig_mat;
    arma::Mat <double> d_mat;
    uint_map_t edge2cl_map;
    uint_set_map_t cl2edge_map;

    mats_init (gr, contig_mat, d_mat, n);
    sets_init (gr, edge2cl_map, cl2edge_map);

    std::set <unsigned int> the_tree;
    int i = 0;
    while (the_tree.size () < n)
    {
        unsigned int efrom = edge2cl_map.at (i),
                     eto = edge2cl_map.at (i);
        if (contig_mat [efrom, eto] > 0)
        {
            unsigned int ishort = find_shortest_connection (from, to, d,
                    edge2cl_map, cl2edge_map, efrom, eto);
            the_tree.insert (ishort);
            unsigned int cli = edge2cl_map.at (efrom),
                         clj = edge2cl_map.at (eto);
            merge_clusters (contig_mat, edge2cl_map, cl2edge_map,
                    i, cli, clj);
            i = 0;
        } else
        {
            i++;
        }
    }
}
