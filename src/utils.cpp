#include "common.h"
#include "utils.h"
#include <algorithm>

// Note that all matrices **CAN** be asymmetrical, and so are always indexed
// (from, to)

/* These routines all work with 2 lots of 2 main maps:
 * 1. vert2index and index2vert maps
 * 2. cl2index and index2cl maps
 *
 * The former map all vertices enumerated in the original from and to vectors to
 * sequential index numbers into the matrices (dists, contig_mat, whatever). The
 * latter are initially direct a->a maps of all indices to themselves. As
 * clusters merge, the values of index2cl maps are updated so that, for example,
 * index2cl(a)->b and index2cl(b)->b. The cl2index map then holds an
 * unordered_set of target indices for each cluster.
 */

bool utils::strfound (const std::string str, const std::string target)
{
    bool found = false;
    if (str.find (target) != std::string::npos)
        found = true;
    return found;
}

size_t utils::sets_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        int2indx_map_t &vert2index_map,
        indx2int_map_t &index2vert_map,
        indx2int_map_t &index2cl_map,
        int2indxset_map_t &cl2index_map)
{
    vert2index_map.clear ();
    index2vert_map.clear ();
    index2cl_map.clear ();
    cl2index_map.clear ();

    intset_t vert_set;
    for (int i = 0; i < from.size (); i++)
    {
        vert_set.emplace (from [i]);
        vert_set.emplace (to [i]);
    }
    int idx = 0;
    for (auto v: vert_set)
    {
        index2vert_map.emplace (idx, v);
        vert2index_map.emplace (v, idx++);
    }

    for (int i = 0; i < from.length (); i++)
    {
        size_t fi = vert2index_map.at (from [i]);
        indxset_t eset; // only one index for each vert at this stage
        eset.insert (fi);
        cl2index_map.emplace (fi, eset);
    }
    for (auto v: vert_set)
    {
        // all verts are their own clusters, so cast indx to cli
        int cli = static_cast <int> (vert2index_map.at (v));
        indxset_t eset;
        if (cl2index_map.find (cli) != cl2index_map.end ())
            eset = cl2index_map.at (cli);
        eset.emplace (cli);
        cl2index_map.emplace (cli, eset);
        index2cl_map.emplace (static_cast <size_t> (cli), cli);
    }
    
    return static_cast <size_t> (vert_set.size ());
}

//' find shortest (or longest) connection between two clusters
//' @param from, to, d the columns of the edge graph
//' @param d_mat distance matrix between all edges (not between clusters!)
//' @param cl2vert_map map of list of all (from, to, d) edges for each cluster
//' @param cfrom Number of cluster which is to be merged
//' @param cto Number of cluster with which it is to be merged
//'
//' @return Index directly into from, to - **NOT** into the actual matrices!
//' @noRd
size_t utils::find_shortest_connection (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const int2indx_map_t &vert2index_map,
        const arma::Mat <double> &d_mat,
        const int2indxset_map_t &cl2index_map,
        const int cfrom,
        const int cto,
        const bool shortest)
{
    if (cl2index_map.find (cfrom) == cl2index_map.end ())
        Rcpp::stop ("cluster index not found");
    if (cl2index_map.find (cto) == cl2index_map.end ())
        Rcpp::stop ("cluster index not found");

    indxset_t index_i = cl2index_map.at (cfrom),
             index_j = cl2index_map.at (cto);

    double dlim = INFINITE_DOUBLE;
    if (!shortest)
        dlim = -dlim;
    size_t short_i = INFINITE_INT, short_j = INFINITE_INT;

    // from and to here are directional, so need to examine both directions
    for (auto i: index_i)
        for (auto j: index_j)
        {
            arma::uword ia = static_cast <arma::uword> (i),
                        ja = static_cast <arma::uword> (j);
            if ((shortest && d_mat (ia, ja) < dlim) ||
                    (!shortest && d_mat (ia, ja) > dlim))
            {
                dlim = d_mat (ia, ja);
                short_i = i;
                short_j = j;
            } else if ((shortest && d_mat (ja, ia) < dlim) ||
                    (!shortest && d_mat (ja, ia) > dlim))
            {
                dlim = d_mat (ja, ia);
                short_i = j;
                short_j = i;
            }
        }
    if (dlim == INFINITE_DOUBLE)
        Rcpp::stop ("no minimal distance; this should not happen");

    // convert short_i and short_j to a single edge 
    // TODO: Make a std::map of vert2dist to avoid this loop
    size_t shortest_edge = INFINITE_INT;
    for (int i = 0; i < from.length (); i++) // int for Rcpp index
    {
        if ((vert2index_map.at (from [i]) == short_i &&
                vert2index_map.at (to [i]) == short_j) ||
            (vert2index_map.at (from [i]) == short_j &&
                vert2index_map.at (to [i]) == short_i))
        {
            shortest_edge = static_cast <size_t> (i);
            break;
        }
    }
    if (shortest_edge == INFINITE_INT)
        Rcpp::stop ("shite");

    return shortest_edge;
}

//' merge two clusters in the contiguity matrix, reducing the size of the matrix
//' by one row and column.
//'
//' @return A logical parameter indicating whether or not the newly formed
//' cluster has any outgoing connections.
//' @noRd
void utils::merge_clusters (
        arma::Mat <int> &contig_mat,
        indx2int_map_t &index2cl_map,
        int2indxset_map_t &cl2index_map,
        const int cluster_from,
        const int cluster_to)
{
    if (cluster_from < 0)
        Rcpp::stop ("cluster_from must be non-negative");
    if (cluster_to < 0)
        Rcpp::stop ("cluster_to must be non-negative");

    arma::uword cfr = static_cast <arma::uword> (cluster_from),
                cto = static_cast <arma::uword> (cluster_to);
    // Set all contig_mat (cluster_from, .) to 1
    for (arma::uword i = 0; i < contig_mat.n_rows; i++)
    {
        if (contig_mat (cfr, i) == 1 || contig_mat (i, cfr) == 1)
        {
            contig_mat (cfr, i) = contig_mat (i, cfr) = 1;
            contig_mat (cto, i) = contig_mat (i, cto) = 1;
        }
    }

    indxset_t idx_from = cl2index_map.at (cluster_from),
              idx_to = cl2index_map.at (cluster_to);

    for (auto i: idx_from)
        for (auto j: idx_to)
        {
            arma::uword ia = static_cast <arma::uword> (i),
                        ja = static_cast <arma::uword> (j);
            contig_mat (ia, ja) = contig_mat (ja, ia) = 1;
        }

    // then re-number all cluster numbers in cl2index 
    cl2index_map.erase (cluster_from);
    cl2index_map.erase (cluster_to);
    for (auto i: idx_from)
        idx_to.insert (i);
    cl2index_map.emplace (cluster_to, idx_to);
    // and in index2cl:
    for (auto i: idx_from)
    {
        index2cl_map.erase (i);
        index2cl_map.emplace (i, cluster_to);
    }
}

//' initial contiguity and distance matrices. The contiguity matrix is between
//' clusters, so is constantly modified, whereas the distance matrix is between
//' edges, so is fixed at load time.
//' @noRd
void utils_slk::mats_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const int2indx_map_t &vert2index_map,
        arma::Mat <int> &contig_mat,
        arma::Mat <double> &d_mat,
        bool shortest)
{
    // arma::uword = unsigned int
    const arma::uword n = static_cast <arma::uword> (vert2index_map.size ());

    contig_mat = arma::zeros <arma::Mat <int> > (n, n);
    //d_mat = arma::zeros <arma::Mat <double> > (n, n);
    d_mat.resize (n, n);
    if (shortest)
        d_mat.fill (INFINITE_DOUBLE);
    else
        d_mat.fill (-INFINITE_DOUBLE);

    for (int i = 0; i < from.length (); i++)
    {
        arma::uword fi = static_cast <arma::uword> (vert2index_map.at (from [i])),
                    ti = static_cast <arma::uword> (vert2index_map.at (to [i]));
        contig_mat (fi, ti) = contig_mat (ti, fi) = 1;
        d_mat (fi, ti) = d_mat (ti, fi) = d [i];
    }
}
