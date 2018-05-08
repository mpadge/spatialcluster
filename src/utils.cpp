#include "common.h"
#include "utils.h"

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
    int i = 0;
    for (auto v: vert_set)
    {
        index2vert_map.emplace (i, v);
        vert2index_map.emplace (v, i++);
    }

    indxset_t eset;
    for (int i = 0; i < from.length (); i++)
    {
        size_t fi = vert2index_map.at (from [i]);
        eset.clear ();
        eset.insert (fi);
        cl2index_map.emplace (fi, eset);
    }
    for (int i = 0; i < to.length (); i++)
    {
        // all verts are their own clusters, so cast indx to cli
        int cli = static_cast <int> (vert2index_map.at (to [i]));
        if (cl2index_map.find (cli) == cl2index_map.end ())
            eset.clear ();
        else
            eset = cl2index_map.at (cli);
        eset.emplace (cli);
        cl2index_map.emplace (cli, eset);
    }
    
    const int n = static_cast <int> (vert_set.size ());
    // Initially assign all verts to clusters of same number:
    for (int i = 0; i < n; i++)
        index2cl_map.emplace (static_cast <size_t> (i), i);

    return static_cast <size_t> (n);
}

//' initial contiguity and distance matrices. The contiguity matrix is between
//' clusters, so is constantly modified, whereas the distance matrix is between
//' edges, so is fixed at load time.
//' @noRd
void utils::mats_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const int2indx_map_t &vert2index_map,
        arma::Mat <int> &contig_mat,
        arma::Mat <double> &d_mat)
{
    // arma::uword = unsigned int
    const arma::uword n = static_cast <arma::uword> (vert2index_map.size ());

    contig_mat = arma::zeros <arma::Mat <int> > (n, n);
    //d_mat = arma::zeros <arma::Mat <double> > (n, n);
    d_mat.resize (n, n);
    d_mat.fill (INFINITE_DOUBLE);

    for (int i = 0; i < from.length (); i++)
    {
        arma::uword fi = static_cast <arma::uword> (vert2index_map.at (from [i])),
                    ti = static_cast <arma::uword> (vert2index_map.at (to [i]));
        contig_mat (fi, ti) = 1;
        d_mat (fi, ti) = d [i];
    }
}

void utils::dmat_full_init (
        const Rcpp::IntegerVector &from, // here, from_full, etc.
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const int2indx_map_t &vert2index_map,
        arma::Mat <double> &d_mat) // here: d_mat_full
{
    //d_mat = arma::zeros <arma::Mat <double> > (n, n);
    const arma::uword n = static_cast <arma::uword> (vert2index_map.size ());
    d_mat.resize (n, n);
    d_mat.fill (INFINITE_DOUBLE);

    for (int i = 0; i < from.length (); i++)
    {
        d_mat [static_cast <arma::uword> (vert2index_map.at (from [i])),
              static_cast <arma::uword> (vert2index_map.at (to [i]))] = d [i];
    }
}

//' find shortest connection between two clusters
//' @param from, to, d the columns of the edge graph
//' @param d_mat distance matrix between all edges (not between clusters!)
//' @param cl2vert_map map of list of all (from, to, d) edges for each cluster
//' @param cfrom Number of cluster which is to be merged
//' @param cto Number of cluster with which it is to be merged
//'
//' @return Index directly into from, to - **NOT** into the actual matrices!
//' @noRd
size_t utils::find_shortest_connection (
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        Rcpp::NumericVector &d,
        int2indx_map_t &vert2index_map,
        arma::Mat <double> &d_mat,
        int2indxset_map_t &cl2index_map,
        const int cfrom,
        const int cto)
{
    if (cl2index_map.find (cfrom) == cl2index_map.end ())
        Rcpp::stop ("cluster index not found");
    if (cl2index_map.find (cto) == cl2index_map.end ())
        Rcpp::stop ("cluster index not found");
    indxset_t index_i = cl2index_map.at (cfrom),
             index_j = cl2index_map.at (cto);

    double dmin = INFINITE_DOUBLE;
    size_t short_i = INFINITE_INT, short_j = INFINITE_INT;

    // from and to here are not directional, so need to examine both directions
    for (auto i: index_i)
        for (auto j: index_j)
        {
            arma::uword ia = static_cast <arma::uword> (i),
                        ja = static_cast <arma::uword> (j);
            if (d_mat (ia, ja) < dmin)
            {
                dmin = d_mat (ia, ja);
                short_i = i;
                short_j = j;
            } else if (d_mat (ja, ia) < dmin)
            {
                dmin = d_mat (ja, ia);
                short_i = j;
                short_j = i;
            }
        }
    if (dmin == INFINITE_DOUBLE)
        Rcpp::stop ("no minimal distance; this should not happen");

    // convert short_i and short_j to a single edge 
    // TODO: Make a std::map of vert2dist to avoid this loop
    size_t shortest = INFINITE_INT;
    for (int i = 0; i < from.length (); i++) // int for Rcpp index
    {
        if (vert2index_map.at (from [i]) == short_i &&
                vert2index_map.at (to [i]) == short_j)
        {
            shortest = static_cast <size_t> (i);
            break;
        }
    }
    if (shortest == INFINITE_INT)
        Rcpp::stop ("shite");

    return shortest;
}

//' merge two clusters in the contiguity matrix, reducing the size of the matrix
//' by one row and column.
//' @noRd
void utils::merge_clusters (
        arma::Mat <int> &contig_mat,
        indx2int_map_t &index2cl_map,
        int2indxset_map_t &cl2index_map,
        int cluster_from,
        int cluster_to)
{
    if (cluster_from < 0)
        Rcpp::stop ("cluster_from must be non-zero");
    if (cluster_to < 0)
        Rcpp::stop ("cluster_to must be non-zero");
    arma::uword cfr = static_cast <arma::uword> (cluster_from),
                cto = static_cast <arma::uword> (cluster_to);
    // Set all contig_mat (cluster_from, .) to 1
    for (arma::uword i = 0; i < contig_mat.n_rows; i++)
    {
        if (contig_mat (cfr, i) == 1 )
        {
            contig_mat (cto, i) = 1;
            contig_mat (i, cto) = 1;
        }
    }

    indxset_t idx_from = cl2index_map.at (cluster_from),
              idx_to = cl2index_map.at (cluster_to);

    // not directonal here, so need both directions:
    for (auto i: idx_from)
        for (auto j: idx_to)
        {
            arma::uword ia = static_cast <arma::uword> (i),
                        ja = static_cast <arma::uword> (j);
            contig_mat (ia, ja) = contig_mat (ja, ia) = 1;
        }

    // then re-number all cluster numbers in cl2index 
    cl2index_map.erase (cluster_from);
    for (auto i: idx_from)
        idx_to.insert (i);
    cl2index_map [cluster_to] = idx_to;
    // and in index2cl:
    for (auto i: idx_from)
        index2cl_map [i] = cluster_to;
}
