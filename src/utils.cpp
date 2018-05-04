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

bool strfound (const std::string str, const std::string target)
{
    bool found = false;
    if (str.find (target) != std::string::npos)
        found = true;
    return found;
}

unsigned int sets_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        uint_map_t &vert2index_map,
        uint_map_t &index2vert_map,
        uint_map_t &index2cl_map,
        uint_set_map_t &cl2index_map)
{
    vert2index_map.clear ();
    index2vert_map.clear ();
    index2cl_map.clear ();
    cl2index_map.clear ();

    std::unordered_set <unsigned int> vert_set;
    for (int i = 0; i < from.size (); i++)
    {
        vert_set.emplace (from [i]);
        vert_set.emplace (to [i]);
    }
    unsigned int i = 0;
    for (auto v: vert_set)
    {
        index2vert_map.emplace (i, v);
        vert2index_map.emplace (v, i++);
    }

    std::unordered_set <unsigned int> eset;
    for (int i = 0; i < from.length (); i++)
    {
        unsigned int fi = vert2index_map.at (from [i]);
        eset.clear ();
        eset.insert (fi);
        cl2index_map.emplace (fi, eset);
    }
    for (int i = 0; i < to.length (); i++)
    {
        unsigned int ti = vert2index_map.at (to [i]);
        if (cl2index_map.find (ti) == cl2index_map.end ())
            eset.clear ();
        else
            eset = cl2index_map.at (ti);
        eset.emplace (ti);
        cl2index_map.emplace (ti, eset);
    }
    
    const unsigned int n = vert_set.size ();
    // Initially assign all verts to clusters of same number:
    for (unsigned int i = 0; i < n; i++)
        index2cl_map.emplace (i, i);

    return n;
}

//' initial contiguity and distance matrices. The contiguity matrix is between
//' clusters, so is constantly modified, whereas the distance matrix is between
//' edges, so is fixed at load time.
//' @noRd
void mats_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const uint_map_t &vert2index_map,
        arma::Mat <unsigned short> &contig_mat,
        arma::Mat <double> &d_mat)
{
    const unsigned int n = vert2index_map.size ();

    contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    //d_mat = arma::zeros <arma::Mat <double> > (n, n);
    d_mat.resize (n, n);
    d_mat.fill (INFINITE_DOUBLE);

    for (int i = 0; i < from.length (); i++)
    {
        unsigned int fi = vert2index_map.at (from [i]),
                     ti = vert2index_map.at (to [i]);
        contig_mat (fi, ti) = 1;
        d_mat (fi, ti) = d [i];
    }
}

void dmat_full_init (
        const Rcpp::IntegerVector &from, // here, from_full, etc.
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        const uint_map_t &vert2index_map,
        arma::Mat <double> &d_mat) // here, d_mat_full
{
    //d_mat = arma::zeros <arma::Mat <double> > (n, n);
    d_mat.resize (vert2index_map.size (), vert2index_map.size ());
    d_mat.fill (INFINITE_DOUBLE);

    for (int i = 0; i < from.length (); i++)
    {
        d_mat [vert2index_map.at (from [i]),
              vert2index_map.at (to [i])] = d [i];
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
int find_shortest_connection (
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        Rcpp::NumericVector &d,
        uint_map_t &vert2index_map,
        arma::Mat <double> &d_mat,
        uint_set_map_t &cl2index_map,
        const unsigned int cfrom,
        const unsigned int cto)
{
    if (cl2index_map.find (cfrom) == cl2index_map.end ())
        Rcpp::stop ("cluster index not found");
    if (cl2index_map.find (cto) == cl2index_map.end ())
        Rcpp::stop ("cluster index not found");
    std::unordered_set <unsigned int> index_i = cl2index_map.at (cfrom),
        index_j = cl2index_map.at (cto);

    double dmin = INFINITE_DOUBLE;
    int short_i = INFINITE_INT, short_j = INFINITE_INT;

    // from and to here are not directional, so need to examine both directions
    for (auto i: index_i)
        for (auto j: index_j)
        {
            if (d_mat (i, j) < dmin)
            {
                dmin = d_mat (i, j);
                short_i = i;
                short_j = j;
            } else if (d_mat (j, i) < dmin)
            {
                dmin = d_mat (j, i);
                short_i = j;
                short_j = i;
            }
        }
    if (dmin == INFINITE_DOUBLE)
        Rcpp::stop ("no minimal distance; this should not happen");

    // convert short_i and short_j to a single edge 
    // TODO: Make a std::map of vert2dist to avoid this loop
    int shortest = INFINITE_INT;
    for (int i = 0; i < from.length (); i++)
    {
        if (vert2index_map.at (from [i]) == short_i &&
                vert2index_map.at (to [i]) == short_j)
        {
            shortest = i;
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
void merge_clusters (
        arma::Mat <unsigned short> &contig_mat,
        uint_map_t &index2cl_map,
        uint_set_map_t &cl2index_map,
        const unsigned int cluster_from,
        const unsigned int cluster_to)
{
    // Set all contig_mat (cluster_from, .) to 1
    for (unsigned int i = 0; i < contig_mat.n_rows; i++)
    {
        if (contig_mat (cluster_from, i) == 1 )
        {
            contig_mat (cluster_to, i) = 1;
            contig_mat (i, cluster_to) = 1;
        }
    }

    std::unordered_set <unsigned int>
        idx_from = cl2index_map.at (cluster_from),
        idx_to = cl2index_map.at (cluster_to);

    // not directonal here, so need both directions:
    for (auto i: idx_from)
        for (auto j: idx_to)
        {
            contig_mat (i, j) = contig_mat (j, i) = 1;
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
