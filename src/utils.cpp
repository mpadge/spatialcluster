#include "common.h"
#include "utils.h"

// Note that all matrices **CAN** be asymmetrical, and so are always indexed
// (from, to)

unsigned int get_n (
    const Rcpp::IntegerVector &from,
    const Rcpp::IntegerVector &to)
{
    std::set <unsigned int> nodes;
    int max_node = 0;
    for (auto i: from)
    {
        nodes.insert (i);
        if (i > max_node)
            max_node = i;
    }
    for (auto i: to)
    {
        nodes.insert (i);
        if (i > max_node)
            max_node = i;
    }
    if (nodes.size () != (static_cast <unsigned int> (max_node) + 1))
        Rcpp::stop ("vertex numbers are discontinuous");

    return nodes.size ();
}

//' initial contiguity and distance matrices. The contiguity matrix is between
//' clusters, so is constantly modified, whereas the distance matrix is between
//' edges, so is fixed at load time.
//' @noRd
void mats_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        arma::Mat <unsigned short> &contig_mat,
        arma::Mat <double> &d_mat,
        const unsigned int n)
{
    contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    //d_mat = arma::zeros <arma::Mat <double> > (n, n);
    d_mat.resize (n, n);
    d_mat.fill (INFINITE_DOUBLE);

    for (int i = 0; i < from.length (); i++)
    {
        contig_mat (from [i], to [i]) = 1;
        d_mat (from [i], to [i]) = d [i];
    }
}

void dmat_full_init (
        const Rcpp::IntegerVector &from, // here, from_full, etc.
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        arma::Mat <double> &d_mat, // here, d_mat_full
        const unsigned int n)
{
    //d_mat = arma::zeros <arma::Mat <double> > (n, n);
    d_mat.resize (n, n);
    d_mat.fill (INFINITE_DOUBLE);

    for (int i = 0; i < from.length (); i++)
    {
        d_mat [from [i], to [i]] = d [i];
    }
}

void sets_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        uint_map_t &vert2cl_map,
        uint_set_map_t &cl2vert_map)
{
    vert2cl_map.clear ();
    cl2vert_map.clear ();

    for (int i = 0; i < from.length (); i++)
    {
        std::set <unsigned int> eset;
        eset.insert (from [i]);
        cl2vert_map.emplace (from [i], eset);
    }
    for (int i = 0; i < to.length (); i++)
    {
        if (cl2vert_map.find (to [i]) == cl2vert_map.end ())
        {
            std::set <unsigned int> eset;
            eset.insert (to [i]);
            cl2vert_map.emplace (to [i], eset);
        } else
        {
            std::set <unsigned int> eset = cl2vert_map.at (to [i]);
            eset.emplace (to [i]);
            cl2vert_map.at (to [i]) = eset;
        }
    }
    
    const unsigned int n = get_n (from, to);
    // Initially assign all verts to clusters of same number:
    for (unsigned int i = 0; i < n; i++)
        vert2cl_map.emplace (i, i);
}

//' find shortest connection between two clusters
//' @param from, to, d the columns of the edge graph
//' @param d_mat distance matrix between all edges (not between clusters!)
//' @param cl2vert_map map of list of all (from, to, d) edges for each cluster
//' @param cfrom Number of cluster which is to be merged
//' @param cto Number of cluster with which it is to be merged
//' @noRd
int find_shortest_connection (
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        Rcpp::NumericVector &d,
        arma::Mat <double> &d_mat,
        uint_set_map_t &cl2vert_map,
        int cfrom,
        int cto)
{
    std::set <unsigned int> verts_i = cl2vert_map.at (cfrom);
    std::set <unsigned int> verts_j = cl2vert_map.at (cto);

    double dmin = INFINITE_DOUBLE;
    int short_i = INFINITE_INT, short_j = INFINITE_INT;

    // from and to here are not direction, so need to examine both directions
    for (auto i: verts_i)
        for (auto j: verts_j)
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
        if (from [i] == short_i && to [i] == short_j)
        {
            shortest = i;
            break;
        }
    }
    if (shortest == INFINITE_INT)
        Rcpp::stop ("shite");

    return shortest;
}
