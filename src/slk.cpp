#include "common.h"
#include "utils.h"
#include "slk.h"

// --------- SINGLE LINKAGE CLUSTER ----------------

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

//' merge two clusters in the contiguity matrix, reducing the size of the matrix
//' by one row and column.
//' @noRd
void merge_clusters (
        arma::Mat <unsigned short> &contig_mat,
        uint_map_t &vert2cl_map,
        uint_set_map_t &cl2vert_map,
        int cluster_from,
        int cluster_to)
{
    // Set all contig_mat (cluster_from, .) to 1
    for (unsigned int j = 0; j < contig_mat.n_rows; j++)
    {
        if (contig_mat (cluster_from, j) == 1 ||
                contig_mat (j, cluster_from) == 1)
        {
            contig_mat (cluster_to, j) = 1;
            contig_mat (j, cluster_to) = 1;
        }
    }

    std::set <unsigned int> verts_from = cl2vert_map.at (cluster_from),
        verts_to = cl2vert_map.at (cluster_to);

    for (auto vi: verts_from)
        for (auto vj: verts_to)
        {
            contig_mat (vi, vj) = contig_mat (vj, vi) = 1;
        }

    // then re-number all cluster numbers in cl2vert 
    cl2vert_map.erase (cluster_from);
    for (auto v: verts_from)
        verts_to.insert (v);
    cl2vert_map.at (cluster_to) = verts_to;
    // and in vert2cl:
    for (auto v: verts_from)
        vert2cl_map [v] = cluster_to;
}

//' rcpp_slk
//'
//' Full-order single linkage cluster redcap algorithm
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_slk (
        const Rcpp::DataFrame gr_full,
        const Rcpp::DataFrame gr)
{
    Rcpp::IntegerVector from_full_ref = gr_full ["from"];
    Rcpp::IntegerVector to_full_ref = gr_full ["to"];
    Rcpp::NumericVector d_full = gr_full ["d"];
    Rcpp::IntegerVector from_ref = gr ["from"];
    Rcpp::IntegerVector to_ref = gr ["to"];
    Rcpp::NumericVector d = gr ["d"];

    // Rcpp classes are always passed by reference, so cloning is necessary to
    // avoid modifying the original data.frames.
    Rcpp::IntegerVector from_full = Rcpp::clone (from_full_ref);
    Rcpp::IntegerVector to_full = Rcpp::clone (to_full_ref);
    Rcpp::IntegerVector from = Rcpp::clone (from_ref);
    Rcpp::IntegerVector to = Rcpp::clone (to_ref);

    // Index vectors are 1-indexed, so
    from_full = from_full - 1;
    to_full = to_full - 1;
    from = from - 1;
    to = to - 1;

    const unsigned int n = get_n (from, to);
    const unsigned int nf = get_n (from_full, to_full);

    arma::Mat <unsigned short> contig_mat;
    arma::Mat <double> d_mat, d_mat_full;
    uint_map_t vert2cl_map;
    uint_set_map_t cl2vert_map;

    mats_init (from, to, d, contig_mat, d_mat, n);
    dmat_full_init (from_full, to_full, d_full, d_mat_full, nf);
    sets_init (from, to, vert2cl_map, cl2vert_map);

    /* The contiguity matrix retains is shape, so is always indexed by the
     * (from, to) vectors. Merging clusters simply switches additional entries
     * from  0 to 1.
     */

    std::unordered_set <unsigned int> the_tree;
    int e = 0; // edge number in gr_full
    while (the_tree.size () < (n - 1)) // tree has n - 1 edges
    {
        int vfrom = from_full (e), vto = to_full (e); // vertex numbers
        if (vert2cl_map.find (vfrom) != vert2cl_map.end () &&
                vert2cl_map.find (vto) != vert2cl_map.end ())
        {
            int cfrom = vert2cl_map.at (vfrom), cto = vert2cl_map.at (vto);
            if (cfrom != cto && contig_mat (vfrom, vto) > 0)
            {
                unsigned int ishort = find_shortest_connection (from, to, d, d_mat,
                        cl2vert_map, cfrom, cto);
                the_tree.insert (ishort);
                merge_clusters (contig_mat, vert2cl_map, cl2vert_map,
                        cfrom, cto);
                e = 0;
            } else
            {
                e++;
            }
        } else
        {
            e++;
        }
        //if (e == from_full.length ())
        //    break;
    }

    std::vector <int> treevec (the_tree.begin (), the_tree.end ());

    return Rcpp::wrap (treevec);
}
