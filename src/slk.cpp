#include "common.h"
#include "utils.h"
#include "slk.h"
#include <algorithm>

// --------- SINGLE LINKAGE CLUSTER ----------------

//' rcpp_slk
//'
//' Full-order single linkage cluster redcap algorithm
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_slk (
        const Rcpp::DataFrame gr_full,
        const Rcpp::DataFrame gr,
        bool shortest)
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

    arma::Mat <int> contig_mat;
    arma::Mat <double> d_mat, d_mat_full;

    // index2cl and cl2index are dynamically updated with cluster memberships;
    // vert2index and index2vert are retained at initial values which map (from,
    // to) vectors to matrix indices. All operations are performed on matrices
    // directly, with membership re-traced at the end via index2vert_map.
    int2indxset_map_t cl2index_map;
    int2indx_map_t vert2index_map;
    indx2int_map_t index2vert_map, index2cl_map;

    size_t n = utils::sets_init (from, to, vert2index_map, index2vert_map,
                          index2cl_map, cl2index_map);

    utils::mats_init (from, to, d, vert2index_map, contig_mat, d_mat,
            shortest);
    utils::dmat_full_init (from_full, to_full, d_full, vert2index_map,
            d_mat_full, shortest);

    /* The contiguity matrix retains is shape, so is always indexed by the
     * (from, to) vectors. Merging clusters simply switches additional entries
     * from  0 to 1.
     */

    indxset_t the_tree;
    size_t e = 0; // edge number in gr_full
    while (the_tree.size () < (n - 1)) // tree has n - 1 edges
    {
        Rcpp::checkUserInterrupt ();

        index_t ifrom = vert2index_map.at (from_full (e)),
                ito = vert2index_map.at (to_full (e));
        if (index2cl_map.find (ifrom) != index2cl_map.end () &&
                index2cl_map.find (ito) != index2cl_map.end ())
        {
            int cfrom = index2cl_map.at (ifrom),
                cto = index2cl_map.at (ito);
            if (cfrom != cto &&
                    contig_mat (static_cast <arma::uword> (ifrom),
                                static_cast <arma::uword> (ito)) > 0)
            {
                size_t ishort = utils::find_shortest_connection (from, to, d,
                        vert2index_map, d_mat, cl2index_map, cfrom, cto,
                        shortest);
                the_tree.insert (ishort);
                utils::merge_clusters (contig_mat, index2cl_map, cl2index_map,
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
        if (e >= from_full.size () || e >= to_full.size ())
        {
            break;
        }
    }

    std::vector <index_t> treevec (the_tree.begin (), the_tree.end ());

    return Rcpp::wrap (treevec);
}
