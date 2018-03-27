#include "common.h"
#include "utils.h"
#include "slk.h"

// --------- SINGLE LINKAGE CLUSTER ----------------

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
