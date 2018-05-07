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

    arma::Mat <unsigned short> contig_mat;
    arma::Mat <double> d_mat, d_mat_full;

    // index2cl and cl2index are dynamically updated with cluster memberships;
    // vert2index and index2vert are retained at initial values which map (from,
    // to) vectors to matrix indices. All operations are performed on matrices
    // directly, with membership re-traced at the end via index2vert_map.
    uint_map_t vert2index_map, index2vert_map, index2cl_map;
    uint_set_map_t cl2index_map;

    unsigned int n = sets_init (from, to, vert2index_map, index2vert_map,
            index2cl_map, cl2index_map);

    mats_init (from, to, d, vert2index_map, contig_mat, d_mat);
    dmat_full_init (from_full, to_full, d_full, vert2index_map, d_mat_full);

    /* The contiguity matrix retains is shape, so is always indexed by the
     * (from, to) vectors. Merging clusters simply switches additional entries
     * from  0 to 1.
     */

    std::unordered_set <unsigned int> the_tree;
    unsigned int e = 0; // edge number in gr_full
    while (the_tree.size () < (n - 1)) // tree has n - 1 edges
    {
        unsigned int ff = static_cast <unsigned int> (from_full (e)),
                     tf = static_cast <unsigned int> (to_full (e));
        unsigned int ifrom = vert2index_map.at (ff),
                     ito = vert2index_map.at (tf);
        if (index2cl_map.find (ifrom) != index2cl_map.end () &&
                index2cl_map.find (ito) != index2cl_map.end ())
        {
            unsigned int cfrom = index2cl_map.at (ifrom),
                         cto = index2cl_map.at (ito);
            if (cfrom != cto && contig_mat (ifrom, ito) > 0)
            {
                if (cl2index_map.find (cfrom) == cl2index_map.end ())
                {
                    Rcpp::Rcout << "cfrom [" << ifrom << "] = " << cfrom <<
                        " not in cl2index_map" << std::endl;
                    Rcpp::stop ("shite");
                }
                if (cl2index_map.find (cto) == cl2index_map.end ())
                {
                    Rcpp::Rcout << "cto [" << ito << "] = " << cto <<
                        " not in cl2index_map" << std::endl;
                    Rcpp::stop ("shite");
                }
                unsigned int ishort = find_shortest_connection (from, to, d,
                        vert2index_map, d_mat, cl2index_map, cfrom, cto);
                the_tree.insert (ishort);
                merge_clusters (contig_mat, index2cl_map, cl2index_map,
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
    }

    std::vector <int> treevec (the_tree.begin (), the_tree.end ());

    return Rcpp::wrap (treevec);
}
