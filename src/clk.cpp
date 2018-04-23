#include "common.h"
#include "utils.h"
#include "clk.h"

// --------- COMPLETE LINKAGE CLUSTER ----------------

void clk_init (CLKDat &clk_dat,
        Rcpp::IntegerVector from_full,
        Rcpp::IntegerVector to_full,
        Rcpp::NumericVector d_full,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d)
{
    unsigned int n = sets_init (from, to, clk_dat.vert2index_map,
            clk_dat.index2vert_map, clk_dat.index2cl_map,
            clk_dat.cl2index_map);
    clk_dat.n = n;

    clk_dat.edges_all.clear ();
    clk_dat.edges_all.reserve (d_full.size ());
    for (auto di: d_full)
        clk_dat.edges_all.push_back (di);
    std::sort (clk_dat.edges_all.begin (), clk_dat.edges_all.end ());

    clk_dat.edges_nn.clear ();
    clk_dat.edges_nn.reserve (d.size ());
    for (auto di: d)
        clk_dat.edges_nn.push_back (di);
    std::sort (clk_dat.edges_nn.begin (), clk_dat.edges_nn.end ());

    // Get set of unique vertices, and store binary tree of edge distances
    std::unordered_set <unsigned int> vert_set;
    for (int i = 0; i < from.size (); i++)
    {
        vert_set.emplace (from [i]);
        vert_set.emplace (to [i]);
    }
    // Construct vert2index_map to map each unique vertex to an index
    unsigned int i = 0;
    for (auto v: vert_set)
        clk_dat.vert2index_map.emplace (v, i++);

    clk_dat.contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    clk_dat.dmax.set_size (n, n);
    clk_dat.dmax.fill (INFINITE_DOUBLE);
    for (int i = 0; i < from.length (); i++)
    {
        unsigned int vf = clk_dat.vert2index_map.at (from [i]),
                     vt = clk_dat.vert2index_map.at (to [i]);
        clk_dat.contig_mat (vf, vt) = 1;
        clk_dat.dmax (vf, vt) = d [i];
    }
}


//' rcpp_clk
//'
//' Full-order complete linkage cluster redcap algorithm
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_clk (
        const Rcpp::DataFrame gr_full,
        const Rcpp::DataFrame gr)
{
    Rcpp::IntegerVector from_full_ref = gr_full ["from"];
    Rcpp::IntegerVector to_full_ref = gr_full ["to"];
    Rcpp::NumericVector d_full_ref = gr_full ["d"];
    Rcpp::IntegerVector from_ref = gr ["from"];
    Rcpp::IntegerVector to_ref = gr ["to"];
    Rcpp::NumericVector d_ref = gr ["d"];

    // Rcpp classes are always passed by reference, so cloning is necessary to
    // avoid modifying the original data.frames.
    Rcpp::IntegerVector from_full = Rcpp::clone (from_full_ref);
    Rcpp::IntegerVector to_full = Rcpp::clone (to_full_ref);
    Rcpp::NumericVector d_full = Rcpp::clone (d_full_ref);
    Rcpp::IntegerVector from = Rcpp::clone (from_ref);
    Rcpp::IntegerVector to = Rcpp::clone (to_ref);
    Rcpp::NumericVector d = Rcpp::clone (d_ref);

    // Index vectors are 1-indexed, so
    from_full = from_full - 1;
    to_full = to_full - 1;
    from = from - 1;
    to = to - 1;

    CLKDat clk_dat;
    clk_init (clk_dat, from_full, to_full, d_full, from, to, d);

    std::vector <int> treevec;

    return Rcpp::wrap (treevec);
}
