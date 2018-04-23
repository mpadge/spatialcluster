#include "common.h"
#include "utils.h"
#include "exact.h"

// --------- EXACT CLUSTER ----------------

void clexact_init (EXDat &clexact_dat,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d)
{
    unsigned int n = sets_init (from, to, clexact_dat.vert2index_map,
            clexact_dat.index2vert_map, clexact_dat.index2cl_map,
            clexact_dat.cl2index_map);
    clexact_dat.n = n;

    clexact_dat.edges_nn.clear ();
    clexact_dat.edges_nn.resize (from.size ());
    for (int i = 0; i < from.size (); i++)
    {
        oneEdge here;
        here.from = from [i];
        here.to = to [i];
        here.dist = d [i];
        clexact_dat.edges_nn [i] = here;
    }
    std::sort (clexact_dat.edges_nn.begin (), clexact_dat.edges_nn.end (),
            edge_sorter);

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
        clexact_dat.vert2index_map.emplace (v, i++);

    clexact_dat.contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    clexact_dat.dmax.zeros (n, n);
    for (int i = 0; i < from.length (); i++)
    {
        unsigned int vf = clexact_dat.vert2index_map.at (from [i]),
                     vt = clexact_dat.vert2index_map.at (to [i]);
        clexact_dat.contig_mat (vf, vt) = 1;
        //clexact_dat.dmax (vf, vt) = d [i]; // NOPE - all dmax = 0 at start!
    }
}

//' clexact_step
//'
//' @param ei The i'th edge of the full sorted list of edge weights
//' @noRd
unsigned int clexact_step (EXDat &clexact_dat, unsigned int i)
{
    unsigned int the_edge = INFINITE_INT;

    return the_edge;
}

//' rcpp_exact
//'
//' Full-order complete linkage cluster redcap algorithm
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_exact (
        const Rcpp::DataFrame gr)
{
    Rcpp::IntegerVector from_ref = gr ["from"];
    Rcpp::IntegerVector to_ref = gr ["to"];
    Rcpp::NumericVector d_ref = gr ["d"];
    Rcpp::IntegerVector from = Rcpp::clone (from_ref);
    Rcpp::IntegerVector to = Rcpp::clone (to_ref);
    Rcpp::NumericVector d = Rcpp::clone (d_ref);

    // Index vectors are 1-indexed, so
    from = from - 1;
    to = to - 1;

    EXDat clexact_dat;
    clexact_init (clexact_dat, from, to, d);

    std::vector <int> treevec;

    return Rcpp::wrap (treevec);
}
