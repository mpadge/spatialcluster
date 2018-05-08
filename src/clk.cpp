#include "common.h"
#include "clk.h"

// --------- COMPLETE LINKAGE CLUSTER ----------------

void clk::clk_init (clk::CLKDat &clk_dat,
        Rcpp::IntegerVector from_full,
        Rcpp::IntegerVector to_full,
        Rcpp::NumericVector d_full,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d)
{
    size_t n = utils::sets_init (from, to, clk_dat.vert2index_map,
            clk_dat.index2vert_map, clk_dat.index2cl_map,
            clk_dat.cl2index_map);
    clk_dat.n = n;

    clk_dat.edges_all.clear ();
    clk_dat.edges_all.resize (static_cast <size_t> (from_full.size ()));
    for (int i = 0; i < from_full.size (); i++)
    {
        utils::OneEdge here;
        here.from = from_full [i];
        here.to = to_full [i];
        here.dist = d_full [i];
        clk_dat.edges_all [static_cast <size_t> (i)] = here;
    }
    // These edges are already passed in sorted form, so no need to explicitly
    // sort here.

    clk_dat.edges_nn.clear ();
    clk_dat.edges_nn.resize (static_cast <size_t> (from.size ()));
    for (int i = 0; i < from.size (); i++)
    {
        utils::OneEdge here;
        here.from = from [i];
        here.to = to [i];
        here.dist = d [i];
        clk_dat.edges_nn [static_cast <size_t> (i)] = here;
    }

    // Get set of unique vertices, and store binary tree of edge distances
    intset_t vert_set;
    for (int i = 0; i < from.size (); i++)
    {
        vert_set.emplace (from [i]);
        vert_set.emplace (to [i]);
    }
    // Construct vert2index_map to map each unique vertex to an index
    index_t i = 0;
    for (auto v: vert_set)
        clk_dat.vert2index_map.emplace (v, i++);

    arma::uword nu = static_cast <arma::uword> (n);
    clk_dat.contig_mat = arma::zeros <arma::Mat <int> > (nu, nu);
    clk_dat.dmax.zeros (nu, nu);
    for (int i = 0; i < from.length (); i++)
    {
        arma::uword vf = static_cast <arma::uword> (
                                        clk_dat.vert2index_map.at (from [i])),
                    vt = static_cast <arma::uword> (
                                        clk_dat.vert2index_map.at (to [i]));
        clk_dat.contig_mat (vf, vt) = 1;
        //clk_dat.dmax (vf, vt) = d [i]; // NOPE - all dmax = 0 at start!
    }
}

//' clk_step
//'
//' @param ei The i'th edge of the full sorted list of edge weights
//' @noRd
size_t clk::clk_step (clk::CLKDat &clk_dat, size_t i)
{
    // find shortest _all edges that connects the two clusters
    utils::OneEdge ei = clk_dat.edges_all [i];
    const size_t u = clk_dat.vert2index_map.at (ei.from),
                 v = clk_dat.vert2index_map.at (ei.to);
    const int cl_u = clk_dat.index2cl_map.at (u),
              cl_v = clk_dat.index2cl_map.at (v);

    // Find shortest edge in MST that connects u and v:
    size_t mmin = INFINITE_INT, lmin = INFINITE_INT, the_edge = INFINITE_INT;
    double dmin = INFINITE_DOUBLE;
    for (size_t j = 0; j < clk_dat.edges_nn.size (); j++)
    {
        utils::OneEdge ej = clk_dat.edges_nn [j];
        size_t m = clk_dat.vert2index_map.at (ej.from),
               l = clk_dat.vert2index_map.at (ej.to);
        if (((clk_dat.index2cl_map.at (m) == cl_u &&
                        clk_dat.index2cl_map.at (l) == cl_v) ||
                    (clk_dat.index2cl_map.at (m) == cl_v &&
                     clk_dat.index2cl_map.at (l) == cl_u)) &&
                ej.dist < dmin)
        {
            the_edge = j;
            mmin = m;
            lmin = l;
            dmin = ej.dist;
        }
    }
    if (dmin == INFINITE_DOUBLE)
        Rcpp::stop ("minimal distance not able to be found");

    utils::merge_clusters (clk_dat.contig_mat,
            clk_dat.index2cl_map,
            clk_dat.cl2index_map,
            clk_dat.index2cl_map.at (mmin),
            clk_dat.index2cl_map.at (lmin));

    for (auto cl: clk_dat.cl2index_map)
    {
        if (cl.first != static_cast <int> (lmin) ||
                cl.first != static_cast <int> (mmin))
        {
            arma::uword clu = static_cast <arma::uword> (cl.first),
                        lu = static_cast <arma::uword> (lmin),
                        mu = static_cast <arma::uword> (mmin);
            const double dl = clk_dat.dmax (clu, lu),
                  dm = clk_dat.dmax (clu, mu);
            double dtemp = dl;
            if (dm < dl)
                dtemp = dm;
            clk_dat.dmax (clu, lu) = dtemp;

            if (clk_dat.contig_mat (clu, lu) == 1 ||
                    clk_dat.contig_mat (clu, mu) == 1)
            {
                clk_dat.contig_mat (clu, lu) = 1;
            }
        }
    } // end for over cl

    return the_edge;
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

    clk::CLKDat clk_dat;
    clk::clk_init (clk_dat, from_full, to_full, d_full, from, to, d);

    std::vector <size_t> treevec;
    for (size_t i = 0; i < clk_dat.edges_all.size (); i++)
    {
        utils::OneEdge ei = clk_dat.edges_all [i];
        arma::uword u = static_cast <arma::uword> (
                                    clk_dat.vert2index_map.at (ei.from)),
                    v = static_cast <arma::uword> (
                                    clk_dat.vert2index_map.at (ei.to));
        if (clk_dat.index2cl_map.at (u) != clk_dat.index2cl_map.at (v) &&
                clk_dat.contig_mat (u, v) == 1 &&
                ei.dist > clk_dat.dmax (u, v))
        {
            size_t the_edge = clk_step (clk_dat, i);
            treevec.push_back (the_edge);
        }
    }

    // treevec here in an index into a **sorted** version of (from, to , d)
    return Rcpp::wrap (treevec);
}
