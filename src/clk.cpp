#include "common.h"
#include "utils.h"
#include "clk.h"

// --------- COMPLETE LINKAGE CLUSTER ----------------

bool edge_sorter (oneEdge const & lhs, oneEdge const & rhs)
{
    return lhs.dist < rhs.dist;
}

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
    clk_dat.edges_all.resize (from_full.size ());
    for (int i = 0; i < from_full.size (); i++)
    {
        oneEdge here;
        here.from = from_full [i];
        here.to = to_full [i];
        here.dist = d_full [i];
        clk_dat.edges_all [i] = here;
    }
    std::sort (clk_dat.edges_all.begin (), clk_dat.edges_all.end (),
            edge_sorter);

    clk_dat.edges_nn.clear ();
    clk_dat.edges_nn.resize (from.size ());
    for (int i = 0; i < from.size (); i++)
    {
        oneEdge here;
        here.from = from [i];
        here.to = to [i];
        here.dist = d [i];
        clk_dat.edges_nn [i] = here;
    }
    std::sort (clk_dat.edges_nn.begin (), clk_dat.edges_nn.end (),
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
        clk_dat.vert2index_map.emplace (v, i++);

    clk_dat.contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    clk_dat.dmax.zeros (n, n);
    for (int i = 0; i < from.length (); i++)
    {
        unsigned int vf = clk_dat.vert2index_map.at (from [i]),
                     vt = clk_dat.vert2index_map.at (to [i]);
        clk_dat.contig_mat (vf, vt) = 1;
        //clk_dat.dmax (vf, vt) = d [i]; // NOPE - all dmax = 0 at start!
    }
}

//' clk_step
//'
//' @param ei The i'th edge of the full sorted list of edge weights
//' @noRd
unsigned int clk_step (CLKDat &clk_dat, unsigned int i)
{
    // find shortest _nn edges that connects the two clusters
    oneEdge ei = clk_dat.edges_all [i];
    unsigned int u = clk_dat.vert2index_map.at (ei.from),
                 v = clk_dat.vert2index_map.at (ei.to);

    // Find shortest edge in MST that connects u and v:
    unsigned int mmin = INFINITE_INT, lmin = INFINITE_INT,
                 the_edge = INFINITE_INT;
    double dmin = INFINITE_DOUBLE;
    for (int j = 0; j < clk_dat.edges_nn.size (); j++)
    {
        oneEdge ej = clk_dat.edges_nn [j];
        unsigned int m = clk_dat.vert2index_map.at (ej.from),
                     l = clk_dat.vert2index_map.at (ej.to);
        if (((clk_dat.index2cl_map.at (m) == u &&
                        clk_dat.index2cl_map.at (l) == v) ||
                    (clk_dat.index2cl_map.at (m) == v &&
                     clk_dat.index2cl_map.at (l) == u)) &&
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

    merge_clusters (clk_dat.contig_mat,
            clk_dat.index2cl_map,
            clk_dat.cl2index_map, mmin, lmin);

    for (auto cl: clk_dat.cl2index_map)
    {
        if (cl.first != lmin || cl.first != mmin)
        {
            const double dl = clk_dat.dmax (cl.first, lmin),
                  dm = clk_dat.dmax (cl.first, mmin);
            double dtemp = dl;
            if (dm < dl)
                dtemp = dm;
            clk_dat.dmax (cl.first, lmin) = dtemp;

            if (clk_dat.contig_mat (cl.first, lmin) == 1 ||
                    clk_dat.contig_mat (cl.first, mmin) == 1)
            {
                clk_dat.contig_mat (cl.first, lmin) = 1;
            } // end if C(c, l) = 1 or C(c, m) = 1 in Guo's terminology
        } // end if cl.first != (cfrom, cto)
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

    CLKDat clk_dat;
    clk_init (clk_dat, from_full, to_full, d_full, from, to, d);

    std::vector <int> treevec;
    for (int i = 0; i < clk_dat.edges_all.size (); i++)
    {
        oneEdge ei = clk_dat.edges_all [i];
        unsigned int m = clk_dat.vert2index_map.at (ei.from),
                     l = clk_dat.vert2index_map.at (ei.to);
        if (clk_dat.index2cl_map.at (l) != clk_dat.index2cl_map.at (m) &&
                clk_dat.contig_mat (l, m) == 1 &&
                ei.dist > clk_dat.dmax (m, l))
        {
            clk_step (clk_dat, i);
        }
    }


    return Rcpp::wrap (treevec);
}
