#include "common.h"
#include "utils.h"
#include "exact.h"

// --------- EXACT CLUSTER ----------------

void clexact_init (EXDat &clexact_dat,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d)
{
    std::unordered_set <unsigned int> vert_set;
    for (int i = 0; i < from.size (); i++)
    {
        vert_set.emplace (from [i]);
        vert_set.emplace (to [i]);
    }
    clexact_dat.n = vert_set.size ();

    unsigned int i = 0;
    for (auto v: vert_set)
    {
        clexact_dat.index2vert_map.emplace (i, v);
        clexact_dat.vert2index_map.emplace (v, i);
        clexact_dat.index2cl_map.emplace (i++, -1);
    }

    clexact_dat.edges.clear ();
    clexact_dat.edges.resize (from.size ());
    for (int i = 0; i < from.size (); i++)
    {
        oneEdge here;
        here.from = from [i];
        here.to = to [i];
        here.dist = d [i];
        clexact_dat.edges [i] = here;
    }

    clexact_dat.index_in_cluster.resize (clexact_dat.n);
    std::fill (clexact_dat.index_in_cluster.begin (),
            clexact_dat.index_in_cluster.end (), false);
}

void clexact_merge (EXDat &clexact_dat, unsigned int m, unsigned int l)
{
}

//' clexact_step
//'
//' All edges are initially in their own clusters. This merges edge#i with the
//' edge having the 
//'
//' @param ei The i'th edge of the sorted list of NN edge weights
//' @noRd
unsigned int clexact_step (EXDat &clexact_dat, unsigned int i)
{
    unsigned int the_edge = INFINITE_INT;

    /*
    oneEdge ei = clexact_dat.edges [i];
    const unsigned int m = clexact_dat.vert2index_map.at (ei.from),
          l = clexact_dat.vert2index_map.at (ei.to),
          cl_m = clexact_dat.index2cl_map.at (m),
          cl_l = clexact_dat.index2cl_map.at (l);

    if (cl_m != cl_l)
    {
        merge_clusters (clexact_dat.contig_mat,
                clexact_dat.index2cl_map,
                clexact_dat.cl2index_map,
                clexact_dat.index2cl_map.at (m),
                clexact_dat.index2cl_map.at (l));
    }
    */

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

    // -------- Assign first edge
    unsigned int clnum = 1, ei = 0;
    oneEdge edge = clexact_dat.edges [ei];
    unsigned int ito = clexact_dat.vert2index_map.at (edge.to),
                 ifrom = clexact_dat.vert2index_map.at (edge.from);
    clexact_dat.index2cl_map.emplace (ito, clnum);
    clexact_dat.index2cl_map.emplace (ifrom, clnum);
    std::unordered_set <unsigned int> cli;
    cli.insert (ito);
    cli.insert (ifrom);
    clexact_dat.cl2index_map.emplace (clnum, cli);

    clexact_dat.index_in_cluster [ito] =
        clexact_dat.index_in_cluster [ifrom] = true;

    clnum++;
    ei++;
    unsigned int nleft = clexact_dat.n - 2; // number of !index_in_cluster
    // -------- finish assign first edge

    while (nleft > 0)
    {
        bool from_in = false, to_in = false;
        edge = clexact_dat.edges [ei];
        ito = clexact_dat.vert2index_map.at (edge.to);
        ifrom = clexact_dat.vert2index_map.at (edge.from);
        if (clexact_dat.index_in_cluster [ito])
            to_in = true;
        if (clexact_dat.index_in_cluster [ifrom])
            from_in = true;

        if (from_in && to_in)
            continue; // both already allocated, so don't do anything

        unsigned int clnum_i = clnum;
        if (from_in)
        {
            clnum_i = clexact_dat.index2cl_map [ifrom];
            clexact_dat.index_in_cluster [ifrom] = true;
        } else if (to_in)
        {
            clnum_i = clexact_dat.index2cl_map [ito];
            clexact_dat.index_in_cluster [ito] = true;
        }

        clexact_dat.index2cl_map [ifrom] =
            clexact_dat.index2cl_map [ito] = clnum_i;
        if (clexact_dat.cl2index_map.find (clnum_i) ==
                clexact_dat.cl2index_map.end ())
            cli.clear ();
        else
            cli = clexact_dat.cl2index_map.at (clnum_i);
        cli.insert (ito);
        cli.insert (ifrom);
        clexact_dat.cl2index_map [clnum_i] = cli;

        nleft--;
        if (clnum_i == clnum)
            clnum++;
    }

    // Then construct vector mapping edges to cluster numbers
    std::vector <int> clvec (from.size ());
    for (unsigned int i = 0; i < from.size (); i++)
    {
        unsigned int index = clexact_dat.vert2index_map.at (from [i]);
        unsigned int clnum_i = clexact_dat.index2cl_map.at (index);
        std::unordered_set <unsigned int> cl_members =
            clexact_dat.cl2index_map.at (clnum_i);
        for (auto cl: cl_members)
        {
            clvec [cl] = clnum_i;
        }
    }
    
    return Rcpp::wrap (clvec);
}
