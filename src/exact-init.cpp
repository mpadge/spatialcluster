#include "common.h"
#include "utils.h"
#include "exact-init.h"

// --------- EXACT CLUSTER ----------------

void exact_init::clexact_init (exact_init::EXDat &clexact_dat,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d)
{
    intset_t vert_set;
    for (int i = 0; i < from.size (); i++)
    {
        vert_set.emplace (from [i]);
        vert_set.emplace (to [i]);
    }
    clexact_dat.n = vert_set.size ();

    index_t i = 0;
    for (auto v: vert_set)
    {
        clexact_dat.index2vert_map.emplace (i, v);
        clexact_dat.vert2index_map.emplace (v, i);
        clexact_dat.index2cl_map.emplace (i++, -1);
    }

    clexact_dat.edges.clear ();
    clexact_dat.edges.resize (static_cast <size_t> (from.size ()));
    for (int i = 0; i < from.size (); i++)
    {
        utils::OneEdge here;
        here.from = from [i];
        here.to = to [i];
        here.dist = d [i];
        clexact_dat.edges [static_cast <size_t> (i)] = here;
    }

    clexact_dat.index_in_cluster.resize (clexact_dat.n);
    std::fill (clexact_dat.index_in_cluster.begin (),
            clexact_dat.index_in_cluster.end (), false);
}


void exact_init::assign_first_edge (exact_init::EXDat &clexact_dat)
{
    int clnum = 0;
    index_t ei = 0;
    utils::OneEdge edge = clexact_dat.edges [ei];
    index_t ito = clexact_dat.vert2index_map.at (edge.to),
            ifrom = clexact_dat.vert2index_map.at (edge.from);

    clexact_dat.index2cl_map [ito] = clnum;
    clexact_dat.index2cl_map [ifrom] = clnum;
    clexact_dat.index_in_cluster [ito] = true;
    clexact_dat.index_in_cluster [ifrom] = true;

    intset_t cli;
    cli.insert (static_cast <int> (ito));
    cli.insert (static_cast <int> (ifrom));
    clexact_dat.cl2index_map.emplace (clnum, cli);

    clexact_dat.vert2cl_map.emplace (edge.to, clnum);
    clexact_dat.vert2cl_map.emplace (edge.from, clnum);
}


//' clexact_step
//'
//' All edges are initially in their own clusters. This merges edge#i with the
//' next closest edge
//'
//' @param ei The i'th edge of the sorted list of NN edge weights
//' @noRd
int exact_init::clexact_step (exact_init::EXDat &clexact_dat,
        const index_t ei, const int clnum)
{
    bool from_in = false, to_in = false;
    utils::OneEdge edge = clexact_dat.edges [ei];
    index_t ito = clexact_dat.vert2index_map.at (edge.to),
            ifrom = clexact_dat.vert2index_map.at (edge.from);
    if (clexact_dat.index_in_cluster [ito])
        to_in = true;
    if (clexact_dat.index_in_cluster [ifrom])
        from_in = true;

    int clnum_i = clnum;

    if (from_in && to_in)
    {
        clnum_i = INFINITE_INT;
    } else
    {
        if (from_in) // then to is not in cluster
            clnum_i = clexact_dat.index2cl_map [ifrom];
        else if (to_in)
            clnum_i = clexact_dat.index2cl_map [ito];

        clexact_dat.index_in_cluster [ito] = true;
        clexact_dat.index_in_cluster [ifrom] = true;
        clexact_dat.index2cl_map [ifrom] = clnum_i;
        clexact_dat.index2cl_map [ito] = clnum_i;

        intset_t cli;
        if (clexact_dat.cl2index_map.find (clnum_i) !=
                clexact_dat.cl2index_map.end ())
            cli = clexact_dat.cl2index_map.at (clnum_i);
        cli.insert (static_cast <int> (ito));
        cli.insert (static_cast <int> (ifrom));
        clexact_dat.cl2index_map [clnum_i] = cli;

        // These values may already be in the map here, but that's okay
        clexact_dat.vert2cl_map.emplace (edge.to, clnum_i);
        clexact_dat.vert2cl_map.emplace (edge.from, clnum_i);
    }

    return clnum_i;
}

//' fill_cl_edges
//'
//' Fill (arma) matrix of strongest/shortest connections between all clusters
//' used to construct the hierarchical relationships
//' @noRd
void exact_init::fill_cl_edges (exact_init::EXDat &clexact_dat,
        arma::Mat <double> &cl_edges, int num_clusters)
{
    int2intset_map_t vert_sets;
    for (int i = 0; i < num_clusters; i++)
    {
        intset_t verts;
        for (auto vi: clexact_dat.vert2cl_map)
        {
            if (vi.second == i)
            {
                verts.emplace (vi.first);
            }
        }
        vert_sets.emplace (i, verts);
    }

    // need a (sparse) matrix of all pairwise edge distances:
    arma::uword nu = utils::to_uword (clexact_dat.n);
    arma::Mat <double> vert_dists (nu, nu);
    for (auto ei: clexact_dat.edges)
    {
        arma::uword i = utils::to_uword (clexact_dat.vert2index_map.at (ei.from)),
                    j = utils::to_uword (clexact_dat.vert2index_map.at (ei.to));
        vert_dists (i, j) = vert_dists (j, i) = ei.dist;
    }

    for (int i = 0; i < (num_clusters - 1); i++)
        for (int j = (i + 1); j < num_clusters; j++)
        {
            intset_t verts_i = vert_sets.at (i),
                     verts_j = vert_sets.at (j);
            double max_d = 0.0;
            for (auto vi: verts_i)
                for (auto vj: verts_j)
                {
                    arma::uword viu = static_cast <arma::uword> (
                            clexact_dat.vert2index_map.at (vi)),
                                vju = static_cast <arma::uword> (
                            clexact_dat.vert2index_map.at (vj));
                    if (vert_dists (viu, vju) > max_d)
                        max_d = vert_dists (viu, vju);
                }
            arma::uword iu = static_cast <arma::uword> (i),
                        ju = static_cast <arma::uword> (j);
            cl_edges (iu, ju) = cl_edges (ju, iu) = max_d;
        }
}


//' rcpp_exact_initial
//'
//' Initial allocation for exact clustering
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_exact_initial (
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

    exact_init::EXDat clexact_dat;
    exact_init::clexact_init (clexact_dat, from, to, d);

    exact_init::assign_first_edge (clexact_dat);
    int clnum = 1; // #1 assigned in assign_first_edge
    index_t ei = 1; // index of next edge to be assigned

    while (clexact_dat.vert2cl_map.size () < clexact_dat.n)
    {
        int clnum_i = clexact_step (clexact_dat, ei, clnum);
        ei++;
        if (clnum_i == clnum)
            clnum++;
    }

    // Then construct the hierarchical relationships among clusters
    arma::uword cu = static_cast <arma::uword> (clnum);
    arma::Mat <double> cl_edges (cu, cu);
    exact_init::fill_cl_edges (clexact_dat, cl_edges, clnum);

    // Then construct vector mapping edges to cluster numbers
    std::vector <int> clvec (clexact_dat.n);
    for (auto ci: clexact_dat.vert2cl_map)
        clvec [static_cast <size_t> (ci.first)] = ci.second;
    
    return Rcpp::wrap (clvec);
}
