#include "common.h"
#include "utils.h"
#include "full-init.h"

// --------- FULL CLUSTER ----------------

void full_init::init (full_init::FullInitDat &clfull_dat,
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
    clfull_dat.n = vert_set.size ();

    index_t i = 0;
    for (auto v: vert_set)
    {
        clfull_dat.index2vert_map.emplace (i, v);
        clfull_dat.vert2index_map.emplace (v, i);
        clfull_dat.index2cl_map.emplace (i++, -1);
    }

    clfull_dat.edges.clear ();
    clfull_dat.edges.resize (static_cast <size_t> (from.size ()));
    for (int i = 0; i < from.size (); i++)
    {
        utils::OneEdge here;
        here.from = from [i];
        here.to = to [i];
        here.dist = d [i];
        clfull_dat.edges [static_cast <size_t> (i)] = here;
    }

    clfull_dat.index_in_cluster.resize (clfull_dat.n);
    std::fill (clfull_dat.index_in_cluster.begin (),
            clfull_dat.index_in_cluster.end (), false);
}


void full_init::assign_first_edge (full_init::FullInitDat &clfull_dat)
{
    int clnum = 0;
    index_t ei = 0;
    utils::OneEdge edge = clfull_dat.edges [ei];
    index_t ito = clfull_dat.vert2index_map.at (edge.to),
            ifrom = clfull_dat.vert2index_map.at (edge.from);

    clfull_dat.index2cl_map [ito] = clnum;
    clfull_dat.index2cl_map [ifrom] = clnum;
    clfull_dat.index_in_cluster [ito] = true;
    clfull_dat.index_in_cluster [ifrom] = true;

    intset_t cli;
    cli.insert (static_cast <int> (ito));
    cli.insert (static_cast <int> (ifrom));
    clfull_dat.cl2index_map.emplace (clnum, cli);

    clfull_dat.vert2cl_map.emplace (edge.to, clnum);
    clfull_dat.vert2cl_map.emplace (edge.from, clnum);
}


//' step
//'
//' All edges are initially in their own clusters. This merges edge#i with the
//' next closest edge
//'
//' @param ei The i'th edge of the sorted list of NN edge weights
//' @noRd
int full_init::step (full_init::FullInitDat &clfull_dat,
        const index_t ei, const int clnum)
{
    bool from_in = false, to_in = false;
    utils::OneEdge edge = clfull_dat.edges [ei];
    index_t ito = clfull_dat.vert2index_map.at (edge.to),
            ifrom = clfull_dat.vert2index_map.at (edge.from);
    if (clfull_dat.index_in_cluster [ito])
        to_in = true;
    if (clfull_dat.index_in_cluster [ifrom])
        from_in = true;

    int clnum_i = clnum;

    if (from_in && to_in)
    {
        clnum_i = INFINITE_INT;
    } else
    {
        if (from_in) // then to is not in cluster
            clnum_i = clfull_dat.index2cl_map [ifrom];
        else if (to_in)
            clnum_i = clfull_dat.index2cl_map [ito];

        clfull_dat.index_in_cluster [ito] = true;
        clfull_dat.index_in_cluster [ifrom] = true;
        clfull_dat.index2cl_map [ifrom] = clnum_i;
        clfull_dat.index2cl_map [ito] = clnum_i;

        intset_t cli;
        if (clfull_dat.cl2index_map.find (clnum_i) !=
                clfull_dat.cl2index_map.end ())
            cli = clfull_dat.cl2index_map.at (clnum_i);
        cli.insert (static_cast <int> (ito));
        cli.insert (static_cast <int> (ifrom));
        clfull_dat.cl2index_map [clnum_i] = cli;

        // These values may already be in the map here, but that's okay
        clfull_dat.vert2cl_map.emplace (edge.to, clnum_i);
        clfull_dat.vert2cl_map.emplace (edge.from, clnum_i);
    }

    return clnum_i;
}

//' fill_cl_edges
//'
//' Fill (arma) matrix of strongest/shortest connections between all clusters
//' used to construct the hierarchical relationships
//' @noRd
void full_init::fill_cl_edges (full_init::FullInitDat &clfull_dat,
        arma::Mat <double> &cl_edges, int num_clusters)
{
    int2intset_map_t vert_sets;
    for (int i = 0; i < num_clusters; i++)
    {
        intset_t verts;
        for (auto vi: clfull_dat.vert2cl_map)
        {
            if (vi.second == i)
            {
                verts.emplace (vi.first);
            }
        }
        vert_sets.emplace (i, verts);
    }

    // need a (sparse) matrix of all pairwise edge distances:
    arma::uword nu = static_cast <arma::uword> (clfull_dat.n);
    arma::Mat <double> vert_dists (nu, nu);
    if (!clfull_dat.shortest)
        vert_dists.fill (INFINITE_DOUBLE);
    for (auto ei: clfull_dat.edges)
    {
        arma::uword i = static_cast <arma::uword> (
                                    clfull_dat.vert2index_map.at (ei.from)),
                    j = static_cast <arma::uword> (
                                    clfull_dat.vert2index_map.at (ei.to));
        vert_dists (i, j) = vert_dists (j, i) = ei.dist;
    }

    for (int i = 0; i < (num_clusters - 1); i++)
        for (int j = (i + 1); j < num_clusters; j++)
        {
            intset_t verts_i = vert_sets.at (i),
                     verts_j = vert_sets.at (j);
            double max_d = 0.0;
            if (!clfull_dat.shortest)
                max_d = INFINITE_DOUBLE; // min covariance
            for (auto vi: verts_i)
                for (auto vj: verts_j)
                {
                    arma::uword viu = static_cast <arma::uword> (
                            clfull_dat.vert2index_map.at (vi)),
                                vju = static_cast <arma::uword> (
                            clfull_dat.vert2index_map.at (vj));
                    if ((clfull_dat.shortest &&
                                vert_dists (viu, vju) > max_d) ||
                        (!clfull_dat.shortest &&
                                vert_dists (viu, vju) < max_d))
                        max_d = vert_dists (viu, vju);
                }
            arma::uword iu = static_cast <arma::uword> (i),
                        ju = static_cast <arma::uword> (j);
            cl_edges (iu, ju) = cl_edges (ju, iu) = max_d;
        }
}


//' rcpp_full_initial
//'
//' Initial allocation for full clustering
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_full_initial (
        const Rcpp::DataFrame gr,
        bool shortest)
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

    full_init::FullInitDat clfull_dat;
    clfull_dat.shortest = shortest;
    full_init::init (clfull_dat, from, to, d);

    full_init::assign_first_edge (clfull_dat);
    int clnum = 1; // #1 assigned in assign_first_edge
    index_t ei = 1; // index of next edge to be assigned

    while (clfull_dat.vert2cl_map.size () < clfull_dat.n)
    {
        int clnum_i = full_init::step (clfull_dat, ei, clnum);
        ei++;
        if (clnum_i == clnum)
            clnum++;
    }

    // Then construct the hierarchical relationships among clusters
    arma::uword cu = static_cast <arma::uword> (clnum);
    arma::Mat <double> cl_edges (cu, cu);
    full_init::fill_cl_edges (clfull_dat, cl_edges, clnum);

    // Then construct vector mapping edges to cluster numbers
    std::vector <int> clvec (clfull_dat.n);
    for (auto ci: clfull_dat.vert2cl_map)
        clvec [static_cast <size_t> (ci.first)] = ci.second;
    
    return Rcpp::wrap (clvec);
}
