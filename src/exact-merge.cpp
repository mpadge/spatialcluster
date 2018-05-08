#include "common.h"
#include "exact-merge.h"

// load data from rcpp_exact_initial into the ExMergeDat struct. The gr data are
// pre-sorted by increasing d.
void ex_merge::init (const Rcpp::DataFrame &gr,
        ex_merge::ExMergeDat &cldat)
{
    Rcpp::IntegerVector from = gr ["from"];
    Rcpp::IntegerVector to = gr ["to"];
    Rcpp::NumericVector d = gr ["d"];
    Rcpp::IntegerVector clnum = gr ["cl"];
    Rcpp::IntegerVector clfrom = gr ["cl_from"];
    Rcpp::IntegerVector clto = gr ["cl_to"];

    const size_t n = static_cast <size_t> (d.size ());
    
    cldat.edges.resize (n);
    int2intset_map_t cl2edge_map;
    size_t edge_count = 0;
    for (int i = 1; i < static_cast <int> (n); i++)
    {
        if (clnum [i] >= 0) // edge in a cluster
        {
            int clnum_i = clnum [i];
            intset_t edgeset;
            if (cl2edge_map.find (clnum_i) != cl2edge_map.end ())
                edgeset = cl2edge_map.at (clnum_i);
            edgeset.emplace (i);
            cl2edge_map [clnum_i] = edgeset;
        } else 
            edge_count++;
    }

    // fill inter-cluster edges
    cldat.edges.resize (edge_count);
    edge_count = 0;
    for (int i = 1; i < static_cast <int> (n); i++)
    {
        if (clnum [i] < 0) // edge not in a cluster
        {
            utils::OneEdge edgei;
            // from and to hold cluster numbers, NOT vertex numbers
            edgei.from = clfrom [i];
            edgei.to = clto [i];
            edgei.dist = d [i];
            cldat.edges [edge_count++] = edgei;
        }
    }

    // Fill intra-cluster data:
    for (int i = 1; i < static_cast <int> (n); i++)
    {
        OneCluster cli;
        intset_t edgeset = cl2edge_map.at (i);
        cli.id = i;
        cli.n = edgeset.size ();
        cli.dist_sum = 0.0;
        cli.dist_max = 0.0;
        for (auto ei: edgeset)
        {
            cli.dist_sum += ei;
            if (ei > cli.dist_max)
                cli.dist_max = ei;
        }
        cldat.clusters.emplace (i, cli);

        cldat.cl_remap.emplace (i, i);
        intset_t members;
        members.emplace (i);
        cldat.cl_members.emplace (i, members);
    }
}

// merge cluster clfrom with clto; clfrom remains as it was but is no longer
// indexed so simply ignored from that point on
ex_merge::OneMerge ex_merge::merge (ex_merge::ExMergeDat &cldat, index_t ei)
{
    const int cl_from_i = cldat.cl_remap.at (cldat.edges [ei].from),
              cl_to_i = cldat.cl_remap.at (cldat.edges [ei].to);

    ex_merge::OneCluster clfrom = cldat.clusters.at (cl_from_i),
                         clto = cldat.clusters.at (cl_to_i);
    clto.n += clfrom.n;
    clto.dist_sum += clfrom.dist_sum;
    if (clfrom.dist_max > clto.dist_max)
        clto.dist_max = clfrom.dist_max;

    std::vector <utils::OneEdge> edges_from = clfrom.edges,
                                 edges_to = clto.edges;
    edges_to.insert (edges_to.end (), edges_from.begin (), edges_from.end ());
    clto.edges.clear ();
    clto.edges.shrink_to_fit ();
    clto.edges = edges_to;

    cldat.clusters.erase (cl_from_i);
    cldat.clusters [cl_to_i] = clto;

    cldat.cl_remap [cl_from_i] = cldat.cl_remap [cl_to_i];
    intset_t members_f = cldat.cl_members.at (cl_from_i),
             members_t = cldat.cl_members.at (cl_to_i);
    members_t.insert (members_f.begin (), members_f.end ());
    cldat.cl_members [cl_to_i] = members_t;
    for (auto m: members_t)
        cldat.cl_remap [m] = cldat.cl_remap [cl_to_i];

    ex_merge::OneMerge the_merge;
    the_merge.cli = cl_from_i;
    the_merge.clj = cl_to_i;
    the_merge.merge_dist = cldat.edges [ei].dist;

    return the_merge;
}

// Each merge joins from to to; from remains unchanged but is no longer indexed.
// Edges nevertheless always refer to original (non-merged) cluster numbers, so
// need to be re-mapped via the cl_remap
void ex_merge::single (ex_merge::ExMergeDat &cldat)
{
    index_t edgei = 0;
    while (cldat.clusters.size () > 1)
    {
        int clfr = cldat.cl_remap.at (cldat.edges [edgei].from),
            clto = cldat.cl_remap.at (cldat.edges [edgei].to);
        if (clfr != clto)
        {
            OneMerge the_merge = merge (cldat, edgei);
            cldat.merges.push_back (the_merge);
        }
        edgei++;
        if (edgei == cldat.edges.size ())
            break;
    }
}

// Successively merge pairs of clusters which yield the lower average
// intra-cluster edge distance
void ex_merge::avg (ex_merge::ExMergeDat &cldat)
{
}

void ex_merge::max (ex_merge::ExMergeDat &cldat)
{
}



//' rcpp_exact_merge
//'
//' Merge clusters generated by rcpp_exact_initial to obtain specified number of
//' clusters.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_exact_merge (
        const Rcpp::DataFrame gr,
        const std::string method)
{
    ex_merge::ExMergeDat clmerge_dat;
    ex_merge::init (gr, clmerge_dat);

    if (utils::strfound (method, "single"))
    {
        ex_merge::single (clmerge_dat);
    } else if (utils::strfound (method, "average"))
    {
        ex_merge::avg (clmerge_dat);
    } else if (utils::strfound (method, "max"))
    {
        ex_merge::max (clmerge_dat);
    } else
        Rcpp::stop ("method not found for exact_merge");

    const size_t n = clmerge_dat.merges.size ();
    Rcpp::NumericMatrix res (static_cast <int> (n), 3);
    for (size_t i = 0; i < n; i++)
    {
        res (i, 0) = clmerge_dat.merges [i].cli;
        res (i, 1) = clmerge_dat.merges [i].clj;
        res (i, 2) = clmerge_dat.merges [i].merge_dist;
    }

    std::vector <std::string> colnames (3);
    colnames [0] = "from";
    colnames [1] = "to";
    colnames [2] = "dist";
    Rcpp::List dimnames (2);
    dimnames (1) = colnames;
    res.attr ("dimnames") = dimnames;

    return res;
}
