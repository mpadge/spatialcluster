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
    Rcpp::IntegerVector clnum = gr ["cluster"];
    Rcpp::IntegerVector clfrom = gr ["cl_from"];
    Rcpp::IntegerVector clto = gr ["cl_to"];

    const size_t n = static_cast <size_t> (d.size ());
    
    cldat.edges.resize (n);
    int2intset_map_t cl2edge_map; // TODO: Delete that!
    std::unordered_map <int, std::unordered_set <double> > cl2dist_map;
    std::unordered_map <std::string, double> edge_dist_map;
    for (int i = 1; i < static_cast <int> (n); i++)
    {
        if (clnum [i] >= 0) // edge in a cluster
        {
            int clnum_i = clnum [i];
            intset_t edgeset;
            std::unordered_set <double> distset;
            if (cl2edge_map.find (clnum_i) != cl2edge_map.end ())
            {
                edgeset = cl2edge_map.at (clnum_i);
                distset = cl2dist_map.at (clnum_i);
            }
            edgeset.emplace (i);
            distset.emplace (d [i]);
            cl2edge_map [clnum_i] = edgeset;
            cl2dist_map [clnum_i] = distset;
        } else 
        { // make set of unordered edge names
            std::string eft = std::to_string (clfrom [i]) + "-" +
                              std::to_string (clto [i]),
                        etf = std::to_string (clto [i]) + "-" +
                              std::to_string (clfrom [i]);
            if (edge_dist_map.find (eft) == edge_dist_map.end () &&
                    edge_dist_map.find (etf) == edge_dist_map.end ())
                edge_dist_map.emplace (eft, d [i]); // d[i] not used here
        }
    }

    cldat.edges.resize (edge_dist_map.size ());
    edge_dist_map.clear ();
    size_t edge_count = 0;
    for (int i = 1; i < static_cast <int> (n); i++)
    {
        if (clnum [i] < 0) // edge not in a cluster
        {
            utils::OneEdge edgei;
            // from and to hold cluster numbers, NOT vertex numbers
            edgei.from = clfrom [i];
            edgei.to = clto [i];
            edgei.dist = d [i];

            std::string eft = std::to_string (edgei.from) + "-" +
                              std::to_string (edgei.to),
                        etf = std::to_string (edgei.to) + "-" +
                              std::to_string (edgei.from);
            if (edge_dist_map.find (eft) == edge_dist_map.end () &&
                    edge_dist_map.find (etf) == edge_dist_map.end ())
            {
                edge_dist_map.emplace (eft, d [i]);
                cldat.edges [edge_count++] = edgei;
            } else if (edge_dist_map.find (etf) != edge_dist_map.end ())
            {
                if ((cldat.shortest && d [i] < edge_dist_map.at (etf)) ||
                        (!cldat.shortest && d [i] > edge_dist_map.at (etf)))
                    edge_dist_map [etf] = d [i];
            } else
            {
                if ((cldat.shortest && d [i] < edge_dist_map.at (eft)) ||
                        (!cldat.shortest && d [i] > edge_dist_map.at (eft)))
                    edge_dist_map [eft] = d [i];
            }
        }
    }
    // Then just loop over cldat.edges to update min distances
    for (auto ei: cldat.edges)
    {
        std::string eft = std::to_string (ei.from) + "-" +
                          std::to_string (ei.to);
        if ((cldat.shortest && edge_dist_map.at (eft) < ei.dist) ||
                (!cldat.shortest && edge_dist_map.at (eft) > ei.dist))
            ei.dist = edge_dist_map.at (eft);
    }

    // Fill intra-cluster data:
    for (auto i: cl2dist_map)
    {
        //intset_t edgeset = i.second;
        std::unordered_set <double> distset = i.second;
        OneCluster cli;
        cli.id = i.first;
        cli.n = distset.size ();
        cli.dist_sum = 0.0;
        cli.dist_max = 0.0;
        if (!cldat.shortest)
            cli.dist_max = INFINITE_DOUBLE;
        for (auto di: distset)
        {
            cli.dist_sum += di;
            // TODO: Ensure that this is correct for !shortest
            if ((cldat.shortest && di > cli.dist_max) ||
                    (!cldat.shortest && di < cli.dist_max))
                cli.dist_max = di;
        }
        cldat.clusters.emplace (i.first, cli);

        cldat.cl_remap.emplace (i.first, i.first);
        intset_t members;
        members.emplace (i.first);
        cldat.cl_members.emplace (i.first, members);
    }
}

// merge cluster clfrom with clto; clfrom remains as it was but is no longer
// indexed so simply ignored from that point on
ex_merge::OneMerge ex_merge::merge_one_single (ex_merge::ExMergeDat &cldat,
        index_t ei)
{
    const int cl_from_i = cldat.cl_remap.at (cldat.edges [ei].from),
              cl_to_i = cldat.cl_remap.at (cldat.edges [ei].to);

    ex_merge::OneCluster clfrom = cldat.clusters.at (cl_from_i),
                         clto = cldat.clusters.at (cl_to_i);
    clto.n += clfrom.n;
    clto.dist_sum += clfrom.dist_sum;
    // TODO: Same as TODO in init fn above
    if ((cldat.shortest && clfrom.dist_max > clto.dist_max) ||
            (!cldat.shortest && clfrom.dist_max < clto.dist_max))
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
void ex_merge::merge_single (ex_merge::ExMergeDat &cldat)
{
    index_t edgei = 0;
    while (cldat.clusters.size () > 1)
    {
        int clfr = cldat.cl_remap.at (cldat.edges [edgei].from),
            clto = cldat.cl_remap.at (cldat.edges [edgei].to);
        if (clfr != clto)
        {
            ex_merge::OneMerge the_merge =
                ex_merge::merge_one_single (cldat, edgei);
            cldat.merges.push_back (the_merge);
        }
        edgei++;
        if (edgei == cldat.edges.size ())
            break;
    }
}

bool ex_merge::avgdist_sorter_incr (const OneAvgDist &lhs,
        const OneAvgDist &rhs)
{
    return lhs.average < rhs.average;
}

bool ex_merge::avgdist_sorter_decr (const OneAvgDist &lhs,
        const OneAvgDist &rhs)
{
    return lhs.average > rhs.average;
}

bool ex_merge::maxdist_sorter_incr (const OneAvgDist &lhs,
        const OneAvgDist &rhs)
{
    return lhs.d < rhs.d;
}

bool ex_merge::maxdist_sorter_decr (const OneAvgDist &lhs,
        const OneAvgDist &rhs)
{
    return lhs.d > rhs.d;
}

void ex_merge::fill_avg_dists (ex_merge::ExMergeDat &cldat,
        ex_merge::AvgDists &cl_dists)
{
    cl_dists.avg_dists.resize (cldat.edges.size ());
    size_t nc = 0;
    std::unordered_set <std::string> edgenames; // TODO: Remove
    for (auto ei: cldat.edges)
    {
        ex_merge::OneAvgDist onedist;
        onedist.cli = ei.from;
        onedist.clj = ei.to;
        onedist.d = ei.dist;
        onedist.di = cldat.clusters [ei.from].dist_sum;
        onedist.dj = cldat.clusters [ei.to].dist_sum;
        onedist.ni = cldat.clusters [ei.from].n;
        onedist.nj = cldat.clusters [ei.to].n;

        onedist.average = (onedist.di + onedist.dj + onedist.d) /
            static_cast <double> (onedist.ni + onedist.nj + 1);

        cl_dists.avg_dists [nc++] = onedist;
    }

    if (cldat.shortest)
        std::sort (cl_dists.avg_dists.begin (), cl_dists.avg_dists.end (),
                &ex_merge::avgdist_sorter_incr);
    else
        std::sort (cl_dists.avg_dists.begin (), cl_dists.avg_dists.end (),
                &ex_merge::avgdist_sorter_decr);
}

// Fill the cli_map and clj_map entries which map cluster numbers onto sets of
// indices in cl_dists.avg_dists
void ex_merge::fill_cl_indx_maps (ex_merge::AvgDists &cl_dists)
{
    cl_dists.cl_map.clear ();
    for (size_t i = 0; i < cl_dists.avg_dists.size (); i++)
    {
        indxset_t indxs;
        const int cli = cl_dists.avg_dists [i].cli;
        if (cl_dists.cl_map.find (cli) != cl_dists.cl_map.end ())
            indxs = cl_dists.cl_map.at (cli);
        indxs.emplace (i);
        cl_dists.cl_map [cli] = indxs;

        indxs.clear ();
        const int clj = cl_dists.avg_dists [i].clj;
        if (cl_dists.cl_map.find (clj) != cl_dists.cl_map.end ())
            indxs = cl_dists.cl_map.at (clj);
        indxs.emplace (i);
        cl_dists.cl_map [clj] = indxs;
    }
}

// Merging is based on AvgDists, which holds all possible pair-wise merges of
// existing clusters. One merge combines the pair in one AvgDists item to make a
// new one. The convention is to merge cli into clj, so cli disappears.
// Importantly, this requires updating all other AvgDists.avg_dists items which
// contain either one of the newly merged pairs. Indices from clusters to
// AvgDists.avg_dists are kept in AvgDists.cli_map and .clj_map. The values of
// the latter are updated to reflect merges, as are the entries of the new cli
// in AvgDists.avg_dists.
ex_merge::OneMerge ex_merge::merge_avg (ex_merge::ExMergeDat &cldat,
        ex_merge::AvgDists &cl_dists)
{
    ex_merge::OneAvgDist the_dist = cl_dists.avg_dists [0];
    const double dtot = the_dist.di + the_dist.dj + the_dist.d;
    const size_t ntot = the_dist.ni + the_dist.nj + 1;
    const double average = dtot / static_cast <double> (ntot);
    const int cli = the_dist.cli,
              clj = the_dist.clj;
    double dmin = INFINITE_DOUBLE; // shortest connecting distance
    if (!cldat.shortest)
        dmin = -dmin;

    indxset_t cli_indx = cl_dists.cl_map.at (cli),
              clj_indx = cl_dists.cl_map.at (clj);
    // update cli_indx & clj_indx entries, and get value of dmin
    for (auto i: clj_indx)
    {
        cl_dists.avg_dists [i].dj = dtot;
        cl_dists.avg_dists [i].nj = ntot;
        if (cl_dists.avg_dists [i].cli == cli)
            cl_dists.avg_dists [i].cli = clj;
        else if (cl_dists.avg_dists [i].clj == cli)
            cl_dists.avg_dists [i].clj = clj;
        if ((cldat.shortest && cl_dists.avg_dists [i].d < dmin) ||
                (!cldat.shortest && cl_dists.avg_dists [i].d > dmin))
            dmin = cl_dists.avg_dists [i].d;
    }
    for (auto i: cli_indx)
    {
        cl_dists.avg_dists [i].di = dtot;
        cl_dists.avg_dists [i].ni = ntot;
        if (cl_dists.avg_dists [i].cli == cli)
            cl_dists.avg_dists [i].cli = clj;
        else if (cl_dists.avg_dists [i].clj == cli)
            cl_dists.avg_dists [i].clj = clj;
        if ((cldat.shortest && cl_dists.avg_dists [i].d < dmin) ||
                (!cldat.shortest && cl_dists.avg_dists [i].d > dmin))
            dmin = cl_dists.avg_dists [i].d;
    }
    // Then update all dmin and average dist values
    for (auto i: clj_indx)
    {
        cl_dists.avg_dists [i].d = dmin;
        cl_dists.avg_dists [i].average =
            (cl_dists.avg_dists [i].di + dtot + dmin) /
            static_cast <double> (cl_dists.avg_dists [i].ni + ntot + 1);
    }
    for (auto i: cli_indx)
    {
        cl_dists.avg_dists [i].d = dmin;
        cl_dists.avg_dists [i].average =
            (cl_dists.avg_dists [i].dj + dtot + dmin) /
            static_cast <double> (cl_dists.avg_dists [i].nj + ntot + 1);
    }
    cl_dists.avg_dists.pop_front ();

    // These can now have reverse-duplicated entries because after merging A->B
    // entries A->C and C->B will become B->C and C->B. There can also be D->A
    // and D->B which will both become D->B.
    std::vector <int> rm;
    std::unordered_set <std::string> edge_names;
    for (size_t i = 0; i < cl_dists.avg_dists.size (); i++)
    {
        std::string cij = std::to_string (cl_dists.avg_dists [i].cli) + "-" +
                          std::to_string (cl_dists.avg_dists [i].clj),
                    cji = std::to_string (cl_dists.avg_dists [i].clj) + "-" +
                          std::to_string (cl_dists.avg_dists [i].cli);
        if (edge_names.find (cij) == edge_names.end () &&
                edge_names.find (cji) == edge_names.end ())
            edge_names.emplace (cij);
        else
            rm.push_back (static_cast <int> (i));
    std::unordered_set <int> merged;
    }
    std::sort (rm.begin (), rm.end (), std::greater <int> ());
    for (auto i: rm)
        cl_dists.avg_dists.erase (cl_dists.avg_dists.begin () + i);

    if (cldat.shortest)
        std::sort (cl_dists.avg_dists.begin (), cl_dists.avg_dists.end (),
                &ex_merge::avgdist_sorter_incr);
    else
        std::sort (cl_dists.avg_dists.begin (), cl_dists.avg_dists.end (),
                &ex_merge::avgdist_sorter_decr);

    // Finally, update the cl_dists.cli_map & clj_map entries
    fill_cl_indx_maps (cl_dists);

    ex_merge::OneMerge the_merge;
    the_merge.cli = cli;
    the_merge.clj = clj;
    the_merge.merge_dist = average;

    return the_merge;
}

// Successively merge pairs of clusters which yield the lower average
// intra-cluster edge distance
void ex_merge::avg (ex_merge::ExMergeDat &cldat)
{
    AvgDists cl_dists;
    ex_merge::fill_avg_dists (cldat, cl_dists);
    ex_merge::fill_cl_indx_maps (cl_dists);

    while (cl_dists.avg_dists.size () > 1)
    {
        ex_merge::OneMerge the_merge = ex_merge::merge_avg (cldat, cl_dists);
        cldat.merges.push_back (the_merge);
    }
}

void ex_merge::fill_max_dists (ex_merge::ExMergeDat &cldat,
        ex_merge::AvgDists &cl_dists)
{
    cl_dists.avg_dists.resize (cldat.edges.size ());
    size_t nc = 0;
    std::unordered_set <std::string> edgenames; // TODO: Remove
    for (auto ei: cldat.edges)
    {
        ex_merge::OneAvgDist onedist;
        onedist.cli = ei.from;
        onedist.clj = ei.to;
        onedist.d = ei.dist;

        cl_dists.avg_dists [nc++] = onedist;
    }

    if (cldat.shortest)
        std::sort (cl_dists.avg_dists.begin (), cl_dists.avg_dists.end (),
                &ex_merge::maxdist_sorter_incr);
    else
        std::sort (cl_dists.avg_dists.begin (), cl_dists.avg_dists.end (),
                &ex_merge::maxdist_sorter_decr);
}

void ex_merge::max (ex_merge::ExMergeDat &cldat)
{
}


ex_merge::OneMerge ex_merge::merge_max (ex_merge::ExMergeDat &cldat,
        ex_merge::AvgDists &cl_dists)
{
    ex_merge::OneAvgDist the_dist = cl_dists.avg_dists [0];
    const double dtot = the_dist.di + the_dist.dj + the_dist.d;
    const size_t ntot = the_dist.ni + the_dist.nj + 1;
    const double average = dtot / static_cast <double> (ntot);
    const int cli = the_dist.cli,
              clj = the_dist.clj;
    double dmin = INFINITE_DOUBLE; // shortest connecting distance
    if (!cldat.shortest)
        dmin = -dmin;

    indxset_t cli_indx = cl_dists.cl_map.at (cli),
              clj_indx = cl_dists.cl_map.at (clj);
    // update cli_indx & clj_indx entries, and get value of dmin
    for (auto i: clj_indx)
    {
        cl_dists.avg_dists [i].dj = dtot;
        cl_dists.avg_dists [i].nj = ntot;
        if (cl_dists.avg_dists [i].cli == cli)
            cl_dists.avg_dists [i].cli = clj;
        else if (cl_dists.avg_dists [i].clj == cli)
            cl_dists.avg_dists [i].clj = clj;
        if ((cldat.shortest && cl_dists.avg_dists [i].d < dmin) ||
                (!cldat.shortest && cl_dists.avg_dists [i].d > dmin))
            dmin = cl_dists.avg_dists [i].d;
    }
    for (auto i: cli_indx)
    {
        cl_dists.avg_dists [i].di = dtot;
        cl_dists.avg_dists [i].ni = ntot;
        if (cl_dists.avg_dists [i].cli == cli)
            cl_dists.avg_dists [i].cli = clj;
        else if (cl_dists.avg_dists [i].clj == cli)
            cl_dists.avg_dists [i].clj = clj;
        if ((cldat.shortest && cl_dists.avg_dists [i].d < dmin) ||
                (!cldat.shortest && cl_dists.avg_dists [i].d > dmin))
            dmin = cl_dists.avg_dists [i].d;
    }
    // Then update all dmin and average dist values
    for (auto i: clj_indx)
    {
        cl_dists.avg_dists [i].d = dmin;
        cl_dists.avg_dists [i].average =
            (cl_dists.avg_dists [i].di + dtot + dmin) /
            static_cast <double> (cl_dists.avg_dists [i].ni + ntot + 1);
    }
    for (auto i: cli_indx)
    {
        cl_dists.avg_dists [i].d = dmin;
        cl_dists.avg_dists [i].average =
            (cl_dists.avg_dists [i].dj + dtot + dmin) /
            static_cast <double> (cl_dists.avg_dists [i].nj + ntot + 1);
    }
    cl_dists.avg_dists.pop_front ();

    // These can now have reverse-duplicated entries because after merging A->B
    // entries A->C and C->B will become B->C and C->B. There can also be D->A
    // and D->B which will both become D->B.
    std::vector <int> rm;
    std::unordered_set <std::string> edge_names;
    for (size_t i = 0; i < cl_dists.avg_dists.size (); i++)
    {
        std::string cij = std::to_string (cl_dists.avg_dists [i].cli) + "-" +
                          std::to_string (cl_dists.avg_dists [i].clj),
                    cji = std::to_string (cl_dists.avg_dists [i].clj) + "-" +
                          std::to_string (cl_dists.avg_dists [i].cli);
        if (edge_names.find (cij) == edge_names.end () &&
                edge_names.find (cji) == edge_names.end ())
            edge_names.emplace (cij);
        else
            rm.push_back (static_cast <int> (i));
    std::unordered_set <int> merged;
    }
    std::sort (rm.begin (), rm.end (), std::greater <int> ());
    for (auto i: rm)
        cl_dists.avg_dists.erase (cl_dists.avg_dists.begin () + i);

    if (cldat.shortest)
        std::sort (cl_dists.avg_dists.begin (), cl_dists.avg_dists.end (),
                &ex_merge::avgdist_sorter_incr);
    else
        std::sort (cl_dists.avg_dists.begin (), cl_dists.avg_dists.end (),
                &ex_merge::avgdist_sorter_decr);

    // Finally, update the cl_dists.cli_map & clj_map entries
    fill_cl_indx_maps (cl_dists);

    ex_merge::OneMerge the_merge;
    the_merge.cli = cli;
    the_merge.clj = clj;
    the_merge.merge_dist = average;


    return the_merge;
}

//' rcpp_exact_merge
//'
//' Merge clusters generated by rcpp_exact_initial to full hierarchy of all
//' possible merges.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_exact_merge (
        const Rcpp::DataFrame gr,
        const std::string linkage,
        const bool shortest)
{
    ex_merge::ExMergeDat clmerge_dat;
    clmerge_dat.shortest = shortest;
    ex_merge::init (gr, clmerge_dat);

    if (utils::strfound (linkage, "single"))
    {
        ex_merge::merge_single (clmerge_dat);
    } else if (utils::strfound (linkage, "average"))
    {
        ex_merge::avg (clmerge_dat);
    } else if (utils::strfound (linkage, "max"))
    {
        ex_merge::max (clmerge_dat);
    } else
        Rcpp::stop ("linkage not found for exact_merge");

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
