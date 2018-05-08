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
    int2intset_map_t cl2edge_map; // TODO: Delete that!
    std::unordered_map <int, std::unordered_set <double> > cl2dist_map;
    size_t edge_count = 0;
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
    for (auto i: cl2dist_map)
    {
        //intset_t edgeset = i.second;
        std::unordered_set <double> distset = i.second;
        OneCluster cli;
        cli.id = i.first;
        cli.n = distset.size ();
        cli.dist_sum = 0.0;
        cli.dist_max = 0.0;
        for (auto di: distset)
        {
            cli.dist_sum += di;
            if (di > cli.dist_max)
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

bool ex_merge::avgdist_sorter (const OneAvgDist &lhs, const OneAvgDist &rhs)
{
    return lhs.average < rhs.average;
}

size_t ex_merge::count_num_clusters (ex_merge::ExMergeDat &cldat,
        const std::unordered_map <std::string, double> &edge_dist_map)
{
    // Number of clusters is probably just cldat.edges.size (), but may be less
    // if multiple edges connect same clusters
    const size_t n = cldat.clusters.size ();
    size_t nc = 0;
    for (int i = 0; i < (n - 1); i++)
        for (int j = (i + 1); j < n; j++)
        {
            // The [] operator requires the unordered_map to be **NON**const
            std::string eij = std::to_string (cldat.clusters [i].id) + "-" +
                              std::to_string (cldat.clusters [j].id);
            if (edge_dist_map.find (eij) != edge_dist_map.end ())
                nc++;
        }

    return nc;
}

void ex_merge::fill_avg_dists (ex_merge::ExMergeDat &cldat,
        const std::unordered_map <std::string, double> &edge_dist_map,
        ex_merge::AvgDists &cl_dists)
{
    size_t nc = 0;
    for (int i = 0; i < (cldat.clusters.size () - 1); i++)
        for (int j = (i + 1); j < cldat.clusters.size (); j++)
        {
            std::string eij = std::to_string (cldat.clusters [i].id) + "-" +
                              std::to_string (cldat.clusters [j].id);
            if (edge_dist_map.find (eij) != edge_dist_map.end ())
            {
                ex_merge::OneAvgDist onedist;
                onedist.cli = cldat.clusters [i].id;
                onedist.clj = cldat.clusters [j].id;
                onedist.di = cldat.clusters [i].dist_sum;
                onedist.dj = cldat.clusters [j].dist_sum;
                onedist.ni = cldat.clusters [i].n;
                onedist.nj = cldat.clusters [j].n;
                onedist.d = edge_dist_map.at (eij);

                onedist.average = (onedist.di + onedist.dj + onedist.d) /
                    static_cast <double> (onedist.ni + onedist.nj + 1);

                cl_dists.avg_dists [nc++] = onedist;
            }
        }

    std::sort (cl_dists.avg_dists.begin (), cl_dists.avg_dists.end (),
            &ex_merge::avgdist_sorter);
}

// Successively merge pairs of clusters which yield the lower average
// intra-cluster edge distance
void ex_merge::avg (ex_merge::ExMergeDat &cldat)
{
    std::unordered_map <std::string, double> edge_dist_map;
    for (auto ei: cldat.edges)
    {
        std::string es = std::to_string (ei.from) + "-" +
                         std::to_string (ei.to);
        if (edge_dist_map.find (es) == edge_dist_map.end ())
            edge_dist_map.emplace (es, ei.dist);
        else if (ei.dist < edge_dist_map.at (es)) // this is not likely, but possible
            edge_dist_map [es] = ei.dist;
    }

    // work out how many unique inter-cluster edges there are. This is probably
    // just cldat.edges.size (), but may be less.
    size_t nc = ex_merge::count_num_clusters (cldat, edge_dist_map);

    AvgDists cl_dists;
    cl_dists.avg_dists.resize (nc);
    ex_merge::fill_avg_dists (cldat, edge_dist_map, cl_dists);

    for (auto i: cl_dists.avg_dists)
        Rcpp::Rcout << i.cli << " -> " << i.clj << ": " <<
            i.average << std::endl;

    Rcpp::Rcout << "---done---" << std::endl;
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
