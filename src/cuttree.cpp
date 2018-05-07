#include "common.h"
#include "cuttree.h"

void fill_edges (std::vector <EdgeComponent> &edges,
    const std::vector <std::string> &from,
    const std::vector <std::string> &to,
    const Rcpp::NumericVector &d)
{
    std::unordered_map <std::string, int> node_map;
    int node_num = 0;
    for (auto f: from)
    {
        if (node_map.find (f) == node_map.end ())
            node_map.emplace (f, node_num++);
    }
    for (auto t: to)
    {
        if (node_map.find (t) == node_map.end ())
            node_map.emplace (t, node_num++);
    }

    for (size_t i = 0; i < edges.size (); i++)
    {
        EdgeComponent this_edge;
        this_edge.d = d [static_cast <int> (i)];
        this_edge.cluster_num = 0;
        this_edge.from = node_map.at (from [i]);
        this_edge.to = node_map.at (to [i]);
        edges [i] = this_edge;
    }
}

// Internal sum of squared deviations of specified cluster number (this is just
// the variance without the scaling by N)
double calc_ss (const std::vector <EdgeComponent> &edges,
        const int cluster_num)
{
    double s2 = 0.0, s = 0.0;
    double count = 0.0;
    for (auto i: edges)
        if (i.cluster_num == cluster_num)
        {
            s += i.d;
            s2 += i.d * i.d;
            count += 1.0;
        }
    //return (s2 - s * s / count) / (count - 1.0); // variance
    return (s2 - s * s / count);
}

size_t cluster_size (const std::vector <EdgeComponent> &edges,
        const int cluster_num)
{
    size_t n = 0;
    for (auto i: edges)
        if (i.cluster_num == cluster_num)
            n++;
    return n;
}

// Build connected component of tree starting from first edge.
std::unordered_set <int> build_one_tree (std::vector <EdgeComponent> &edges)
{
    std::unordered_set <int> tree;

    tree.emplace (edges [0].from);
    tree.emplace (edges [0].to);

    bool done = false;
    while (!done)
    {
        bool added = false;
        for (size_t j = 1; j < edges.size (); j++)
        {
            if (tree.find (edges [j].from) != tree.end () &&
                    tree.find (edges [j].to) == tree.end ())
            {
                tree.emplace (edges [j].to);
                added = true;
            } else if (tree.find (edges [j].to) != tree.end () &&
                    tree.find (edges [j].from) == tree.end ())
            {
                tree.emplace (edges [j].from);
                added = true;
            }
        }
        done = !added;
    }
    return tree;
}

TwoSS sum_component_ss (const std::vector <EdgeComponent> &edges,
        const std::unordered_set <int> &tree)
{
    double sa = 0.0, sa2 = 0.0, sb = 0.0, sb2 = 0.0, na = 0.0, nb = 0.0;
    for (auto e: edges)
    {
        if (tree.find (e.from) != tree.end ())
        {
            sa += e.d;
            sa2 += e.d * e.d;
            na += 1.0;
        } else
        {
            sb += e.d;
            sb2 += e.d * e.d;
            nb += 1.0;
        }
    }
    TwoSS res;
    // res.ss1 = (sa2 - sa * sa / na) / (na - 1.0); // variance
    res.ss1 = (sa2 - sa * sa / na);
    res.ss2 = (sb2 - sb * sb / nb);
    res.n1 = static_cast <int> (na);
    res.n2 = static_cast <int> (nb);
    return res;
}

// Find the component split of edges in cluster_num which yields the lowest sum
// of internal variance.
BestCut find_min_cut (const std::vector <EdgeComponent> &edges,
        const int cluster_num)
{
    size_t n = cluster_size (edges, cluster_num);

    // fill component vector
    std::vector <EdgeComponent> cluster_edges;
    cluster_edges.reserve (n);
    for (auto e: edges)
        if (e.cluster_num == cluster_num)
            cluster_edges.push_back (e);

    std::vector <EdgeComponent> edges_copy;

    // Remove each edge in turn
    BestCut the_cut;
    the_cut.pos = the_cut.n1 = the_cut.n2 = INFINITE_INT;
    the_cut.ss1 = the_cut.ss2 = INFINITE_DOUBLE;
    the_cut.ss_diff = 0.0; // default, coz search is over max ss_diff
    double ssmin = INFINITE_DOUBLE;
    // TODO: Rewrite this to just erase and re-insert a single edge each time
    for (int i = 0; i < n; i++)
    {
        edges_copy.resize (0);
        edges_copy.shrink_to_fit ();
        edges_copy.resize (n);
        std::copy (cluster_edges.begin (), cluster_edges.end (),
                edges_copy.begin ());
        edges_copy.erase (edges_copy.begin () + i);
        std::unordered_set <int> tree = build_one_tree (edges_copy);
        // only include groups with >= MIN_CLUSTER_SIZE members
        if (tree.size () >= MIN_CLUSTER_SIZE &&
                tree.size () < (edges_copy.size () - MIN_CLUSTER_SIZE - 1))
        {
            TwoSS ss = sum_component_ss (edges_copy, tree);
            if ((ss.ss1 + ss.ss2) < ssmin)
            {
                ssmin = ss.ss1 + ss.ss2;
                the_cut.pos = i;
                the_cut.ss1 = ss.ss1;
                the_cut.ss2 = ss.ss2;

                the_cut.n1 = ss.n1;
                the_cut.n2 = ss.n2;

                the_cut.nodes.clear ();
                for (auto t: tree)
                    the_cut.nodes.emplace (t);
            }
        }
    }

    if (the_cut.ss1 != INFINITE_DOUBLE)
        the_cut.ss_diff = calc_ss (edges, cluster_num) - the_cut.ss1 - the_cut.ss2;

    return the_cut;
}

//' rcpp_cut_tree
//'
//' Cut tree into specified number of clusters by minimising internal cluster
//' variance.
//'
//' @param tree tree to be processed
//'
//' @return Vector of cluster IDs for each tree edge
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_cut_tree (const Rcpp::DataFrame tree, const int ncl)
{
    Rcpp::StringVector from_in = tree ["from"];
    Rcpp::StringVector to_in = tree ["to"];
    Rcpp::NumericVector dref = tree ["d"];

    std::vector <std::string> from =
        Rcpp::as <std::vector <std::string> > (from_in); // implicit clone
    std::vector <std::string> to =
        Rcpp::as <std::vector <std::string> > (to_in);

    std::vector <EdgeComponent> edges (static_cast <size_t> (dref.size ()));
    fill_edges (edges, from, to, dref);

    BestCut the_cut = find_min_cut (edges, 0);
    std::vector <double> ss_diff, ss1, ss2;
    ss_diff.push_back (the_cut.ss_diff); // ss0 - ss1 - ss2
    ss1.push_back (the_cut.ss1);
    ss2.push_back (the_cut.ss2);
    // map from index into ss vectors to actual cluster numbers
    std::unordered_map <size_t, int> cluster_map;
    cluster_map.emplace (0, 0);

    int num_clusters = 1;
    // This loop fills the three vectors (ss_diff, ss1, ss2), as well as the
    // cluster_map.
    while (num_clusters < ncl)
    {
        auto mp = std::max_element (ss_diff.begin (), ss_diff.end ());
        long int maxi_int = std::distance (ss_diff.begin (), mp);
        // could assert non-negative here, but no need
        size_t maxi = static_cast <size_t> (maxi_int);
        if (cluster_map.find (maxi) == cluster_map.end ())
            Rcpp::Rcout << "ss_diff has no max element in cluster_map" << std::endl;
        int clnum = cluster_map.at (maxi);
        // maxi is index of cluster to be split

        if (ss_diff [maxi] == 0.0) // no further cuts possible
        {
            break;
        }
        
        the_cut = find_min_cut (edges, clnum);
        // Break old clnum into 2:
        int count = 0;
        for (auto &e: edges)
        {
            if (e.cluster_num == clnum)
            {
                if (count == the_cut.pos)
                    e.cluster_num = INFINITE_INT;
                else if (the_cut.nodes.find (e.from) == the_cut.nodes.end ())
                    e.cluster_num = num_clusters;
            }
            count++;
        }
        // find new best cut of now reduced cluster
        the_cut = find_min_cut (edges, clnum);

        ss_diff [maxi] = the_cut.ss_diff;
        ss1 [maxi] = the_cut.ss1;
        ss2 [maxi] = the_cut.ss2;
        // and also of new cluster
        the_cut = find_min_cut (edges, num_clusters);

        ss_diff.push_back (the_cut.ss_diff);
        ss1.push_back (the_cut.ss1);
        ss2.push_back (the_cut.ss2);

        cluster_map.emplace (num_clusters, ss_diff.size () - 1);

        num_clusters++;
    }

    Rcpp::IntegerVector res (edges.size ());
    for (int i = 0; i < edges.size (); i++)
    {
        if (edges [static_cast <size_t> (i)].cluster_num == INFINITE_INT)
            res [i] = NA_INTEGER;
        else
            res [i] = edges [static_cast <size_t> (i)].cluster_num;
    }
    return res;
}
