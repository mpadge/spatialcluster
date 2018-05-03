#include "common.h"
#include "cuttree.h"

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

    int oldj = 1, j = 1;
    bool done = false;
    while (!done)
    {
        if (tree.find (edges [j].from) != tree.end () &&
                tree.find (edges [j].to) == tree.end ())
        {
            tree.emplace (edges [j].to);
            j++;
        } else if (tree.find (edges [j].to) != tree.end () &&
                tree.find (edges [j].from) == tree.end ())
        {
            tree.emplace (edges [j].from);
            j++;
        }
        if (j == oldj)
            done = true;
        else if (j >= edges.size ())
            j = 1;
        oldj = j;
    }
    return tree;
}

TwoSS sum_component_ss (const std::vector <EdgeComponent> &edges,
        const std::unordered_set <int> &tree)
{
    double sa = 0.0, sa2 = 0.0, sb = 0.0, sb2 = 0.0, na = 0.0, nb = 0.0;
    for (auto e:edges)
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
    return res;
}

// Find the component split of edges in cluster_num which yields the lowest sum
// of internal variance. NOTE: This modifies edges!
BestCut find_min_cut (std::vector <EdgeComponent> &edges,
        const int cluster_num)
{
    size_t n = cluster_size (edges, cluster_num);

    // fill component vector
    std::vector <EdgeComponent> cluster_edges;
    cluster_edges.reserve (n);
    std::unordered_set <int> cluster_numbers_set;
    for (auto i: edges)
    {
        if (i.cluster_num == cluster_num)
        {
            cluster_edges.push_back (i);
        }
        cluster_numbers_set.emplace (i.cluster_num);
    }
    const int total_clusters = static_cast <int> (cluster_numbers_set.size ());

    std::vector <EdgeComponent> edges_copy;

    // Remove each edge in turn
    std::unordered_set <int> tree_out;
    BestCut the_cut;
    the_cut.pos = INFINITE_INT;
    the_cut.ss1 = the_cut.ss2 = INFINITE_DOUBLE;
    double ssmin = INFINITE_DOUBLE;
    for (int i = 0; i < n; i++)
    {
        edges_copy.resize (0);
        edges_copy.shrink_to_fit ();
        edges_copy.resize (n);
        std::copy (cluster_edges.begin (), cluster_edges.end (),
                edges_copy.begin ());
        edges_copy.erase (edges_copy.begin () + i);
        std::unordered_set <int> tree = build_one_tree (edges_copy);
        TwoSS ss = sum_component_ss (edges_copy, tree);
        if ((ss.ss1 + ss.ss2) < ssmin)
        {
            ssmin = ss.ss1 + ss.ss2;
            the_cut.pos = i;
            the_cut.ss1 = ss.ss1;
            the_cut.ss2 = ss.ss2;

            tree_out.clear ();
            for (auto t: tree)
                tree_out.emplace (t);
        }
    }

    // Break old cluster_num into 2:
    int count = 0;
    for (auto i: edges)
    {
        if (i.cluster_num == cluster_num)
        {
            if (count == the_cut.pos)
                i.cluster_num = INFINITE_INT;
            else if (tree_out.find (i.from) == tree_out.end ())
                i.cluster_num = total_clusters + 1;
        }
        count++;
    }

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

    std::vector <EdgeComponent> edges (dref.size ());
    for (size_t i = 0; i < dref.size (); i++)
    {
        EdgeComponent this_edge;
        this_edge.d = dref [i];
        this_edge.cluster_num = 0;
        this_edge.from = node_map.at (from [i]);
        this_edge.to = node_map.at (to [i]);
        edges [i] = this_edge;
    }

    const double var_full = calc_ss (edges, 0);

    std::vector <double> cluster_ss;
    cluster_ss.push_back (var_full);
    double var_tot = var_full;
    BestCut the_cut = find_min_cut (edges, 0);

    int num_clusters = 1;
    while (num_clusters < ncl)
    {
        //size_t mini = *std::min_element (cluster_ss.begin (),
        //        cluster_ss.end ()); // cluster number to be split
        //the_cut = find_min_cut (edges, mini);
        //cluster_ss.insert (cluster_ss.begin () + mini, the_cut.ss1 + the_cut.ss2);

        num_clusters++;
    }

    Rcpp::IntegerVector res;
    return res;
}
