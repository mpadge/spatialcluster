#include "common.h"
#include "cuttree.h"

double simple_variance (const std::vector <double> &x)
{
    double s2 = 0.0, s = 0.0;
    for (auto i: x)
    {
        s += i;
        s2 += i * i;
    }
    const double n = static_cast <double> (x.size ());
    return (s2 - s * s / n) / (n - 1.0);
}

// Internal variance of specified cluster number
double calc_variance (const std::vector <EdgeComponent> &edges,
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
    return (s2 - s * s / count) / (count - 1.0);
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

double sum_component_variances (const std::vector <EdgeComponent> &edges,
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
    sa2 = (sa2 - sa * sa / na) / (na - 1.0);
    sb2 = (sb2 - sb * sb / nb) / (nb - 1.0);
    return sa2 + sb2;
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
    double vmin = INFINITE_DOUBLE;
    int imin = INFINITE_INT, size_at_min = INFINITE_INT;
    std::unordered_set <int> tree_out;
    for (int i = 0; i < n; i++)
    {
        edges_copy.resize (0);
        edges_copy.shrink_to_fit ();
        edges_copy.resize (n);
        std::copy (cluster_edges.begin (), cluster_edges.end (),
                edges_copy.begin ());
        edges_copy.erase (edges_copy.begin () + i);
        std::unordered_set <int> tree = build_one_tree (edges_copy);
        double v = sum_component_variances (edges_copy, tree);
        if (v < vmin)
        {
            vmin = v;
            imin = i;
            size_at_min = tree.size ();
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
            if (count == imin)
                i.cluster_num = INFINITE_INT;
            else if (tree_out.find (i.from) == tree_out.end ())
                i.cluster_num = total_clusters + 1;
        }
        count++;
    }

    BestCut the_cut;
    the_cut.variance = vmin;
    the_cut.pos = imin;
    the_cut.n1 = size_at_min;
    the_cut.n2 = n - size_at_min - 1;

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

    const double var_full = calc_variance (edges, 0);

    std::vector <double> cluster_variances;
    cluster_variances.push_back (var_full);
    double var_tot = var_full;
    BestCut the_cut = find_min_cut (edges, 0);

    int num_clusters = 1;
    while (num_clusters < ncl)
    {
        size_t mini = *std::min_element (cluster_variances.begin (),
                cluster_variances.end ()); // cluster number to be split
        the_cut = find_min_cut (edges, mini);
        cluster_variances.insert (cluster_variances.begin () + mini,
                the_cut.variance);

        num_clusters++;
    }

    Rcpp::IntegerVector res;
    return res;
}
