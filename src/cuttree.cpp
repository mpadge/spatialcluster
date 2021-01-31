#include "common.h"
#include "cuttree.h"

void cuttree::fill_edges (cuttree::TreeDat &tree,
        const std::vector <int> &from,
        const std::vector <int> &to,
        Rcpp::NumericVector &d)
{
    std::unordered_map <int, int> vert2index_map;
    intset_t vert_set;
    for (size_t i = 0; i < from.size (); i++)
    {
        vert_set.emplace (from [i]);
        vert_set.emplace (to [i]);
    }

    int vert_num = 0;
    for (auto v: vert_set)
    {
        vert2index_map.emplace (v, vert_num++);
    }

    for (size_t i = 0; i < tree.edges.size (); i++)
    {
        cuttree::EdgeComponent this_edge;
        this_edge.d = d [static_cast <int> (i)];
        this_edge.cluster_num = 0;
        this_edge.from = vert2index_map.at (from [i]);
        this_edge.to = vert2index_map.at (to [i]);
        tree.edges [i] = this_edge;
    }
}

// Internal sum of squared deviations of specified cluster number (this is just
// the variance without the scaling by N)
double cuttree::calc_ss (const std::vector <cuttree::EdgeComponent> &edges,
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

// Internal mean covariance of specified cluster number 
double cuttree::calc_covsum (const std::vector <cuttree::EdgeComponent> &edges,
        const int cluster_num)
{
    double s = 0.0;
    double count = 0.0;
    for (auto i: edges)
        if (i.cluster_num == cluster_num)
        {
            s += i.d;
            count += 1.0;
        }
    return s / count;
}

size_t cuttree::cluster_size (const std::vector <cuttree::EdgeComponent> &edges,
        const int cluster_num)
{
    size_t n = 0;
    for (auto i: edges)
        if (i.cluster_num == cluster_num)
            n++;
    return n;
}

// Build connected component of tree starting from first edge.
std::unordered_set <int> cuttree::build_one_tree (
        std::vector <cuttree::EdgeComponent> &edges)
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

cuttree::TwoSS cuttree::sum_component_ss (
        const std::vector <cuttree::EdgeComponent> &edges,
        const std::unordered_set <int> &tree_edges,
        const bool shortest)
{
    double sa = 0.0, sa2 = 0.0, sb = 0.0, sb2 = 0.0, na = 0.0, nb = 0.0;
    for (auto e: edges)
    {
        if (tree_edges.find (e.from) != tree_edges.end ())
        {
            sa += e.d;
            if (shortest)
                sa2 += e.d * e.d;
            na += 1.0;
        } else
        {
            sb += e.d;
            if (shortest)
                sb2 += e.d * e.d;
            nb += 1.0;
        }
    }

    cuttree::TwoSS res;
    // res.ss1 = (sa2 - sa * sa / na) / (na - 1.0); // variance
    if (shortest)
    {
        res.ss1 = (sa2 - sa * sa / na);
        res.ss2 = (sb2 - sb * sb / nb);
    } else
    { // covariances are mean values, *NOT* sums like SS values
        res.ss1 = sa / na;
        res.ss2 = sb / nb;
    }
    res.n1 = static_cast <int> (na);
    res.n2 = static_cast <int> (nb);
    return res;
}

// Find the component split of edges in cluster_num which yields the lowest sum
// of internal variance.
cuttree::BestCut cuttree::find_min_cut (
        const TreeDat &tree,
        const int cluster_num,
        const bool shortest)
{
    size_t n = cuttree::cluster_size (tree.edges, cluster_num);

    // fill component vector
    std::vector <cuttree::EdgeComponent> cluster_edges;
    cluster_edges.reserve (n);
    for (auto e: tree.edges)
        if (e.cluster_num == cluster_num)
            cluster_edges.push_back (e);

    std::vector <cuttree::EdgeComponent> edges_copy;

    // Remove each edge in turn
    cuttree::BestCut the_cut;
    the_cut.pos = the_cut.n1 = the_cut.n2 = INFINITE_INT;
    the_cut.ss1 = the_cut.ss2 = INFINITE_DOUBLE;
    the_cut.ss_diff = 0.0; // default, coz search is over max ss_diff
    double ssmin = INFINITE_DOUBLE;
    // TODO: Rewrite this to just erase and re-insert a single edge each time
    for (int i = 0; i < static_cast <int> (n); i++)
    {
        edges_copy.resize (0);
        edges_copy.shrink_to_fit ();
        edges_copy.resize (n);
        std::copy (cluster_edges.begin (), cluster_edges.end (),
                edges_copy.begin ());
        edges_copy.erase (edges_copy.begin () + i);
        std::unordered_set <int> tree_edges = cuttree::build_one_tree (edges_copy);
        // only include groups with >= MIN_CLUSTER_SIZE members
        if (tree_edges.size () >= cuttree::MIN_CLUSTER_SIZE &&
                tree_edges.size () < (edges_copy.size () -
                                cuttree::MIN_CLUSTER_SIZE - 1))
        {
            cuttree::TwoSS ss;
            ss = cuttree::sum_component_ss (edges_copy, tree_edges, shortest);

            if ((ss.ss1 + ss.ss2) < ssmin) // applies to both distances & cov
            {
                ssmin = ss.ss1 + ss.ss2;
                the_cut.pos = i;
                the_cut.ss1 = ss.ss1;
                the_cut.ss2 = ss.ss2;

                the_cut.n1 = ss.n1;
                the_cut.n2 = ss.n2;

                the_cut.nodes.clear ();
                for (auto te: tree_edges)
                    the_cut.nodes.emplace (te);
            }
        }
    }

    if (the_cut.ss1 < INFINITE_DOUBLE)
        the_cut.ss_diff = cuttree::calc_ss (tree.edges, cluster_num) -
            the_cut.ss1 - the_cut.ss2;

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
Rcpp::IntegerVector rcpp_cut_tree (const Rcpp::DataFrame tree, const int ncl,
        const bool shortest)
{
    Rcpp::IntegerVector from_in = tree ["from"];
    Rcpp::IntegerVector to_in = tree ["to"];
    Rcpp::NumericVector dref = tree ["d"];

    std::vector <int> from = Rcpp::as <std::vector <int> > (from_in);
    std::vector <int> to = Rcpp::as <std::vector <int> > (to_in);

    cuttree::TreeDat tree_dat;
    tree_dat.edges.resize (static_cast <size_t> (dref.size ()));
    cuttree::fill_edges (tree_dat, from, to, dref);

    cuttree::BestCut the_cut = cuttree::find_min_cut (tree_dat, 0, shortest);
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
            Rcpp::Rcout << "ss_diff has no max element in cluster_map" <<
                std::endl;
        int clnum = cluster_map.at (maxi);
        // maxi is index of cluster to be split

        if (ss_diff [maxi] == 0.0) // no further cuts possible
        {
            break;
        }
        
        the_cut = cuttree::find_min_cut (tree_dat, clnum, shortest);
        // Break old clnum into 2:
        int count = 0;
        for (auto &e: tree_dat.edges)
        {
            if (e.cluster_num == clnum)
            {
                if (count == the_cut.pos)
                    e.cluster_num = INFINITE_INT;
                else if (the_cut.nodes.find (e.from) == the_cut.nodes.end ())
                    e.cluster_num = num_clusters;
                count++;
            }
        }
        // find new best cut of now reduced cluster
        the_cut = cuttree::find_min_cut (tree_dat, clnum, shortest);

        ss_diff [maxi] = the_cut.ss_diff;
        ss1 [maxi] = the_cut.ss1;
        ss2 [maxi] = the_cut.ss2;
        // and also of new cluster
        the_cut = cuttree::find_min_cut (tree_dat, num_clusters, shortest);

        ss_diff.push_back (the_cut.ss_diff);
        ss1.push_back (the_cut.ss1);
        ss2.push_back (the_cut.ss2);

        cluster_map.emplace (num_clusters, ss_diff.size () - 1);

        num_clusters++;
    }

    Rcpp::IntegerVector res (tree_dat.edges.size ());
    for (int i = 0; i < static_cast <int> (tree_dat.edges.size ()); i++)
    {
        if (tree_dat.edges [static_cast <size_t> (i)].cluster_num == INFINITE_INT)
            res [i] = NA_INTEGER;
        else
            res [i] = tree_dat.edges [static_cast <size_t> (i)].cluster_num;
    }
    return res;
}
