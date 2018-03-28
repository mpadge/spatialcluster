#include "common.h"
#include "utils.h"
#include "bst.h"
#include "alk.h"

// --------- AVERAGE LINKAGE CLUSTER ----------------
//
// TODO: Remove all usage of the vert2cl and cl2vert pairs, and replace
// throughout with edgewt2clpair_map.

void edge_tree_init (Edge_tree * edge_tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d)
{
    unsigned int n = sets_init (from, to, edge_tree->vert2index_map,
            edge_tree->index2vert_map, edge_tree->vert2cl_map,
            edge_tree->cl2vert_map);
    edge_tree->n = n;

    std::unordered_set <unsigned int> vert_set;

    for (int i = 0; i < from.size (); i++)
    {
        vert_set.emplace (from [i]);
        vert_set.emplace (to [i]);
        if (i == 0)
            edge_tree->tree = treeNewNode (d [0]);
        else
            treeInsertNode (edge_tree->tree, d [i]);
    }
    unsigned int i = 0;
    for (auto v: vert_set)
        edge_tree->vert2index_map.emplace (v, i++);

    for (int i = 0; i < from.size (); i++)
    {
        edge_tree->edgewt2clpair_map.emplace (d [i],
                std::make_pair (edge_tree->vert2index_map.at (from [i]),
                    edge_tree->vert2index_map.at (to [i])));
    }

    edge_tree->contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    edge_tree->num_edges = arma::ones <arma::Mat <unsigned short> > (n, n);
    edge_tree->avg_dist.set_size (n, n);
    //edge_tree->avg_dist.fill (INFINITE_DOUBLE);
    edge_tree->avg_dist.fill (0.0);
    edge_tree->dmat.set_size (n, n);
    edge_tree->dmat.fill (INFINITE_DOUBLE);
    for (int i = 0; i < from.length (); i++)
    {
        unsigned int vf = edge_tree->vert2index_map.at (from [i]),
                     vt = edge_tree->vert2index_map.at (to [i]);
        edge_tree->contig_mat (vf, vt) = 1;
        edge_tree->num_edges (vf, vt) = 1;
        //edge_tree->avg_dist (vf, vt) = 0.0;
        edge_tree->avg_dist (vf, vt) = d [i];
        edge_tree->dmat (vf, vt) = d [i];
    }
}

void edge_tree_step (Edge_tree * edge_tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d,
        std::unordered_set <unsigned int> &the_tree)
{
    // Step through to find the minimal-distance edge that (i) connects
    // different clusters, (ii) represents contiguous clusters, and (iii) has
    // distance greater than the average dist between those 2 clusters.
    // T used to step through successive min values:
    Tree <double> * T = treeMinTree (edge_tree->tree);
    double edge_dist = T->data;
    std::pair <unsigned int, unsigned int> pr =
        edge_tree->edgewt2clpair_map.at (edge_dist);
    unsigned int l = edge_tree->vert2index_map.at (pr.first),
                 m = edge_tree->vert2index_map.at (pr.second);
    while (l == m || edge_tree->contig_mat (l, m) == 0 ||
            edge_dist < edge_tree->avg_dist (l, m))
    {
        T = treeSuccessorInOrder (T);
        edge_dist = T->data;
        pr = edge_tree->edgewt2clpair_map.at (edge_dist);
        l = edge_tree->vert2index_map.at (pr.first);
        m = edge_tree->vert2index_map.at (pr.second);
    }

    int ishort = find_shortest_connection (from, to, d,
            edge_tree->vert2index_map, edge_tree->dmat,
            edge_tree->cl2vert_map, l, m);
    // ishort is an index into (from, to)
    the_tree.insert (ishort);
    merge_clusters (edge_tree->contig_mat, edge_tree->vert2cl_map,
            edge_tree->cl2vert_map, l, m);

    /* Cluster numbers start off here the same as vertex numbers, and so are
     * initially simple indices into the vert-by-vert matrices (contig_mat,
     * dmat, avg_dist, num_edges). As cluster form, numbers merge to one of the
     * pre-existing ones, so are still indexed into these same matrices which do
     * not change size. Cluster merging simply means that previous rows and
     * columns of these matrices will no longer be indexed, and all new indices
     * are derived from constantly updated values of vert2cl_map and
     * cl2vert_map.
     */
    for (auto cl: edge_tree->cl2vert_map)
    {
        if (cl.first != l || cl.first != m)
        {
            std::set <unsigned int> verts_i =
                edge_tree->cl2vert_map.at (cl.first);

            const double tempd_l = edge_tree->avg_dist (cl.first, l),
                   tempd_m = edge_tree->avg_dist (cl.first, m);
            const unsigned int nedges_l = edge_tree->num_edges (cl.first, l),
                  nedges_m = edge_tree->num_edges (cl.first, m);

            edge_tree->avg_dist (cl.first, l) =
                (tempd_l * nedges_l + tempd_m * nedges_m) /
                (nedges_l + nedges_m);
            edge_tree->num_edges (cl.first, l) = nedges_l + nedges_m;

            if (edge_tree->contig_mat (cl.first, l) == 1 ||
                    edge_tree->contig_mat (cl.first, m) == 1)
            {
                edge_tree->contig_mat (cl.first, l) = 1;
                if (tempd_l > 0.0)
                    treeDeleteNode (edge_tree->tree, tempd_l);
                if (tempd_m > 0.0)
                    treeDeleteNode (edge_tree->tree, tempd_m);
                
                treeInsertNode (edge_tree->tree,
                        edge_tree->avg_dist (cl.first, l));
            } // end if C(c, l) = 1 or C(c, m) = 1 in Guo's terminology
        } // end if cl.first != (cfrom, cto)
    } // end for over cl
}


//' rcpp_alk
//'
//' Full-order average linkage cluster redcap algorithm
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_alk (
        const Rcpp::DataFrame gr)
{
    Rcpp::IntegerVector from_ref = gr ["from"];
    Rcpp::IntegerVector to_ref = gr ["to"];
    Rcpp::NumericVector d = gr ["d"];
    // Rcpp classes are always passed by reference, so cloning is necessary to
    // avoid modifying the original data.frames.
    Rcpp::IntegerVector from = Rcpp::clone (from_ref);
    Rcpp::IntegerVector to = Rcpp::clone (to_ref);
    // Index vectors are 1-indexed, so
    from = from - 1;
    to = to - 1;

    Edge_tree edge_tree;
    edge_tree_init (&edge_tree, from, to, d);
    const unsigned int n = edge_tree.n;

    std::unordered_set <unsigned int> the_tree;
    while (the_tree.size () < (n - 1)) // tree has n - 1 edges
    {
        Rcpp::checkUserInterrupt ();
        edge_tree_step (&edge_tree, from, to, d, the_tree);
    }
    treeClear (edge_tree.tree);

    std::vector <int> treevec (the_tree.begin (), the_tree.end ());

    return Rcpp::wrap (treevec);
}
