#include "common.h"
#include "utils.h"
#include "bst.h"
#include "alk.h"

#include <math.h> // isnan

// --------- AVERAGE LINKAGE CLUSTER ----------------

/* TODO: To get this working, I need to properly trace the vert2cl and cl2vert
 * maps to index the matrices. The loops over vertices must be done the same as
 * in `utils/find_shortest_connection`, and the vertex sets updated as in
 * `slk/merge_clusters`
 */

void edge_tree_init (Edge_tree * edge_tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d)
{
    for (int i = 0; i < from.size (); i++)
    {
        if (i == 0)
            edge_tree->tree = treeNewNode (d [0]);
        else
            treeInsertNode (edge_tree->tree, d [i]);
        edge_tree->edgewt2id_map.emplace (d [i], i);
        edge_tree->id2edgewt_map.emplace (i, d [i]);
    }

    unsigned int n = get_n (from, to);
    edge_tree->n = n;
    edge_tree->contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    edge_tree->num_edges = arma::ones <arma::Mat <unsigned short> > (n, n);
    edge_tree->avg_dist.set_size (n, n);
    //edge_tree->avg_dist.fill (INFINITE_DOUBLE);
    edge_tree->avg_dist.fill (0.0);
    edge_tree->dmat.set_size (n, n);
    edge_tree->dmat.fill (INFINITE_DOUBLE);
    for (int i = 0; i < from.length (); i++)
    {
        edge_tree->contig_mat (from [i], to [i]) = 1;
        edge_tree->num_edges (from [i], to [i]) = 1;
        edge_tree->dmat (from [i], to [i]) = d [i];
        //edge_tree->avg_dist (from [i], to [i]) = 0.0;
        edge_tree->avg_dist (from [i], to [i]) = d [i];
    }

    sets_init (from, to, edge_tree->vert2cl_map, edge_tree->cl2vert_map);
}

void edge_tree_step (Edge_tree * edge_tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d,
        std::unordered_set <unsigned int> &the_tree)
{
    double edge_dist = treeMin (edge_tree->tree);
    unsigned int edge_i = edge_tree->edgewt2id_map.at (edge_dist);
    unsigned int vfrom = from [edge_i], vto = to [edge_i];

    int cfrom = edge_tree->vert2cl_map.at (vfrom),
        cto = edge_tree->vert2cl_map.at (vto);
    double min_cl_dist = edge_tree->avg_dist (cfrom, cto);
    if (edge_tree->avg_dist (cto, cfrom) < min_cl_dist)
        min_cl_dist = edge_tree->avg_dist (cto, cfrom);

    // Step through to find the minimal-distance edge that (i) connects
    // different clusters, (ii) represents contiguous clusters, and (iii) has
    // distance greater than the average dist between those 2 clusters.
    // T used to step through successive min values:
    Tree <double> * T = treeMinTree (edge_tree->tree);
    while (cfrom == cto ||
            edge_tree->contig_mat (vfrom, vto) == 0 ||
            d [edge_i] < min_cl_dist)
    {
        edge_dist = T->data;
        if (edge_tree->edgewt2id_map.find (edge_dist) == 
                edge_tree->edgewt2id_map.end ())
            Rcpp::stop ("shite, that shouldn't happen");
        edge_i = edge_tree->edgewt2id_map.at (edge_dist);
        vfrom = from [edge_i];
        vto = to [edge_i];
        cfrom = edge_tree->vert2cl_map.at (vfrom);
        cto = edge_tree->vert2cl_map.at (vto);
        min_cl_dist = edge_tree->avg_dist (cfrom, cto);
        if (edge_tree->avg_dist (cto, cfrom) < min_cl_dist)
            min_cl_dist = edge_tree->avg_dist (cto, cfrom);

        T = treeSuccesorInOrder (T); // pointer to node with next shortest d
    }

    int ishort = find_shortest_connection (from, to, d, edge_tree->dmat,
            edge_tree->cl2vert_map, cfrom, cto);
    the_tree.insert (ishort);
    merge_clusters (edge_tree->contig_mat, edge_tree->vert2cl_map,
            edge_tree->cl2vert_map, cfrom, cto);

    // Then update inter-cluster avg_dists, noting that cfrom is no longer part
    // of cl2vert_map, but remains a valid index into avg_dist and num_edges
    std::set <unsigned int> verts_to = edge_tree->cl2vert_map.at (cto);

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
        if (cl.first != cto)
        {
            std::set <unsigned int> verts_i =
                edge_tree->cl2vert_map.at (cl.first);

            const double tempd_f = edge_tree->avg_dist (cl.first, cfrom),
                   tempd_t = edge_tree->avg_dist (cl.first, cto);
            const unsigned int nedges_f = edge_tree->num_edges (cl.first, cfrom),
                  nedges_t = edge_tree->num_edges (cl.first, cto);

            edge_tree->avg_dist (cl.first, cfrom) =
                (tempd_f * nedges_f + tempd_t * nedges_t) /
                (nedges_f + nedges_t);
            edge_tree->num_edges (cl.first, cfrom) = nedges_f + nedges_t;

            if (edge_tree->contig_mat (cl.first, cfrom) == 1 ||
                    edge_tree->contig_mat (cfrom, cl.first) == 1 ||
                    edge_tree->contig_mat (cl.first, cto) == 1 ||
                    edge_tree->contig_mat (cto, cl.first) == 1)
            {
                edge_tree->contig_mat (cl.first, cfrom) = 1;
                edge_tree->contig_mat (cfrom, cl.first) = 1;
                edge_tree->contig_mat (cl.first, cto) = 1;
                edge_tree->contig_mat (cto, cl.first) = 1;
                
                //treeDeleteNode (edge_tree->tree, tempd_f);
                //treeDeleteNode (edge_tree->tree, tempd_t);
                Tree <double> * T2 = treeGetNode (edge_tree->tree, tempd_f);
                if (T2 != nullptr)
                {
                    double thisd = T2->data;
                    treeDeleteNode (edge_tree->tree, thisd);
                }

                T2 = treeGetNode (edge_tree->tree, tempd_t);
                if (T2 != nullptr)
                {
                    double thisd = T2->data;
                    treeDeleteNode (edge_tree->tree, thisd);
                }


                // finally, add the new dist to the bst
                treeInsertNode (edge_tree->tree,
                        edge_tree->avg_dist (cl.first, cfrom));
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
    const unsigned int n = get_n (from, to);

    Edge_tree edge_tree;
    edge_tree_init (&edge_tree, from, to, d);

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
