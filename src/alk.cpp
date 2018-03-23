#include "common.h"
#include "utils.h"
#include "bst.h"
#include "alk.h"

// --------- AVERAGE LINKAGE CLUSTER ----------------

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
    edge_tree->num_edges = arma::zeros <arma::Mat <unsigned short> > (n, n);
    edge_tree->avg_dist.set_size (n, n);
    edge_tree->avg_dist.fill (INFINITE_DOUBLE);
    edge_tree->dmat.set_size (n, n);
    edge_tree->dmat.fill (INFINITE_DOUBLE);
    for (int i = 0; i < from.length (); i++)
    {
        edge_tree->contig_mat (from [i], to [i]) = 1;
        edge_tree->num_edges (from [i], to [i]) = 1;
        edge_tree->dmat (from [i], to [i]) = d [i];
        edge_tree->avg_dist (from [i], to [i]) = 0.0;
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
    unsigned int e = 0;
    // T used to step through successive min values:
    Tree <double> * T = treeMinTree (edge_tree->tree);
    while (cfrom == cto ||
            edge_tree->contig_mat (vfrom, vto) == 0 ||
            d [edge_i] < min_cl_dist)
    {
        e++;
        edge_dist = T->data;
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

    // Then merge clusters and update inter-cluster avg_dists
    std::unordered_set <double> removed_dists;
    for (auto cl: edge_tree->cl2vert_map)
    {
        if (cl.first != cfrom && cl.first != cto)
        {
            edge_tree->avg_dist (cl.first, cfrom) =
                (edge_tree->avg_dist (cl.first, cfrom) *
                 edge_tree->num_edges (cl.first, cfrom) +
                 edge_tree->avg_dist (cl.first, cto) *
                 edge_tree->num_edges (cl.first, cto)) /
                (edge_tree->num_edges (cl.first, cfrom) +
                 edge_tree->num_edges (cl.first, cto));
            edge_tree->num_edges (cl.first, cfrom) =
                edge_tree->num_edges (cl.first, cfrom) + 
                edge_tree->num_edges (cl.first, cto);
            edge_tree->num_edges (cl.first, cto) =
                edge_tree->num_edges (cl.first, cfrom);

            if (edge_tree->contig_mat (cl.first, cfrom) == 1 ||
                    edge_tree->contig_mat (cl.first, cto) == 1)
            {
                edge_tree->contig_mat (cl.first, cfrom) = 1;
                edge_tree->contig_mat (cl.first, cto) = 1;

                std::set <unsigned int> verts_from = 
                    edge_tree->cl2vert_map.at (cl.first),
                    verts_to = edge_tree->cl2vert_map.at (cto);
                for (auto vi: verts_from)
                    for (auto vj: verts_to)
                    {
                        double tempd = edge_tree->avg_dist (vi, vj);
                        // TODO: tempd should never be INFINITE_DOUBLE here -
                        // check out why that happens and FIX!
                        if (tempd < INFINITE_DOUBLE &&
                                removed_dists.find (tempd) == removed_dists.end ())
                        {
                            removed_dists.emplace (tempd);
                            Rcpp::Rcout << "removing " << tempd << std::endl;
                            //treeDeleteNode (edge_tree->tree, tempd);
                        } else if (removed_dists.find (tempd) == removed_dists.end ())
                            Rcpp::Rcout << "---" << std::endl;
                    }
                verts_to = edge_tree->cl2vert_map.at (cfrom);
                for (auto vi: verts_from)
                    for (auto vj: verts_to)
                    {
                        double tempd = edge_tree->avg_dist (vi, vj);
                        if (tempd < INFINITE_DOUBLE &&
                                removed_dists.find (tempd) == removed_dists.end ())
                        {
                            removed_dists.emplace (tempd);
                            Rcpp::Rcout << "removing " << tempd << std::endl;
                            //treeDeleteNode (edge_tree->tree, tempd);
                        } else if (removed_dists.find (tempd) == removed_dists.end ())
                            Rcpp::Rcout << "---" << std::endl;
                    }
                // finally, add the new dist to the bst
                treeInsertNode (edge_tree->tree,
                        edge_tree->avg_dist (cl.first, cfrom));
            }
        }
    }
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

    std::unordered_set <unsigned int> the_tree;
    edge_tree_step (&edge_tree, from, to, d, the_tree);

    /*
    Rcpp::Rcout << "tree has [" << treeSize (edge_tree.tree) << "] nodes " <<
        "; and maps have [" << edge_tree.edgewt2id_map.size () << ", " <<
        edge_tree.id2edgewt_map.size () << "] entries" << std::endl;
    Rcpp::Rcout << "tree min = " << treeMin (edge_tree.tree) << std::endl;

    Tree <double> * T = treeSuccesorInOrder (edge_tree.tree);
    Rcpp::Rcout << "followed by ";
    for (int i = 0; i < 5; i++)
    {
        Rcpp::Rcout << T->data << " -> ";
        T = treeSuccesorInOrder (T);
    }
    Rcpp::Rcout << T->data << std::endl;
    */

    treeClear (edge_tree.tree);

    Rcpp::IntegerVector res;
    return res;
}
