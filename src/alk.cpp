#include "common.h"
#include "utils.h"
#include "alk.h"

// --------- AVERAGE LINKAGE CLUSTER ----------------

void alk_init (ALKDat &alk_dat,
        BinarySearchTree &tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d)
{
    unsigned int n = sets_init (from, to, alk_dat.vert2index_map,
            alk_dat.index2vert_map, alk_dat.index2cl_map,
            alk_dat.cl2index_map);
    alk_dat.n = n;

    std::unordered_set <unsigned int> vert_set;

    // Get set of unique vertices, and store binary tree of edge distances
    for (int i = 0; i < from.size (); i++)
    {
        vert_set.emplace (from [i]);
        vert_set.emplace (to [i]);
        tree.insert (d [i]);
    }
    // Construct vert2index_map to map each unique vertex to an index
    unsigned int i = 0;
    for (auto v: vert_set)
        alk_dat.vert2index_map.emplace (v, i++);

    // Construct idx2edgewt_map and edgewt2idx_pair_map
    for (int i = 0; i < from.size (); i++)
    {
        unsigned int fi = alk_dat.vert2index_map.at (from [i]),
                     ti = alk_dat.vert2index_map.at (to [i]);
        alk_dat.edgewt2idx_pair_map.emplace (d [i], std::make_pair (fi, ti));

        // idx2edgewt_map.second is an unordered set that needs to be expanded
        std::unordered_set <double> wtset;
        if (alk_dat.idx2edgewt_map.find (fi) ==
                alk_dat.idx2edgewt_map.end ())
        {
            wtset.emplace (d [i]);
            alk_dat.idx2edgewt_map.emplace (fi, wtset);
        } else
        {
            wtset = alk_dat.idx2edgewt_map.at (fi);
            wtset.emplace (d [i]);
            alk_dat.idx2edgewt_map [fi] = wtset;
        }

        // repeat same for "to" vertex
        if (alk_dat.idx2edgewt_map.find (ti) ==
                alk_dat.idx2edgewt_map.end ())
        {
            wtset.emplace (d [i]);
            alk_dat.idx2edgewt_map.emplace (ti, wtset);
        } else
        {
            wtset = alk_dat.idx2edgewt_map.at (ti);
            wtset.emplace (d [i]);
            alk_dat.idx2edgewt_map [ti] = wtset;
        }
    }

    alk_dat.contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    alk_dat.num_edges = arma::ones <arma::Mat <unsigned short> > (n, n);
    alk_dat.avg_dist.set_size (n, n);
    //alk_dat.avg_dist.fill (INFINITE_DOUBLE);
    alk_dat.avg_dist.fill (0.0);
    alk_dat.dmat.set_size (n, n);
    alk_dat.dmat.fill (INFINITE_DOUBLE);
    for (int i = 0; i < from.length (); i++)
    {
        unsigned int vf = alk_dat.vert2index_map.at (from [i]),
                     vt = alk_dat.vert2index_map.at (to [i]);
        alk_dat.contig_mat (vf, vt) = 1;
        alk_dat.num_edges (vf, vt) = 1;
        //alk_dat.avg_dist (vf, vt) = 0.0;
        alk_dat.avg_dist (vf, vt) = d [i];
        alk_dat.dmat (vf, vt) = d [i];
    }
}

// update both idx2edgewt and edgewt2idx maps to reflect merging of cluster m
// into cluster l (using Guo's original notation there). The cl2index
// and index2cl maps are updated in `merge_clusters`
void update_edgewt_maps (ALKDat &alk_dat,
        unsigned int m, unsigned int l)
{
    std::unordered_set <double> wtsl = alk_dat.idx2edgewt_map.at (l),
        wtsm = alk_dat.idx2edgewt_map.at (m);
    for (auto w: wtsm)
        wtsl.insert (w);
    alk_dat.idx2edgewt_map.erase (m);
    alk_dat.idx2edgewt_map.erase (l);
    alk_dat.idx2edgewt_map.emplace (l, wtsl);

    for (auto w: wtsl)
    {
        std::pair <unsigned int, unsigned int> pr =
            alk_dat.edgewt2idx_pair_map.at (w);
        bool update = true;
        if (pr.first == m)
            pr.first = l;
        else if (pr.second == m)
            pr.second = l;
        else
            update = false;
        if (update)
        {
            alk_dat.edgewt2idx_pair_map [w] = pr;
        }
    }

    // Any edgewt2idx pairs with entries of m also have to be re-mapped to l
    // TODO: Is there a better way to do this?
    std::unordered_set <double> wts;
    for (auto w: alk_dat.edgewt2idx_pair_map)
    {
        if (w.second.first == m || w.second.second == m)
            wts.emplace (w.first);
    }
    if (wts.size () > 0)
    {
        for (auto w: wts)
        {
            std::pair <unsigned int, unsigned int> pr =
                alk_dat.edgewt2idx_pair_map.at (w);
            if (pr.first == m)
                pr.first = l;
            if (pr.second == m)
                pr.second = l;
            alk_dat.edgewt2idx_pair_map [w] = pr;
        }
    }
}

int alk_step (ALKDat &alk_dat,
        BinarySearchTree &tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d)
{
    // Step through to find the minimal-distance edge that (i) connects
    // different clusters, (ii) represents contiguous clusters, and (iii) has
    // distance greater than the average dist between those 2 clusters.
    // T used to step through successive min values:
    double edge_dist = tree.treeMin();
    tree_node * node = tree.getRoot ();
    node = tree.getNode (node, edge_dist);
    std::pair <unsigned int, unsigned int> pr =
        alk_dat.edgewt2idx_pair_map.at (edge_dist);
    unsigned int l = pr.first, m = pr.second;
    while (l == m || alk_dat.contig_mat (l, m) == 0 ||
            edge_dist < alk_dat.avg_dist (l, m))
    {
        node = tree.nextHi (node);
        if (node == nullptr)
            Rcpp::stop ("can not go past highest node");
        edge_dist = node->data;
        pr = alk_dat.edgewt2idx_pair_map.at (edge_dist);
        l = pr.first;
        m = pr.second;
    }
    
    int ishort = find_shortest_connection (from, to, d,
            alk_dat.vert2index_map, alk_dat.dmat,
            alk_dat.cl2index_map, m, l);
    // ishort is return value; an index into (from, to)
    merge_clusters (alk_dat.contig_mat,
            alk_dat.index2cl_map,
            alk_dat.cl2index_map, m, l);
    update_edgewt_maps (alk_dat, m, l);

    /* Cluster numbers start off here the same as vertex numbers, and so are
     * initially simple indices into the vert-by-vert matrices (contig_mat,
     * dmat, avg_dist, num_edges). As cluster form, numbers merge to one of the
     * pre-existing ones, so are still indexed into these same matrices which do
     * not change size. Cluster merging simply means that previous rows and
     * columns of these matrices will no longer be indexed, and all new indices
     * are derived from constantly updated values of index2cl_map and
     * cl2index_map.
     */

    for (auto cl: alk_dat.cl2index_map)
    {
        if (cl.first != l || cl.first != m)
        {
            const double tempd_l = alk_dat.avg_dist (cl.first, l),
                   tempd_m = alk_dat.avg_dist (cl.first, m);
            const unsigned int nedges_l = alk_dat.num_edges (cl.first, l),
                  nedges_m = alk_dat.num_edges (cl.first, m);

            alk_dat.avg_dist (cl.first, l) =
                (tempd_l * nedges_l + tempd_m * nedges_m) /
                static_cast <double> (nedges_l + nedges_m);
            alk_dat.num_edges (cl.first, l) = nedges_l + nedges_m;

            if (alk_dat.contig_mat (cl.first, l) == 1 ||
                    alk_dat.contig_mat (cl.first, m) == 1)
            {
                alk_dat.contig_mat (cl.first, l) = 1;
                if (tempd_l > 0.0)
                    tree.remove (tempd_l);
                if (tempd_m > 0.0)
                    tree.remove (tempd_m);
                
                double tempd = alk_dat.avg_dist (cl.first, l);
                if (tempd > 0.0)
                {
                    tree.insert (tempd);
                    alk_dat.edgewt2idx_pair_map [tempd] =
                        std::make_pair (cl.first, l);

                    std::unordered_set <double> wtset;
                    if (alk_dat.idx2edgewt_map.find (cl.first) ==
                            alk_dat.idx2edgewt_map.end ())
                        wtset.clear ();
                    else
                        wtset = alk_dat.idx2edgewt_map.at (cl.first);
                    wtset.emplace (tempd);
                    alk_dat.idx2edgewt_map.erase (cl.first);
                    alk_dat.idx2edgewt_map.emplace (cl.first, wtset);
                }
            } // end if C(c, l) = 1 or C(c, m) = 1 in Guo's terminology
        } // end if cl.first != (cfrom, cto)
    } // end for over cl

    return ishort;
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

    ALKDat alk_dat;
    BinarySearchTree tree;
    alk_init (alk_dat, tree, from, to, d);
    const unsigned int n = alk_dat.n;

    std::unordered_set <unsigned int> the_tree;
    while (the_tree.size () < (n - 1)) // tree has n - 1 edges
    {
        Rcpp::checkUserInterrupt ();
        int ishort = alk_step (alk_dat, tree, from, to, d);
        the_tree.insert (ishort);
    }

    std::vector <int> treevec (the_tree.begin (), the_tree.end ());

    return Rcpp::wrap (treevec);
}
