#include "common.h"
#include "utils.h"

unsigned int get_n (
    const Rcpp::IntegerVector &from,
    const Rcpp::IntegerVector &to)
{
    std::set <unsigned int> nodes;
    int max_node = 0;
    for (auto i: from)
    {
        nodes.insert (i);
        if (i > max_node)
            max_node = i;
    }
    for (auto i: to)
    {
        nodes.insert (i);
        if (i > max_node)
            max_node = i;
    }
    if (nodes.size () != (static_cast <unsigned int> (max_node) + 1))
        Rcpp::stop ("vertex numbers are discontinuous");

    return nodes.size ();
}

//' initial contiguity and distance matrices. The contiguity matrix is between
//' clusters, so is constantly modified, whereas the distance matrix is between
//' edges, so is fixed at load time.
//' @noRd
void mats_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        arma::Mat <unsigned short> &contig_mat,
        arma::Mat <double> &d_mat,
        const unsigned int n)
{
    contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);
    //d_mat = arma::zeros <arma::Mat <double> > (n, n);
    d_mat.resize (n, n);
    d_mat.fill (INFINITE_DOUBLE);

    for (int i = 0; i < from.length (); i++)
    {
        contig_mat (from [i], to [i]) = 1;
        contig_mat (to [i], from [i]) = 1;
        d_mat (from [i], to [i]) = d [i];
        d_mat (to [i], from [i]) = d [i];
    }
}

void dmat_full_init (
        const Rcpp::IntegerVector &from, // here, from_full, etc.
        const Rcpp::IntegerVector &to,
        const Rcpp::NumericVector &d,
        arma::Mat <double> &d_mat, // here, d_mat_full
        const unsigned int n)
{
    //d_mat = arma::zeros <arma::Mat <double> > (n, n);
    d_mat.resize (n, n);
    d_mat.fill (INFINITE_DOUBLE);

    for (int i = 0; i < from.length (); i++)
    {
        d_mat [from [i], to [i]] = d [i];
        d_mat [to [i], from [i]] = d [i];
    }
}

void sets_init (
        const Rcpp::IntegerVector &from,
        const Rcpp::IntegerVector &to,
        uint_map_t &vert2cl_map,
        uint_set_map_t &cl2vert_map)
{
    vert2cl_map.clear ();
    cl2vert_map.clear ();

    for (int i = 0; i < from.length (); i++)
    {
        std::set <unsigned int> eset;
        eset.insert (from [i]);
        cl2vert_map.emplace (from [i], eset);
    }
    const unsigned int n = get_n (from, to);
    // Initially assign all verts to clusters of same number:
    for (unsigned int i = 0; i < n; i++)
        vert2cl_map.emplace (i, i);
}
