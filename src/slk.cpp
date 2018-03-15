#include "common.h"
#include "slk.h"

void contig_mat_init (const Rcpp::DataFrame &gr,
        arma::Mat <unsigned short> &contig_mat)
{
    Rcpp::IntegerVector from = gr ["from"];
    Rcpp::IntegerVector to = gr ["to"];

    const int n = to.length ();

    contig_mat = arma::zeros <arma::Mat <unsigned short> > (n, n);

    for (int i = 0; i < n; i++)
    {
        contig_mat [from [i], to [i]] = 1;
        contig_mat [to [i], from [i]] = 1;
    }
}


//' rcpp_slk
//'
//' Full-order single linkage cluster redcap algorithm
//'
//' @noRd
void rcpp_slk (const Rcpp::DataFrame &grfull,
        Rcpp::DataFrame &gr)
{
    Rcpp::StringVector from = gr ["from"];
    Rcpp::StringVector to = gr ["to"];
    Rcpp::StringVector d = gr ["d"];

    arma::Mat <unsigned short> contig_mat;
    contig_mat_init (gr, contig_mat);
}
