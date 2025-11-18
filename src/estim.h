#ifndef ESTIM_H
#define ESTIM_H

Rcpp::List newAC_PGD(Rcpp::List& d1AC, arma::mat& Aold, arma::mat& Qmatrix,
                     arma::cube& Cold, double ss);

Rcpp::List newAD_PGD(Rcpp::List& d1AD, arma::mat& Aold, arma::mat& Qmatrix,
                     arma::cube& Dold, double ss);

Rcpp::List newAC_MD(Rcpp::List& d1AC, arma::mat& Aold, arma::mat& Qmatrix,
                    arma::cube& Cold, double& ssA, double& ssC);

Rcpp::List newAD_MD(Rcpp::List& d1AD, arma::mat& Aold, arma::mat& Qmatrix,
                    arma::cube& Dold, double& ssA, double& ssC);

Rcpp::List newAD_MDp(Rcpp::List& d1AD, arma::mat& Aold, arma::mat& Qmatrix,
                     arma::cube& Dold, double& ssA, double& ssC);

Rcpp::List newAD_MD_hess(Rcpp::List& d1AD,arma::mat& Aold, arma::mat& Qmatrix,
                         arma::cube& Dold, double& ssA, double& ssC);

Rcpp::List newAD_MD_adam(Rcpp::List& d1AD,arma::mat& Aold, arma::mat& Qmatrix,
                         arma::cube& Dold, double& ssA, double& ssC,
                         int iter, arma::vec& mt, arma::vec& vt, Rcpp::List& control);

arma::vec newM(arma::vec& mu, arma::mat& R, arma::mat& Z, double& ss);

arma::mat newL(arma::vec& mu, arma::mat& L, arma::mat& Z, double& ss, bool cor);

arma::mat newL_adam(arma::vec& mu, arma::mat& L, arma::mat& Z, double& ss, bool cor,
                    int iter, arma::vec& mt, arma::vec& vt, Rcpp::List& control);

void newpMR(arma::mat& Z, arma::mat& posM, arma::cube& posR, int iter);

arma::mat newZ_ULA(arma::mat& Y, arma::mat& PI, arma::mat& Z, arma::mat& A, arma::cube& C,
                   arma::vec& mu, arma::mat& L,
                   arma::mat& isd, double& h);

arma::mat newZ_MALA(arma::mat& Y, arma::mat& PI,arma::mat& Z, arma::mat& A, arma::cube& C,
                    arma::vec& mu, arma::mat& L, double& ar,
                    arma::mat& isd, double& h, arma::vec& knots, const unsigned int degree, const std::string& basis);

arma::mat newZ_RWMH(arma::mat& Y, arma::mat& PI,arma::mat& Z, arma::mat& A, arma::cube& C,
                    arma::vec& mu, arma::mat& L, double& ar,
                    const double h, arma::vec& knots, const unsigned int degree, const std::string& basis);

arma::mat newG_MD(Rcpp::List& d1G, arma::mat& Gold, double& ss);

arma::mat newG_PL2(Rcpp::List& d1G, arma::mat& Gold, double& ss);

arma::mat newG_PL2_adam(Rcpp::List& d1G,arma::mat& Gold, double& ss,
                        int iter, arma::vec& mt, arma::vec& vt, Rcpp::List& control);

arma::mat newZ_ULA_aCDM(arma::mat& Y, arma::mat& PI, arma::mat& Z,
                        arma::mat& G, arma::vec& mu, arma::mat& L, double& h);

arma::mat newZ_RWMH_aCDM(arma::mat& Y, arma::mat& PI, arma::mat& Z,
                         arma::mat& G, arma::vec& mu, arma::mat& L, const double h, double& ar);

arma::mat newZ_MALA_aCDM(arma::mat& Y, arma::mat& PI, arma::mat& Z,
                         arma::mat& G, arma::vec& mu, arma::mat& L, double& h, double& ar);

double fy_gapmCDM_IS(arma::mat& Y, arma::mat& A, arma::cube& C,
                     arma::vec& mu, arma::mat& L, arma::mat& Z,
                     arma::mat& pM, arma::cube& pR, Rcpp::List& control);

double fy_aCDM_IS(arma::mat& Y, arma::mat& G,
                  arma::vec& mu, arma::mat& L, arma::mat& Z,
                  arma::mat& pM, arma::cube& pR, Rcpp::List& control);

#endif
