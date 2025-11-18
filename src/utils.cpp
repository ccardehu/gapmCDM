#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <splines2Armadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(splines2)]]

using namespace Rcpp;
using namespace arma;

arma::mat check_inv_sympd(arma::mat& Lnew, arma::mat& Lold){
  arma::mat Rnew = Lnew*Lnew.t();
  arma::mat check(arma::size(Rnew));
  bool success = arma::inv_sympd(check,Rnew);
  if(success){
    return Lnew;
  } else {
    return Lold;
  }
}

arma::mat vR2Rc(arma::mat& R){
  arma::mat Is = arma::diagmat(arma::sqrt(1/arma::diagvec(R)));
  arma::mat out = Is * R * Is;
  return(out);
}

arma::mat Z2U(arma::mat& Z){
  const int n = Z.n_rows;
  const int p = Z.n_cols;
  arma::mat out(n,p);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      out(i,j) = R::pnorm(Z(i,j), 0.0, 1.0, 1, 0);
    }
  }
  return(out);
}

// [[Rcpp::export]]
arma::cube D2C(arma::cube& D){
  const int p = D.n_rows;
  const int q = D.n_slices;
  arma::cube out(arma::size(D));
  for(int j = 0; j < q; j++){
    for(int i = 0; i < p; i++){
      out.slice(j).row(i) = arma::cumsum(D.slice(j).row(i));
    }
  }
  return(out);
}

arma::mat cumsumMat(arma::vec& expd){
  const int tp = expd.n_elem;
  arma::vec outv(tp);
  for(int i = 1; i < tp; i++){
    outv(i-1) = arma::accu(expd(arma::span(i,tp-1)));
  }
  arma::mat out(tp,tp);
  out.each_row() = outv.t();
  return(out);
}

arma::vec ProxD(arma::vec& y){
  // Taken from https://arxiv.org/pdf/1309.1541.pdf
  // See also https://angms.science/doc/CVX/Proj_simplex.pdf
  // and https://arxiv.org/pdf/1101.6081.pdf
  const int m = y.size() ;
  const arma::vec u = arma::sort(y, "descend") ;
  arma::vec rho(m), out1(m), out2(m);
  for(int i = 0; i < m; i++){
    double ii = i+1;
    out1(i) = u(i) + (1/ii)* (1 - arma::accu(u(arma::span(0,i)))) ;
  }
  arma::uvec idx = arma::find(out1 > 0) ;
  arma::uword id = idx.max() ;
  double lambda = out1(id) - u(id) ;
  for(int i = 0; i < m; i++){
    out2(i) = (y(i) + lambda > 0) ? (y(i) + lambda) : 0 ;
  }
  return(out2);
}

Rcpp::List QP_Rcpp(Rcpp::NumericMatrix& Dmat,
                   Rcpp::NumericVector& dvec,
                   Rcpp::NumericMatrix& Amat,
                   Rcpp::NumericVector& bvec,
                   int meq = 0L,
                   bool factor = false) {

  Rcpp::Environment quadprogR = Rcpp::Environment::namespace_env("quadprog");
  Rcpp::Function solveQP = quadprogR["solve.QP"];
  Rcpp::List out = solveQP(Rcpp::_["Dmat"] = Dmat,
                           Rcpp::_["dvec"] = dvec,
                           Rcpp::_["Amat"] = Amat,
                           Rcpp::_["bvec"] = bvec,
                           Rcpp::_["meq"] = meq,
                           Rcpp::_["factorized"] = factor);

  return out;
}

arma::vec QP_Arma(arma::mat& Dmat,
                  arma::vec& dvec,
                  arma::mat& Amat,
                  arma::vec& bvec,
                  const int meq = 0L,
                  bool factor = false) {

  int n = dvec.n_elem;
  int m = bvec.n_elem;
  Rcpp::NumericMatrix DmatT = Rcpp::NumericMatrix(n,n,Dmat.begin());
  Rcpp::NumericVector dvecT = Rcpp::NumericVector(dvec.begin(),dvec.end());
  Rcpp::NumericMatrix AmatT = Rcpp::NumericMatrix(n,m,Amat.begin());
  Rcpp::NumericVector bvecT = Rcpp::NumericVector(bvec.begin(),bvec.end());
  Rcpp::List list = QP_Rcpp(DmatT,dvecT,AmatT,bvecT,meq,factor);
  Rcpp::NumericVector out0 = list["solution"];
  arma::vec out = Rcpp::as<arma::vec>(out0);
  return out;
}

arma::vec ProjSim(arma::vec& y){
  const int np = y.n_elem;
  arma::mat Dmat = arma::eye(np,np);
  arma::rowvec vtmp(np, arma::fill::value(-1.0));
  arma::mat AmatT = arma::join_cols(arma::eye(np,np),vtmp);
  arma::mat Amat = AmatT.t();
  arma::vec bvec(np + 1);
  bvec.tail(1) = -1.0;
  arma::vec out = QP_Arma(Dmat,y,Amat,bvec,0L,false);
  return(out);
  return(bvec);
}

arma::mat ProxL(arma::mat& L){
  int q = L.n_rows;
  arma::mat Lout(q,q);
  for(int i = 0; i < q; i++){
    double cross = arma::norm(L.row(i),2L);
    if(cross > 0){
      Lout.row(i) = L.row(i)/cross;
    } else {
      Rcpp::stop("L has norm zero");
    }
  }
  return(Lout);
}

arma::mat fixL(arma::mat& Lnew, arma::mat& Lold){
  int q = Lnew.n_rows;
  for(int i = 0; i < q; i++){
    if(!Lnew.row(i).is_finite()) Lnew.row(i) = Lold.row(i);
    arma::rowvec tmpr = arma::abs(Lnew.row(i) - Lold.row(i));
    if(arma::any(tmpr > 0.2)) Lnew.row(i) = Lold.row(i);
  }
  if(!(Lnew*Lnew.t()).is_sympd()) Lnew = Lold;
  return(Lnew);
}

void cube2eye(arma::cube& C){
  const int q = C.n_rows;
  C.fill(0.0);
  C.each_slice() += arma::eye(q,q);
}

arma::mat rmvNorm(const int n, arma::vec& mu, arma::mat& R){
  const int q = R.n_cols;
  arma::mat Y(n,q), eigm(q,q);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < q; j++){
      Y(i,j) = R::rnorm(0,1);
    }
  }
  arma::vec eigv(q);
  arma::eig_sym(eigv,eigm,R);
  arma::mat teigm = eigm.t();
  teigm.each_col() %= arma::sqrt(eigv);
  return arma::repmat(mu, 1, n).t() + Y * (eigm * teigm);
}

arma::mat rmvNorm_IS(arma::mat& posM, arma::cube& posR){
  const int q = posR.n_cols;
  const int n = posM.n_rows;
  arma::mat Y(n,q), eigm(q,q);
  arma::vec eigv(q);
  for(int i = 0; i < n; i++){
    arma::mat posRi = posR.slice(i);
    arma::eig_sym(eigv,eigm,posRi);
    arma::mat teigm = eigm.t();
    teigm.each_col() %= arma::sqrt(eigv);
    for(int j = 0; j < q; j++){
      Y(i,j) = R::rnorm(0,1);
    }
    Y.row(i) = posM.row(i) + Y.row(i) * (eigm * teigm);
  }
  return(Y);
}

arma::mat rmvStNorm(const int n, const int q){
  arma::mat Y(n,q) ;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < q; j++){
      Y(i,j) = R::rnorm(0,1);
    }
  }
  return (Y);
}

double dmvNorm(arma::vec& y, arma::vec& mu, arma::mat& R, const bool log = true){
  const int q = y.size() ;
  double lR = arma::log_det_sympd(R);
  double out = -0.5*(q*std::log(2*M_PI) + lR + arma::as_scalar((y - mu).t()*(arma::solve(R, y - mu,arma::solve_opts::likely_sympd)))) ;
  if(log) return(out);
  else return(std::exp(out));
}

arma::mat dmvStNorm(arma::mat& Z){
  const int n = Z.n_rows;
  const int q = Z.n_cols;
  arma::mat out(n,q) ;
  for(int j = 0; j < q; j++){
    for(int m = 0; m < n; m++){
      out(m,j) = R::dnorm(Z(m,j),0.0,1.0,0);
    }
  }
  return (out);
}

arma::mat Z2Zo(arma::mat& Z, arma::vec& mu, arma::mat& Lt){
  const int n = Z.n_rows;
  arma::mat iL = arma::inv(arma::trimatu(Lt));
  arma::mat muV = arma::repmat(mu, 1, n).t();
  arma::mat out = (Z - muV) * iL;
  return (out);
}

arma::mat UA(arma::mat& U, arma::mat& Apat){
  const int n = U.n_rows;
  const int K = Apat.n_rows;
  arma::mat out(n,K);
  for(int k = 0; k < K; k++){
    arma::mat tAk(arma::size(U));
    tAk.each_row() = Apat.row(k);
    out.col(k) = arma::prod(arma::pow(U,tAk) % arma::pow(1.0 - U, 1.0 - tAk), 1);
  }
  return(out);
}

arma::mat dUA(arma::mat& U, arma::rowvec& Ak){
  arma::mat tAk(arma::size(U));
  tAk.each_row() = Ak;
  return (tAk/U - (1-tAk)/(1-U));
}

// [[Rcpp::export]]
Rcpp::List genpar(const int p, const int q, const int tp,
                  const double probSparse, arma::mat& Qmatrix,
                  const std::string& basis){
  arma::mat A(p,q);
  arma::cube C(p,tp,q);
  arma::cube D(p,tp,q);
  arma::vec sparse(p);
  for(int j = 0; j < p; j++){
    if(q > 1){
      if(arma::accu(Qmatrix) == p*q){
        sparse(j) = R::rbinom(1,probSparse);
        if(sparse(j) == 1){
          arma::uvec iA = arma::regspace<uvec>(0,q-1);
          int i2A = arma::as_scalar(Rcpp::RcppArmadillo::sample(iA,1,false,1/iA.size()));
          A(j,i2A) = 1;
        } else {
          arma::rowvec tA(q+1);
          tA.tail(1) = 1;
          Rcpp::NumericVector tARcpp = Rcpp::runif(q-1);
          tA(span(1,tARcpp.size())) = arma::sort(Rcpp::as<arma::rowvec>(tARcpp));
          tA = arma::diff(tA);
          A.row(j) = tA;
        }
      } else {
        arma::rowvec Qj = Qmatrix.row(j);
        arma::uvec iQj = arma::find(Qj == 1);
        if(iQj.n_elem > 1){
          arma::rowvec tA(iQj.n_elem + 1);
          tA.tail(1) = 1;
          Rcpp::NumericVector tARcpp = Rcpp::runif(iQj.n_elem-1);
          tA(span(1,tARcpp.size())) = arma::sort(Rcpp::as<arma::rowvec>(tARcpp));
          tA = arma::diff(tA);
          arma::rowvec tAj(q);
          tAj(iQj) = tA;
          A.row(j) = tAj;
        } else {
          A.row(j) = Qj;
        }
      }
    } else {
      A.fill(1.0);
    }

    if(basis == "is"){
      for(int k = 0; k < q; k++){
        arma::rowvec tC(tp+1);
        tC.tail(1) = 1;
        Rcpp::NumericVector tCRcpp = Rcpp::runif(tp-1);
        tC(span(1,tCRcpp.size())) = arma::sort(Rcpp::as<arma::rowvec>(tCRcpp));
        tC = arma::diff(tC);
        C.slice(k).row(j) = tC;
      }
    } else {
      for(int k = 0; k < q; k++){
        arma::rowvec tD(tp+1);
        tD.tail(1) = 1;
        Rcpp::NumericVector tDRcpp = Rcpp::runif(tp-1);
        tD(span(1,tDRcpp.size())) = arma::sort(Rcpp::as<arma::rowvec>(tDRcpp));
        tD = arma::diff(tD);
        D.slice(k).row(j) = tD;
      }
    }
  }
  if(basis != "is") C = D2C(D);
  return Rcpp::List::create(Rcpp::Named("A") = A,
                            Rcpp::Named("C") = C,
                            Rcpp::Named("D") = D);
}

// [[Rcpp::export]]
arma::mat genpar_aCDM(arma::mat& Qmatrix, const double maxG0, const double maxSl){
  const int p = Qmatrix.n_rows;
  const int q = Qmatrix.n_cols;
  arma::mat G(p,q+1);
  for(int j = 0; j < p; j++){
    arma::rowvec Qj = Qmatrix.row(j);
    arma::uvec iQj = arma::find(Qj == 1) + 1;
    arma::uvec iGj(iQj.n_elem + 1);
    iGj.subvec(1, iQj.n_elem) = iQj;
    double tmpGuess = R::runif(0.0,maxG0);
    double tmpSlip = R::runif(0.0,maxSl);
    if(iQj.n_elem > 1){
      arma::rowvec tG(iGj.n_elem);
      tG.tail(1) = 1.0 - tmpGuess - tmpSlip;
      Rcpp::NumericVector tGRcpp = Rcpp::runif(iQj.n_elem-1,0.0,1.0 - tmpGuess - tmpSlip);
      tG(arma::span(1,tGRcpp.size())) = arma::sort(Rcpp::as<arma::rowvec>(tGRcpp));
      tG = arma::diff(tG);
      arma::rowvec t2G(tG.n_elem + 1);
      t2G(0) = tmpGuess;
      t2G.subvec(1, tG.n_elem) = tG;
      arma::rowvec tGj(q+1);
      tGj(iGj) = t2G;
      G.row(j) = tGj;
    } else {
      arma::rowvec tG(iGj.n_elem);
      tG(0) = tmpGuess;
      tG(1) = 1.0 - tmpGuess - tmpSlip;
      arma::rowvec tGj(q+1);
      tGj(iGj) = tG;
      G.row(j) = tGj;
    }
  }
  return(G);
}

// [[Rcpp::export]]
arma::cube SpU_isp(arma::mat& U, arma::vec& knots, const unsigned int deg){
  const arma::vec bouK = {0.0,1.0};
  const int n = U.n_rows;
  const int q = U.n_cols;
  const int np = knots.size() + deg ;
  arma::mat out1(n,q*np);
  arma::mat out2(n,q*np);
  arma::cube out(n,q*np,2);
  for(int j = 0; j < q; j++){
    splines2::ISpline iSp(U.col(j),knots,deg,bouK);
    out1.cols(j*np,(j+1)*np-1) = iSp.basis(false);
    out2.cols(j*np,(j+1)*np-1) = iSp.derivative(1,false);
  }
  out.slice(0) = out1;
  out.slice(1) = out2;
  return(out)  ;
}

// [[Rcpp::export]]
arma::cube SpU_bsp(arma::mat& U, arma::vec& knots, const unsigned int deg){
  const arma::vec bouK = {0.0,1.0};
  const int n = U.n_rows;
  const int q = U.n_cols;
  const int np = knots.size() + deg ;
  arma::mat out1(n,q*np);
  arma::mat out2(n,q*np);
  arma::cube out(n,q*np,2);
  for(int j = 0; j < q; j++){
    splines2::BSpline bSp(U.col(j),knots,deg,bouK);
    out1.cols(j*np,(j+1)*np-1) = bSp.basis(false);
    out2.cols(j*np,(j+1)*np-1) = bSp.derivative(1,false);
  }
  out.slice(0) = out1;
  out.slice(1) = out2;
  return(out)  ;
}

// [[Rcpp::export]]
arma::mat prob(arma::mat& A, arma::cube& C, arma::mat& ism){
  const int n = ism.n_rows;
  const int p = A.n_rows;
  const int q = C.n_slices;
  const int np = C.n_cols;
  arma::mat out(n,p, arma::fill::zeros);
  for(int i = 0; i < p; i++){
    for(int j = 0; j < q; j++){
      out.col(i) += A(i,j) * ism.cols(j*np,(j+1)*np-1)*(C.slice(j).row(i)).t();
    }
  }
  return(out);
}

arma::mat fyz(arma::mat& Y, arma::mat& PI){
  const int n = Y.n_rows;
  const int p = Y.n_cols;
  arma::mat out(n,p,arma::fill::zeros);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      if(std::isnan(Y(i,j))) continue;
      out(i,j) = (Y(i,j) > 0) ? std::max(arma::datum::log_min, std::log(PI(i,j))) : std::max(arma::datum::log_min, std::log(1-PI(i,j))) ;
    }
  }
  return(out);
}

arma::mat fz(arma::mat& Z, arma::vec& mu, arma::mat& R){
  const int n = Z.n_rows;
  arma::vec out(n);
  for(int i = 0; i < n; i++){
    arma::vec zrow = Z.row(i).t();
    out(i) = std::max(arma::datum::log_min, dmvNorm(zrow, mu, R, true)) ;
  }
  return(out);
}

arma::mat fz_IS(arma::mat& Z, arma::mat& pM, arma::cube& pR){
  const int n = Z.n_rows;
  arma::vec out(n);
  for(int i = 0; i < n; i++){
    arma::vec zrow = Z.row(i).t();
    arma::vec mu = pM.row(i).t();
    arma::mat R = pR.slice(i);
    out(i) = std::max(arma::datum::log_min, dmvNorm(zrow, mu, R, true)) ;
  }
  return(out);
}

// // [[Rcpp::export]]
// double fy_gapmCDM(arma::mat& Y, arma::mat& A, arma::cube& C,
//                   arma::vec& mu, arma::mat& R, Rcpp::List& control){
//
//   const unsigned int degree = control["degree"];
//   const int nsim = control["nsim"];
//   arma::vec knots = control["knots"];
//   const bool verbose = control["verbose"];
//   const std::string basis = Rcpp::as<std::string>(control["basis"]);
//
//   const int n = Y.n_rows;
//   if(verbose) Rcpp::Rcout << "\n Calculating marginal log-likelihood ... ";
//   double mllk = 0;
//   arma::mat Eobj(n,nsim);
//   arma::mat Zsim(n,R.n_cols);
//   arma::mat isMo(n,R.n_cols * C.n_cols);
//   arma::cube spObj(n,R.n_cols * C.n_cols,2);
//   for(int ii = 0; ii < nsim; ii++){
//     if (ii % 2 == 0) Rcpp::checkUserInterrupt();
//     Zsim = rmvNorm(n,mu,R);
//     arma::mat Usim = Z2U(Zsim);
//     if(basis == "is"){
//       spObj = SpU_isp(Usim,knots,degree);
//       isMo = spObj.slice(0);
//     } else {
//       spObj = SpU_bsp(Usim,knots,degree);
//       isMo = spObj.slice(0);
//     }
//     arma::mat piH = prob(A,C,isMo);
//     Eobj.col(ii) = arma::sum(fyz(Y,piH),1);
//   }
//   arma::vec maxEobj(n);
//   for(int m = 0; m < n; m++){
//     maxEobj(m) = arma::max(Eobj.row(m));
//   }
//   Eobj.each_col() -= maxEobj;
//   arma::mat EEobj = arma::exp(Eobj);
//   arma::vec V1 = arma::log(arma::sum(EEobj,1));
//   mllk = arma::accu(V1) + arma::accu(maxEobj) - n*std::log(nsim);
//   if(verbose) Rcpp::Rcout << "\r Calculating marginal log-likelihood ... (m-llk: " << std::to_string(mllk) << ")";
//   return(mllk);
// }

Rcpp::List d1AC(arma::mat& Y, arma::mat& PI, arma::mat& ism,
                arma::mat& A, arma::cube& C){
  arma::mat PIm1PI = PI % (1-PI);
  arma::mat YmPI = (Y - PI)/(PIm1PI.clamp(arma::datum::eps, 1.0));
  const int p = YmPI.n_cols;
  const int q = C.n_slices;
  const int np = C.n_cols;
  const int tp = A.size() + C.size();
  arma::uvec ia = arma::regspace<arma::uvec>(0, A.size());
  arma::umat iA(ia.begin(), A.n_rows, A.n_cols, false);
  arma::uvec ic = arma::regspace<arma::uvec>(A.size(), C.size() + A.size());
  arma::ucube iC(ic.begin(), C.n_rows, C.n_cols, C.n_slices, false);
  arma::vec out(tp);
  for(int i = 0; i < p; i++){
    arma::vec YmPIcol = YmPI.col(i);
    arma::uvec noNA = arma::find_finite(YmPIcol);
    for(int j = 0; j < q; j++){
      arma::uvec jcols = arma::regspace<uvec>(j*np,(j+1)*np-1);
      arma::uword idA = iA(i,j);
      arma::uvec idC = iC.slice(j).row(i).t();
      out(idA) = arma::as_scalar(YmPIcol.elem(noNA).t() * (ism(noNA,jcols)*(C.slice(j).row(i)).t()));
      arma::vec A2 = ism(noNA,jcols).t() * YmPIcol.elem(noNA);
      out(idC) = A(i,j) * A2 ; //A2.t()
    }
  }
  return Rcpp::List::create(Rcpp::Named("grad") = out,
                            Rcpp::Named("iA") = iA,
                            Rcpp::Named("iC") = iC);
}

// Rcpp::List d2AC(arma::mat& YmPI, arma::mat& ism,
//                 arma::mat& A, arma::cube& C){
//   const int p = YmPI.n_cols;
//   const int q = C.n_slices;
//   const int np = C.n_cols;
//   const int tp = A.size() + C.size();
//   arma::uvec ia = arma::regspace<arma::uvec>(0, A.size());
//   arma::umat iA(ia.begin(), A.n_rows, A.n_cols, false);
//   arma::uvec ic = arma::regspace<arma::uvec>(A.size(), C.size() + A.size());
//   arma::ucube iC(ic.begin(), C.n_rows, C.n_cols, C.n_slices, false);
//   arma::mat Asq = arma::pow(A,2);
//   arma::mat out(tp,tp);
//   for(int i = 0; i < p; i++){
//     for(int j = 0; j < q; j++){
//       arma::uword idA = iA(i,j);
//       arma::uvec idC = iC.slice(j).row(i).t();
//       out(idA,idA) = arma::as_scalar((ism.cols(j*np,(j+1)*np-1)*(C.slice(j).row(i)).t()).t()*arma::diagmat(YmPI.col(i))*(ism.cols(j*np,(j+1)*np-1)*(C.slice(j).row(i)).t()));
//       arma::mat A2 = ism.cols(j*np,(j+1)*np-1).t() * arma::diagmat(YmPI.col(i)) * ism.cols(j*np,(j+1)*np-1);
//       out(idC,idC) =  Asq(i,j) * A2 ;
//     }
//   }
//   return Rcpp::List::create(Rcpp::Named("hess") = out,
//                             Rcpp::Named("iA") = iA,
//                             Rcpp::Named("iC") = iC);
// }

// arma::mat d1CdD(arma::vec& d){
//   const int tp = d.n_elem;
//   arma::vec expd = arma::exp(d);
//   arma::mat d1(tp,tp);
//   arma::vec N = arma::cumsum(expd);
//   double D = arma::accu(expd);
//
//   arma::mat expDM(tp,tp);
//   expDM.each_col() = expd;
//   arma::mat Nm(tp,tp);
//   Nm.each_row() = N.t();
//   arma::mat DmN = D - Nm;
//   Nm %= expDM;
//   DmN %= expDM;
//   arma::uvec iLm = arma::trimatl_ind(arma::size(d1),-1);
//   arma::uvec iUm = arma::trimatu_ind(arma::size(d1));
//   d1(iUm) = DmN(iUm);
//   d1(iLm) = -Nm(iLm);
//
//   return(d1/std::pow(D,2));
// }

arma::mat d1CdD(arma::vec& d){
  int tp = d.n_elem;
  arma::mat out(tp,tp);
  arma::uvec iO = arma::trimatl_ind(arma::size(out));
  out(iO).fill(1.0);
  return(out);
}

Rcpp::List dCdD(arma::vec& d){
  const int tp = d.n_elem;
  arma::vec expd = arma::exp(d);
  arma::mat d1(tp,tp);
  arma::vec N = arma::cumsum(expd);
  double D = arma::accu(expd);

  // For d1CdD
  arma::mat expDM(tp,tp);
  expDM.each_col() = expd;
  arma::mat Nm(tp,tp);
  Nm.each_row() = N.t();
  arma::mat DmN = D - Nm;
  Nm %= expDM;
  DmN %= expDM;
  arma::uvec iLm = arma::trimatl_ind(arma::size(d1),-1);
  arma::uvec iUm = arma::trimatu_ind(arma::size(d1));
  d1(iUm) = DmN(iUm);
  d1(iLm) = -Nm(iLm);

  // // For d2CdD
  arma::mat NmDmN = Nm - DmN;
  NmDmN %= expDM;
  arma::cube d2(tp,tp,tp);
  arma::mat tM2(tp,tp,arma::fill::value(22));
  tM2(iLm).fill(11);

  for(int i = 0; i < tp; i++){
    arma::mat tM1(tp,i,arma::fill::value(22));
    arma::uvec tiUm = arma::trimatu_ind(arma::size(tM1));
    tM1(tiUm).fill(11);
    arma::mat tM2a = tM2.cols(i,tp-1);
    arma::mat tmat = arma::join_rows(tM1,tM2a);
    tmat.shed_row(i);
    arma::rowvec ri = -d1.row(i)*(expd(i)-(D-expd(i)));
    tmat.insert_rows(i,ri);
    arma::uvec i11 = arma::find(tmat == 11);
    arma::uvec i22 = arma::find(tmat == 22);
    tmat(i22) = -2.0*d1(i22)*expd(i);
    tmat(i11) = NmDmN(i11) * expd(i);
    d2.slice(i) = tmat;
    d2.slice(i).tail_cols(1).fill(0.0);
  }
  return Rcpp::List::create(Rcpp::Named("d1CdD") = d1/std::pow(D,2),
                            Rcpp::Named("d2CdD") = d2/std::pow(D,3));
}

Rcpp::List d1AD(arma::mat& Y, arma::mat& PI, arma::mat& ism,
                arma::mat& A, arma::cube& D){
  arma::mat PIm1PI = PI % (1-PI);
  arma::mat YmPI = (Y - PI)/(PIm1PI.clamp(arma::datum::eps, 1.0));
  const int p = YmPI.n_cols;
  const int q = D.n_slices;
  const int np = D.n_cols;
  const int tp = A.size() + D.size();
  arma::cube C = D2C(D);
  arma::uvec ia = arma::regspace<arma::uvec>(0, A.size());
  arma::umat iA(ia.begin(), A.n_rows, A.n_cols, false);
  arma::uvec ic = arma::regspace<arma::uvec>(A.size(), C.size() + A.size());
  arma::ucube iC(ic.begin(), C.n_rows, C.n_cols, C.n_slices, false);
  arma::vec out(tp);
  for(int i = 0; i < p; i++){
    arma::vec YmPIcol = YmPI.col(i);
    arma::uvec noNA = arma::find_finite(YmPIcol);
    for(int j = 0; j < q; j++){
      arma::uvec jcols = arma::regspace<uvec>(j*np,(j+1)*np-1);
      arma::uword idA = iA(i,j);
      arma::uvec idC = iC.slice(j).row(i).t();
      out(idA) = arma::as_scalar(YmPIcol.elem(noNA).t() * (ism(noNA,jcols)*(C.slice(j).row(i)).t()));
      arma::vec Dij = D.slice(j).row(i).t();
      arma::mat dD1 = d1CdD(Dij).t();
      arma::vec V1 = ism(noNA,jcols).t() * YmPIcol.elem(noNA);
      arma::vec V2 = (dD1 * V1);
      out(idC) = A(i,j) * V2 ;
    }
  }
  return Rcpp::List::create(Rcpp::Named("gs") = out,
                            Rcpp::Named("iA") = iA,
                            Rcpp::Named("iD") = iC);
}

Rcpp::List dAD(arma::mat& Y, arma::mat& PI, arma::mat& ism,
               arma::mat& A, arma::cube& D){
  arma::mat PIm1PI = PI % (1-PI);
  arma::mat YmPI = (Y - PI)/(PIm1PI.clamp(arma::datum::eps, 1.0));
  arma::mat YmPI2 = (-(arma::pow(PI,2) - 2.0*(PI%Y) + Y)/(arma::pow(PI-1,2)%arma::pow(PI,2)));
  const int p = YmPI.n_cols;
  const int q = D.n_slices;
  const int np = D.n_cols;
  const int tp = A.size() + D.size();
  arma::cube C = D2C(D);
  arma::uvec ia = arma::regspace<arma::uvec>(0, A.size());
  arma::umat iA(ia.begin(), A.n_rows, A.n_cols, false);
  arma::uvec ic = arma::regspace<arma::uvec>(A.size(), C.size() + A.size());
  arma::ucube iC(ic.begin(), C.n_rows, C.n_cols, C.n_slices, false);
  arma::vec out(tp);
  arma::mat Asq = arma::pow(A,2);
  arma::mat out2(tp,tp);
  for(int i = 0; i < p; i++){
    arma::vec YmPIcol = YmPI.col(i);
    arma::vec YmPI2col = YmPI2.col(i);
    arma::uvec noNA = arma::find_finite(YmPIcol);
    for(int j = 0; j < q; j++){
      arma::uvec jcols = arma::regspace<uvec>(j*np,(j+1)*np-1);
      arma::uword idA = iA(i,j);
      arma::uvec idC = iC.slice(j).row(i).t();
      out(idA) = arma::as_scalar(YmPIcol.elem(noNA).t() * (ism(noNA,jcols)*(C.slice(j).row(i)).t()));
      arma::vec Dij = D.slice(j).row(i).t();
      Rcpp::List dD = dCdD(Dij);
      arma::mat dD1 = dD["d1CdD"];
      arma::cube dD2 = dD["d2CdD"];
      arma::mat M1(arma::size(dD1));
      arma::vec V1 = ism(noNA,jcols).t() * YmPIcol.elem(noNA);
      arma::vec V2 = (dD1 * V1);
      for(int i = 0; i < np; i++){
        M1 += dD2.col_as_mat(i)*V1(i);
      }
      arma::mat M2 = (dD1 * ism(noNA,jcols).t()) * arma::diagmat(YmPI2col.elem(noNA)) * (ism(noNA,jcols) * dD1.t());
      out(idC) = A(i,j) * V2 ;
      out2(idC,idC) = Asq(i,j) * M2 +  A(i,j) * M1;
    }
  }
  return Rcpp::List::create(Rcpp::Named("gs") = out,
                            Rcpp::Named("hs") = out2,
                            Rcpp::Named("iA") = iA,
                            Rcpp::Named("iD") = iC);
}

arma::cube d1PIdZ(arma::mat& A, arma::cube& C, arma::mat& isd, arma::mat& Z, arma::mat& sM){
  const int p = A.n_rows;
  const int q = C.n_slices;
  const int n = isd.n_rows;
  const int np = C.n_cols;
  arma::mat dZ = dmvStNorm(Z) * sM;
  arma::cube out(n,q,p);
  for(int i = 0; i < p; i++){
    for(int j = 0; j < q; j++){
      out.slice(i).col(j) = A(i,j) * (isd.cols(j*np,(j+1)*np-1) * C.slice(j).row(i).t()) % dZ.col(j) ;
    }
  }
  return(out);
}

arma::cube d1PIdZ_aCDM(arma::mat& G, arma::mat& Z, arma::mat& sM){
  const int p = G.n_rows;
  const int q = Z.n_cols;
  const int n = Z.n_rows;
  arma::mat dZ = dmvStNorm(Z) * sM;
  arma::cube out(n,q,p);
  arma::uvec iGi = arma::regspace<arma::uvec>(1, G.n_cols-1);
  for(int i = 0; i < p; i++){
    arma::rowvec Gi = G.row(i);
    arma::mat tG(arma::size(Z));
    tG.each_row() = Gi(iGi).t();
    out.slice(i) = (tG % dZ) ;
  }
  return(out);
}

arma::mat d1PostZ(arma::mat& YmPI, arma::mat& Z, arma::mat& isd,
                  arma::mat& A, arma::cube& C, arma::vec& mu, arma::mat& L){
  const int n = YmPI.n_rows;
  const int p = YmPI.n_cols;
  const int q = C.n_slices;
  arma::mat out(n,q);
  arma::mat Ri = L * L.t();
  arma::mat sM = arma::diagmat(arma::sqrt(arma::diagvec(Ri)));
  arma::mat iR = arma::inv_sympd(Ri);
  arma::cube dpidz = d1PIdZ(A,C,isd,Z,sM);
  for(int i = 0; i < p; i++){
    arma::mat temp = dpidz.slice(i).each_col() % YmPI.col(i);
    for(int j = 0; j < q; j++){
      for(int m = 0; m < n; m++){
        if(std::isnan(temp(m,j))) continue;
        out(m,j) += temp(m,j);
      }
    }
  }
  arma::mat Zcent(Z);
  Zcent.each_row() -= mu.t();
  out -= Zcent*iR;
  return(out);
}

// [[Rcpp::export]]
arma::mat prob_aCDM(arma::mat& G, arma::mat& U){
  const int n = U.n_rows;
  const int p = G.n_rows;
  arma::vec c1(n, arma::fill::value(1.0));
  arma::mat U1 = arma::join_rows(c1,U);
  arma::mat out(n,p, arma::fill::zeros);
  for(int i = 0; i < p; i++){
    arma::rowvec Gi = G.row(i);
    out.col(i) = U1 * Gi.t();
  }
  return(out);
}


// // [[Rcpp::export]]
// double fy_aCDM(arma::mat& Y, arma::mat& G, arma::mat& Qmatrix, arma::mat& Apat,
//                arma::vec& mu, arma::mat& R, Rcpp::List& control){
//
//   const int nsim = control["nsim"];
//   const bool verbose = control["verbose"];
//
//   const int n = Y.n_rows;
//   if(verbose) Rcpp::Rcout << "\n Calculating marginal log-likelihood ... ";
//   double mllk = 0;
//   arma::mat Eobj(n,nsim);
//   arma::mat Zsim(n,R.n_cols);
//   for(int ii = 0; ii < nsim; ii++){
//     if (ii % 2 == 0) Rcpp::checkUserInterrupt();
//     Zsim = rmvNorm(n,mu,R);
//     arma::mat Usim = Z2U(Zsim);
//     arma::mat piH = prob_aCDM(G,Usim,Apat);
//     Eobj.col(ii) = arma::sum(fyz(Y,piH),1);
//   }
//   arma::vec maxEobj(n);
//   for(int m = 0; m < n; m++){
//     maxEobj(m) = arma::max(Eobj.row(m));
//   }
//   Eobj.each_col() -= maxEobj;
//   arma::mat EEobj = arma::exp(Eobj);
//   arma::vec V1 = arma::log(arma::sum(EEobj,1));
//   mllk = arma::accu(V1) + arma::accu(maxEobj) - n*std::log(nsim);
//   if(verbose) Rcpp::Rcout << "\r Calculating marginal log-likelihood ... (m-llk: " << std::to_string(mllk) << ")";
//   return(mllk);
// }

Rcpp::List d1G(arma::mat& Y, arma::mat& U, arma::mat& PI, arma::mat& G){
  arma::mat PIm1PI = PI % (1-PI);
  arma::mat YmPI = (Y - PI)/(PIm1PI.clamp(arma::datum::eps, 1.0));
  const int p = Y.n_cols;
  const int n = Y.n_rows;
  arma::vec c1(n, arma::fill::value(1.0));
  arma::mat U1 = arma::join_rows(c1,U);
  arma::mat tout(arma::size(G));
  for(int i = 0; i < p; i++){
    arma::vec YmPIcol = YmPI.col(i);
    arma::uvec noNA = arma::find_finite(YmPIcol);
    tout.row(i) = YmPIcol.elem(noNA).t() * U1.rows(noNA);
  }
  arma::vec out = arma::vectorise(tout);
  arma::uvec ig = arma::regspace<arma::uvec>(0, tout.size());
  arma::umat iG(ig.begin(), tout.n_rows, tout.n_cols, false);
  return Rcpp::List::create(Rcpp::Named("grad") = out,
                            Rcpp::Named("iG") = iG);
}

arma::mat d1PostZ_aCDM(arma::mat& Y, arma::mat& PI, arma::mat& Z,
                       arma::mat& G, arma::vec& mu, arma::mat& L){
  arma::mat PIm1PI = PI % (1-PI);
  arma::mat YmPI = (Y - PI)/(PIm1PI.clamp(arma::datum::eps, 1.0));
  const int p = Y.n_cols;
  const int n = Y.n_rows;
  const int q = mu.n_elem;
  arma::mat out(n,q);
  arma::mat Ri = L * L.t();
  arma::mat sM = arma::diagmat(arma::sqrt(arma::diagvec(Ri)));
  arma::cube dpidz = d1PIdZ_aCDM(G,Z,sM);
  for(int i = 0; i < p; i++){
    arma::mat temp = dpidz.slice(i).each_col() % YmPI.col(i);
    for(int j = 0; j < q; j++){
      for(int m = 0; m < n; m++){
        if(std::isnan(temp(m,j))) continue;
        out(m,j) += temp(m,j);
      }
    }
  }
  arma::mat Lt = L.t();
  arma::mat Zo = Z2Zo(Z, mu, Lt) ;
  out -= Zo;
  return(out);
}
