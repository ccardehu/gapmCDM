#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR

#include <RcppArmadillo.h>
#include "utils.h"
#include "estim.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List gapmCDM_sim_rcpp(const int n, const int q, const int p,
                            arma::mat& A, arma::cube& C,
                            arma::vec& mu, arma::mat& R,
                            Rcpp::List& control){
  arma::vec knots = control["knots"];
  const unsigned int degree = control["degree"];
  const std::string basis = Rcpp::as<std::string>(control["basis"]);
  arma::mat Y(n,p);
  arma::mat Z = rmvNorm(n,mu,R);
  arma::mat U = Z2U(Z);
  arma::mat spM(n,q*C.n_cols);
  arma::cube spObj(n,q*C.n_cols,2);
  if(basis == "is"){
    spObj = SpU_isp(U,knots,degree);
    spM = spObj.slice(0);
  } else {
    spObj = SpU_bsp(U,knots,degree);
    spM = spObj.slice(0);
  }
  arma::mat PI = prob(A,C,spM);
  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      Y(i,j) = R::rbinom(1,PI(i,j));
    }
  }
  return Rcpp::List::create(Rcpp::Named("Y") = Y,
                            Rcpp::Named("Z") = Z,
                            Rcpp::Named("U") = U,
                            Rcpp::Named("spM") = spM,
                            Rcpp::Named("PI") = PI) ;
}

// [[Rcpp::export]]
Rcpp::List gapmCDM_fit_rcpp(arma::mat& Y, arma::mat& A, arma::cube& C, arma::cube& D,
                            arma::vec& mu, arma::mat& R, arma::mat& Z, Rcpp::List& control){

  const int burnin = control["burn.in"];
  const int iterlim = control["iter.lim"];
  const int tunelim = control["tune.lim"];
  const int window = control["window"];
  const unsigned int degree = control["degree"];
  const double epsStop = control["stop.eps"];
  const double gTune = control["tune.gamma"];
  const bool verbose = control["verbose"];
  const int every = control["verbose.every"];
  const bool stopFLAG = control["stop.atconv"];
  const bool traceFLAG = control["return.trace"];
  const bool corFLAG = control["cor.R"];
  const std::string sampler = Rcpp::as<std::string>(control["sampler"]);
  const std::string basis = Rcpp::as<std::string>(control["basis"]);
  const std::string algo = Rcpp::as<std::string>(control["algorithm"]);
  const double h = control["h"];
  const double g = control["gamma"];
  const double gA = control["gamma.A"];
  const double gC = control["gamma.CD"];
  const double gR = control["gamma.R"];
  const double sap = control["sa.power"];
  arma::vec knots = control["knots"];
  arma::mat Q = control["Qmatrix"];

  const int n = Y.n_rows;
  const int q = R.n_cols;
  const int np = knots.size() + degree ;

  if(q == 1){
    A.resize(Y.n_cols,q);
    A.fill(1.0);
  }

  arma::mat L = arma::chol(R,"lower");
  arma::mat Ln(L);
  arma::mat PI(arma::size(Y));
  arma::mat Zout(arma::size(Z)), Zn(Z);
  arma::mat Aout(arma::size(A)), An(A);
  arma::cube Cout(arma::size(C)), Cn(C);
  arma::cube Dout(arma::size(D)), Dn(D);
  arma::vec Mout(arma::size(mu)), Mn(mu);
  arma::mat Lout(arma::size(L));
  arma::mat Rout(arma::size(R));
  arma::mat spM(n,q*np), spD(n,q*np);
  arma::cube spObj(n,q*np,2);
  arma::cube pR(q,q,n), pRout(q,q,n);
  cube2eye(pR);
  arma::mat pM(arma::size(Z)), pMout(arma::size(Z));
  arma::uvec Rld;
  if(q > 1){
    Rld = arma::trimatl_ind(arma::size(R),-1);
  } else {
    Rld.resize(1);
    Rld.fill(0);
  }

  int tp = 0;
  tp = A.size() + C.size() + mu.size() + Rld.size();
  arma::mat theta(window,tp);
  arma::mat patrace(iterlim/10,tp);
  arma::mat lltrace(iterlim/10,2);
  arma::vec artrace(iterlim/10);
  const int etp = A.size() + C.size();
  const int rtp = q*(q+1)/2 - 1;
  arma::vec mt(etp,fill::zeros), vt(etp,fill::zeros);
  arma::vec mtr(rtp,fill::zeros), vtr(rtp,fill::zeros);
  int iter = 1;
  double ar = 0;
  double ssA;
  double ssC;
  double ssR;
  double ssZ;

  for(int ii = 1; ii <= iterlim + tunelim; ii++){
    if (ii % 2 == 0) Rcpp::checkUserInterrupt();
    arma::mat U = Z2U(Z);
    if(ii <= tunelim){
      ssA = gTune/n * gA;
      ssC = gTune/n * gC;
      ssR = gTune/n * gR;
      if(verbose & (ii % every == 0)) Rcpp::Rcout << "\r Iteration: " << std::setw(5) << ii << " (tuning phase)";
      ar = 0 ;
      ssZ = h ;
    } else {
      ssA = g/n * gA * std::pow(iter,sap) ;
      ssC = g/n * gC * std::pow(iter,sap) ;
      ssR = g/n * gR * std::pow(iter,sap) ;
      if(sampler == "ULA"){
        ssZ = h * std::pow(iter,-.333) ;
      } else if(sampler == "MALA"){
        ssZ = h * std::pow(q,-.333) ;
      } else if(sampler == "RWMH"){
        ssZ = h * std::pow(q,-1) ;
      }
    }
    if(basis == "is"){
      spObj = SpU_isp(U,knots,degree);
      spM = spObj.slice(0);
      spD = spObj.slice(1);
      PI = prob(A,C,spM);
      Rcpp::List gac = d1AC(Y,PI,spM,A,C);
      Rcpp::List ACn = newAC_MD(gac,A,Q,C,ssA,ssC);
      An = Rcpp::as<arma::mat>(ACn["A"]);
      Cn = Rcpp::as<arma::cube>(ACn["C"]);
    } else {
      spObj = SpU_bsp(U,knots,degree);
      spM = spObj.slice(0);
      spD = spObj.slice(1);
      PI = prob(A,C,spM);
      Rcpp::List gad = d1AD(Y,PI,spM,A,D);
      Rcpp::List ADn = newAD_MDp(gad,A,Q,D,ssA,ssC);
       An = Rcpp::as<arma::mat>(ADn["A"]);
      Dn = Rcpp::as<arma::cube>(ADn["D"]);
      Cn = D2C(Dn);
    }
    // arma::vec Mn = newM(mu,R,Z,ssAC);

    if(q > 1){
      if(algo == "GD"){
        Ln = newL(mu,L,Z,ssR,corFLAG);
      } else if(algo == "ADAM"){
        Ln = newL_adam(mu,L,Z,ssR,corFLAG,iter,mtr,vtr,control);
      } else if(algo == "mixed"){
        if(ii <= tunelim + 0.5*iterlim){
          Ln = newL(mu,L,Z,ssR,corFLAG);
        } else {
          Ln = newL_adam(mu,L,Z,ssR,corFLAG,iter,mtr,vtr,control);
        }
      }
    } else {
      Ln = L;
    }

    if(ii <= tunelim){
      Zn = Z;
    } else {
      if(sampler == "ULA"){
        Zn = newZ_ULA(Y,PI,Z,A,C,mu,L,spD,ssZ);
      } else if(sampler == "MALA"){
        Zn = newZ_MALA(Y,PI,Z,A,C,mu,L,ar,spD,ssZ,knots,degree,basis);
      } else if(sampler == "RWMH"){
        Zn = newZ_RWMH(Y,PI,Z,A,C,mu,L,ar,ssZ,knots,degree,basis);
      }
      newpMR(Zn,pM,pR,iter);
    }
    A = An;
    C = Cn;
    D = Dn;
    // mu = Mn;
    L = Ln;
    R = Ln*Ln.t();
    Z = Zn;

    for(int t1 = window - 1; t1 > 0; t1--){
      theta.row(t1) = theta.row(t1 - 1);
    }
    arma::rowvec ai = arma::vectorise(A).t();
    arma::rowvec ci = arma::vectorise(C).t();
    arma::rowvec mi = mu.t();
    arma::rowvec ri = R(Rld).t();
    arma::rowvec input = arma::join_rows(ai,ci,mi,ri);
    theta.row(0) = input;
    arma::mat dtheta = arma::diff(theta,1,0);
    double maxvalue = arma::max(arma::abs(arma::vectorise(dtheta)));

    if(iter % 10 == 0){
      patrace.row(iter/10-1) = input;
      double fzll = arma::accu(fz(Z,mu,R));
      lltrace(iter/10-1,0) = fzll;
      double fyzll = arma::accu(fyz(Y,PI));
      lltrace(iter/10-1,1) = fyzll;
      double cdll = fyzll + fzll;
      artrace(iter/10-1) = ar/(n*iter);

      if(verbose & (iter % every == 0)){
        if(iter >= burnin){
          if(sampler == "ULA"){
            Rcpp::Rcout << "\r Iteration: " << std::setw(5) << iter << " (estimation phase, CD-llk: " << std::to_string(cdll) <<
              ", max (abs) change: " << std::to_string(maxvalue) << ")";
          } else {
            Rcpp::Rcout << "\r Iteration: " << std::setw(5) << iter << " (estimation phase, CD-llk: " << std::to_string(cdll) <<
              ", max (abs) change: " << std::to_string(maxvalue) << ", AR (" << sampler << "): " << std::to_string(ar/(n*iter)) << ")";
          }
        } else {
          if(sampler == "ULA"){
            Rcpp::Rcout << "\r Iteration: " << std::setw(5) << iter << " (burn-in phase, CD-llk: " << std::to_string(cdll) <<
              ", max (abs) change: " << std::to_string(maxvalue) << ")";
          } else {
            Rcpp::Rcout << "\r Iteration: " << std::setw(5) << iter << " (burn-in phase, CD-llk: " << std::to_string(cdll) <<
              ", max (abs) change: " << std::to_string(maxvalue) << ", AR (" << sampler << "): " << std::to_string(ar/(n*iter)) << ")";
          }
        }
      }

      if(stopFLAG & (maxvalue < epsStop)) {
        iter++;
        break;
      }
    }

    if(iter >= burnin){
      Aout += A;
      Cout += C;
      Dout += D;
      Mout += mu;
      Lout += L;
      Rout += R;
      Zout += Z;
      pMout += pM;
      pRout += pR;
    }
    if(ii > tunelim) iter++;
  }
  Aout /= (iter - burnin);
  Cout /= (iter - burnin);
  Dout /= (iter - burnin);
  Mout /= (iter - burnin);
  Lout /= (iter - burnin);
  Rout /= (iter - burnin);
  Zout /= (iter - burnin);
  pMout /= (iter - burnin);
  pRout /= (iter - burnin);

  arma::mat Lout_ = arma::chol(Rout, "lower");
  if(corFLAG) Lout_ = ProxL(Lout_);

  const int nsim = control["nsim"];
  double llk(0), BIC(0), AIC(0);
  if(nsim != 0){
    llk = fy_gapmCDM_IS(Y,Aout,Cout,Mout,Lout_,Zout,
                        pMout,pRout,control);
    BIC = -2*llk + std::log(n)*Aout.n_elem;
    AIC = -2*llk + 2*Aout.n_elem;
  }
  arma::mat Uout = Z2U(Zout);
  spObj = SpU_bsp(Uout,knots,degree);
  spM = spObj.slice(0);
  PI = prob(Aout,Cout,spM);

  if(traceFLAG){
    return Rcpp::List::create(Rcpp::Named("A") = Aout,
                              Rcpp::Named("C") = Cout,
                              Rcpp::Named("D") = Dout,
                              Rcpp::Named("mu") = Mout,
                              Rcpp::Named("R") = Rout,
                              Rcpp::Named("Z") = Zout,
                              Rcpp::Named("U") = Uout,
                              Rcpp::Named("PI") = PI,
                              Rcpp::Named("llk") = llk,
                              Rcpp::Named("BIC") = BIC,
                              Rcpp::Named("AIC") = AIC,
                              Rcpp::Named("cdllk.trace") = lltrace,
                              Rcpp::Named("ar.trace") = artrace,
                              Rcpp::Named("theta.trace") = patrace,
                              Rcpp::Named("posMu") = pMout,
                              Rcpp::Named("posR") = pRout);
  } else {
    return Rcpp::List::create(Rcpp::Named("A") = Aout,
                              Rcpp::Named("C") = Cout,
                              Rcpp::Named("D") = Dout,
                              Rcpp::Named("mu") = Mout,
                              Rcpp::Named("R") = Rout,
                              Rcpp::Named("Z") = Zout,
                              Rcpp::Named("U") = Uout,
                              Rcpp::Named("PI") = PI,
                              Rcpp::Named("llk") = llk,
                              Rcpp::Named("BIC") = BIC,
                              Rcpp::Named("AIC") = AIC,
                              Rcpp::Named("posMu") = pMout,
                              Rcpp::Named("posR") = pRout);
  }
}


// [[Rcpp::export]]
Rcpp::List gapmCDM_cv_rcpp(arma::mat& Ytrain, arma::mat& Ytest, arma::mat& A, arma::cube& C,
                           arma::vec& mu, arma::mat& R, arma::mat& Z, Rcpp::List& control){

  const int nsim = control["nsim"];
  const unsigned int degree = control["degree"];
  const double h = control["h"];
  arma::vec knots = control["knots"];
  const bool verbose = control["verbose"];
  const int every = control["verbose.every"];
  const std::string sampler = Rcpp::as<std::string>(control["sampler"]);
  const std::string basis = Rcpp::as<std::string>(control["basis"]);

  const int n = Ytrain.n_rows;
  const int p = Ytrain.n_cols;
  const int q = R.n_cols;
  const int np = knots.size() + degree ;

  arma::mat Zn(arma::size(Z));
  arma::mat piH(n,p);
  arma::uvec isNA = arma::find_nan(Ytrain);
  const int nmis = isNA.n_elem;

  arma::mat spM(n,q*np), spD(n,q*np);
  arma::cube spObj(n,q*np,2);
  arma::mat L = arma::chol(R,"lower");

  double ar = 0;
  if(verbose) Rcpp::Rcout << "\n CV iteration: " ;

  for(int ii = 1; ii <= nsim; ii++){
    if (ii % 2 == 0) Rcpp::checkUserInterrupt();
    arma::mat U = Z2U(Z);
    if(basis == "is"){
      spObj = SpU_isp(U,knots,degree);
      spM = spObj.slice(0);
      spD = spObj.slice(1);
    } else {
      spObj = SpU_bsp(U,knots,degree);
      spM = spObj.slice(0);
      spD = spObj.slice(1);
    }
    arma::mat PI = prob(A,C,spM);
    if(sampler == "ULA"){
      double ssZ = h * std::pow(ii,-.33) ;
      Zn = newZ_ULA(Ytrain,PI,Z,A,C,mu,L,spD,ssZ);
      if(verbose & (ii % every == 0)) Rcpp::Rcout << "\r CV iteration: " << ii  << " ... ";
    } else if(sampler == "MALA"){
      double ssZ = h * std::pow(q,-.33) ;
      Zn = newZ_MALA(Ytrain,PI,Z,A,C,mu,L,ar,spD,ssZ,knots,degree,basis);
      if(verbose & (ii % every == 0)) Rcpp::Rcout << "\r CV iteration: " << ii  << ", AR (" << sampler << "): " << std::to_string(ar/(n*ii)) << " ... ";
    } else if(sampler == "RWMH"){
      double ssZ = h * std::pow(q,-1) ;
      Zn = newZ_RWMH(Ytrain,PI,Z,A,C,mu,L,ar,ssZ,knots,degree,basis);
      if(verbose & (ii % every == 0)) Rcpp::Rcout << "\r CV iteration: " << ii  << ", AR (" << sampler << "): " << std::to_string(ar/(n*ii)) << " ... ";
    }
    piH += PI;
    Z = Zn;
  }
  piH /= nsim;
  double out = 0;
  for(int ii = 0; ii < nmis; ii++){
    arma::uword id = isNA(ii);
    out += Ytest(id)*std::max(arma::datum::log_min, std::log(piH(id))) + (1-Ytest(id))*std::max(arma::datum::log_min, std::log(1-piH(id)));
  }
  out /= -nmis;
  if(verbose){
    if(sampler == "ULA"){
      Rcpp::Rcout << "\r CV iteration: " << nsim << " ... [CV-error: " << std::to_string(out) << "]\n" ;
    } else {
      Rcpp::Rcout << "\r CV iteration: " << nsim  << ", AR (" << sampler << "): " << std::to_string(ar/(n*nsim)) << " ... ";
    }
  }
  arma::vec outYhat = piH(isNA);
  arma::vec outYtst = Ytest(isNA);
  return Rcpp::List::create(Rcpp::Named("Yhat") = Rcpp::wrap(outYhat),
                            Rcpp::Named("Yobs") = Rcpp::wrap(outYtst),
                            Rcpp::Named("CV.error") = out);
}

// [[Rcpp::export]]
Rcpp::List apmCDM_sim_rcpp(const int n,
                           arma::mat& G,
                           arma::vec& mu, arma::mat& R){
  int p = G.n_rows;
  arma::mat Y(n,p);
  arma::mat Z = rmvNorm(n,mu,R);
  arma::mat U = Z2U(Z);
  arma::mat PI = prob_aCDM(G,U);
  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      Y(i,j) = R::rbinom(1,PI(i,j));
    }
  }
  return Rcpp::List::create(Rcpp::Named("Y") = Y,
                            Rcpp::Named("Z") = Z,
                            Rcpp::Named("U") = U,
                            Rcpp::Named("PI") = PI) ;
}

// [[Rcpp::export]]
Rcpp::List apmCDM_fit_rcpp(arma::mat& Y, arma::mat& G, arma::mat& Qmatrix,
                           arma::vec& mu, arma::mat& R, arma::mat& Z, Rcpp::List& control){

  const int burnin = control["burn.in"];
  const int iterlim = control["iter.lim"];
  const int tunelim = control["tune.lim"];
  const int window = control["window"];
  const double epsStop = control["stop.eps"];
  const double gTune = control["tune.gamma"];
  const bool verbose = control["verbose"];
  const int every = control["verbose.every"];
  const bool stopFLAG = control["stop.atconv"];
  const bool traceFLAG = control["return.trace"];
  const bool corFLAG = control["cor.R"];
  const bool slipFLAG = control["allow.slip"];
  const std::string sampler = Rcpp::as<std::string>(control["sampler"]);
  const std::string algo = Rcpp::as<std::string>(control["algorithm"]);
  const double h = control["h"];
  const double g = control["gamma"];
  const double damp = control["damp.factor"];
  double gG = control["gamma.G"];
  double gM = control["gamma.mu"];
  double gR = control["gamma.R"];
  const double sap = control["sa.power"];

  const int n = Y.n_rows;
  const int p = Y.n_cols;
  const int q = R.n_cols;

  arma::mat L = arma::chol(R,"lower");
  arma::mat Ln(L);
  arma::mat PI(arma::size(Y));
  arma::mat Zout(arma::size(Z)), Zn(Z);
  arma::mat Gout(arma::size(G)), Gn(arma::size(G));
  arma::vec Mout(arma::size(mu)), Mn(arma::size(mu));
  arma::mat Lout(arma::size(L));
  arma::mat Rout(arma::size(R));
  arma::cube pR(q,q,n), pRout(q,q,n);
  cube2eye(pR);
  arma::mat pM(arma::size(Z)), pMout(arma::size(Z));

  arma::uvec Rld;
  if(corFLAG){
    if(q > 1){
      Rld = arma::trimatl_ind(arma::size(R),-1);
    } else {
      Rld.resize(1);
      Rld.fill(0);
    }
  } else {
    if(q > 1){
      Rld = arma::trimatl_ind(arma::size(R),0);
    } else {
      Rld.resize(1);
      Rld.fill(0);
    }
  }

  arma::vec c1(p, arma::fill::value(1.0));
  arma::mat Qmatrix1 = join_rows(c1,Qmatrix);
  arma::uvec iG = arma::find(Qmatrix1 != 0);

  const int tp = iG.size() + mu.size() + Rld.size();
  arma::mat theta(window,tp);
  arma::mat patrace(iterlim/10,tp);
  arma::mat lltrace(iterlim/10,2);
  arma::vec artrace(iterlim/10);
  arma::vec mt(G.size(),fill::zeros), vt(G.size(),fill::zeros);
  const int rtp = Rld.size();
  arma::vec mtr(rtp,fill::zeros), vtr(rtp,fill::zeros);
  int iter = 1;
  double ar = 0;
  double ssG;
  double ssM;
  double ssR;
  double ssZ;

  for(int ii = 1; ii <= iterlim + tunelim; ii++){
    if (ii % 2 == 0) Rcpp::checkUserInterrupt();
    arma::mat U = Z2U(Z);
    if(ii <= tunelim){
      ssG = gTune/n * gG;
      ssM = gTune/n * gM;
      ssR = gTune/n * gR;
      if(verbose & (ii % every == 0)) Rcpp::Rcout << "\r Iteration: " << std::setw(5) << ii << " (tuning phase)";
      ar = 0 ;
      ssZ = h ;
    } else {
      ssG = g/n * gG * std::pow(iter,sap) ;
      ssM = g/n * gM * std::pow(iter,sap) ;
      ssR = g/n * gR * std::pow(iter,sap) ;
      if(sampler == "ULA"){
        ssZ = h * std::pow(iter,-.333) ;
      } else if(sampler == "MALA"){
        ssZ = h * std::pow(q,-.333) ;
      } else if(sampler == "RWMH"){
        ssZ = h * std::pow(q,-1) ;
      }
    }
    PI = prob_aCDM(G,U);
    Rcpp::List dG = d1G(Y,U,PI,G);
    if(slipFLAG){
      if(iter < 0.1*burnin){
        Gn = newG_MD(dG,G,ssG);
      } else {
        if(algo == "GD"){
          ssG *= damp;
          Gn = newG_PL2(dG,G,ssG);
        } else if(algo == "ADAM") {
          ssG *= damp;
          Gn = newG_PL2_adam(dG,G,ssG,iter,mt,vt,control);
        } else if(algo == "mixed") {
          if(ii <= tunelim + .5*iterlim){
            ssG *= damp;
            Gn = newG_PL2(dG,G,ssG);
          } else {
            ssG *= damp;
            Gn = newG_PL2_adam(dG,G,ssG,iter,mt,vt,control);
          }
        }
      }
    } else {
      Gn = newG_MD(dG,G,ssG);
    }

    Mn = newM(mu,R,Z,ssM);

    if(q > 1){
      if(algo == "GD"){
        Ln = newL(mu,L,Z,ssR,corFLAG);
      } else if(algo == "ADAM"){
        Ln = newL_adam(mu,L,Z,ssR,corFLAG,iter,mtr,vtr,control);
      } else if(algo == "mixed"){
        if(iter <= 0.5*iterlim){
          Ln = newL(mu,L,Z,ssR,corFLAG);
        } else {
          Ln = newL_adam(mu,L,Z,ssR,corFLAG,iter,mtr,vtr,control);
        }
      }
      if(iter < 0.1*burnin) Ln = ProxL(Ln);
    } else {
      Ln = L;
    }

    if(ii <= tunelim){
      Zn = Z;
    } else {
      if(sampler == "ULA"){
        Zn = newZ_ULA_aCDM(Y,PI,Z,G,mu,L,ssZ);
      } else if(sampler == "MALA"){
        Zn = newZ_MALA_aCDM(Y,PI,Z,G,mu,L,ssZ,ar);
      } else if(sampler == "RWMH"){
        Zn = newZ_RWMH_aCDM(Y,PI,Z,G,mu,L,ssZ,ar);
      }
      newpMR(Zn,pM,pR,iter);
    }

    G = Gn;
    mu = Mn;
    L = Ln;
    R = Ln*Ln.t();
    Z = Zn;

    for(int t1 = window - 1; t1 > 0; t1--){
      theta.row(t1) = theta.row(t1 - 1);
    }
    arma::vec Gpar = G(iG);
    arma::rowvec gi = Gpar.t();
    arma::rowvec mi = mu.t();
    arma::rowvec ri = R(Rld).t();
    arma::rowvec input = arma::join_rows(gi,mi,ri);
    theta.row(0) = input;
    arma::mat dtheta = arma::diff(theta,1,0);
    double maxvalue = arma::max(arma::abs(arma::vectorise(dtheta)));

    if(iter % 10 == 0){
      patrace.row(iter/10-1) = input;
      double fzll = arma::accu(fz(Z,mu,R));
      lltrace(iter/10-1,0) = fzll;
      double fyzll = arma::accu(fyz(Y,PI));
      lltrace(iter/10-1,1) = fyzll;
      double cdll = fyzll + fzll;
      artrace(iter/10-1) = ar/(n*iter);

      if(verbose & (iter % every == 0)){
        if(iter >= burnin){
          if(sampler == "ULA"){
            Rcpp::Rcout << "\r Iteration: " << std::setw(5) << iter << " (estimation phase, CD-llk: " << std::to_string(cdll) <<
              ", max (abs) change: " << std::to_string(maxvalue) << ")";
          } else {
            Rcpp::Rcout << "\r Iteration: " << std::setw(5) << iter << " (estimation phase, CD-llk: " << std::to_string(cdll) <<
              ", max (abs) change: " << std::to_string(maxvalue) << ", AR (" << sampler << "): " << std::to_string(ar/(n*iter)) << ")";
          }
        } else {
          if(sampler == "ULA"){
            Rcpp::Rcout << "\r Iteration: " << std::setw(5) << iter << " (burn-in phase, CD-llk: " << std::to_string(cdll) <<
              ", max (abs) change: " << std::to_string(maxvalue) << ")";
          } else {
            Rcpp::Rcout << "\r Iteration: " << std::setw(5) << iter << " (burn-in phase, CD-llk: " << std::to_string(cdll) <<
              ", max (abs) change: " << std::to_string(maxvalue) << ", AR (" << sampler << "): " << std::to_string(ar/(n*iter)) << ")";
          }
        }
      }

      if(stopFLAG & (maxvalue < epsStop)) {
        iter++;
        break;
      }
    }

    if(iter >= burnin){
      Gout += G;
      Mout += mu;
      Lout += L;
      Rout += R;
      Zout += Z;
      pMout += pM;
      pRout += pR;
    }
    if(ii > tunelim) iter++;
  }

  Gout /= (iter - burnin);
  Mout /= (iter - burnin);
  Lout /= (iter - burnin);
  Rout /= (iter - burnin);
  Zout /= (iter - burnin);
  pMout /= (iter - burnin);
  pRout /= (iter - burnin);

  arma::mat Lout_ = arma::chol(Rout, "lower");
  if(corFLAG) Lout_ = ProxL(Lout_);

  const int nsim = control["nsim"];
  double llk(0), BIC(0), AIC(0);
  if(nsim != 0){
    llk = fy_aCDM_IS(Y,Gout,Mout,Lout_,Zout,
                     pMout,pRout,control);
    BIC = -2*llk + std::log(n)*tp;
    AIC = -2*llk + 2*tp;
  }

  arma::mat Uout = Z2U(Zout);
  PI = prob_aCDM(Gout,Uout);

  if(traceFLAG){
    return Rcpp::List::create(Rcpp::Named("G") = Gout,
                              Rcpp::Named("mu") = Mout,
                              Rcpp::Named("R") = Rout,
                              Rcpp::Named("Z") = Zout,
                              Rcpp::Named("U") = Uout,
                              Rcpp::Named("PI") = PI,
                              Rcpp::Named("llk") = llk,
                              Rcpp::Named("BIC") = BIC,
                              Rcpp::Named("AIC") = AIC,
                              Rcpp::Named("cdllk.trace") = lltrace,
                              Rcpp::Named("ar.trace") = artrace,
                              Rcpp::Named("theta.trace") = patrace,
                              Rcpp::Named("posMu") = pMout,
                              Rcpp::Named("posR") = pRout);
  } else {
    return Rcpp::List::create(Rcpp::Named("G") = Gout,
                              Rcpp::Named("mu") = Mout,
                              Rcpp::Named("R") = Rout,
                              Rcpp::Named("Z") = Zout,
                              Rcpp::Named("U") = Uout,
                              Rcpp::Named("PI") = PI,
                              Rcpp::Named("llk") = llk,
                              Rcpp::Named("BIC") = BIC,
                              Rcpp::Named("AIC") = AIC,
                              Rcpp::Named("posMu") = pMout,
                              Rcpp::Named("posR") = pRout);
  }
}

// [[Rcpp::export]]
Rcpp::List apmCDM_cv_rcpp(arma::mat& Ytrain, arma::mat& Ytest, arma::mat& G,
                          arma::vec& mu, arma::mat& R, arma::mat& Z, Rcpp::List& control){

  const int nsim = control["nsim"];
  const double h = control["h"];
  const bool verbose = control["verbose"];
  const int every = control["verbose.every"];
  const std::string sampler = Rcpp::as<std::string>(control["sampler"]);

  const int n = Ytrain.n_rows;
  const int p = Ytrain.n_cols;
  const int q = mu.n_elem;

  arma::mat Zn(arma::size(Z));
  arma::mat piH(n,p);
  arma::uvec isNA = arma::find_nan(Ytrain);
  const int nmis = isNA.n_elem;

  double ar = 0;
  if(verbose) Rcpp::Rcout << "\n CV iteration: " ;

  for(int ii = 1; ii <= nsim; ii++){
    if (ii % 2 == 0) Rcpp::checkUserInterrupt();
    arma::mat U = Z2U(Z);
    arma::mat PI = prob_aCDM(G,U);
    if(sampler == "ULA"){
      double ssZ = h * std::pow(ii,-.33) ;
      Zn = newZ_ULA_aCDM(Ytrain,PI,Z,G,mu,R,ssZ);
      if(verbose & (ii % every == 0)) Rcpp::Rcout << "\r CV iteration: " << ii  << " ... ";
    } else if(sampler == "MALA"){
      double ssZ = h * std::pow(q,-.33) ;
      Zn = newZ_MALA_aCDM(Ytrain,PI,Z,G,mu,R,ssZ,ar);
      if(verbose & (ii % every == 0)) Rcpp::Rcout << "\r CV iteration: " << ii  << ", AR (" << sampler << "): " << std::to_string(ar/(n*ii)) << " ... ";
    } else if(sampler == "RWMH"){
      double ssZ = h * std::pow(q,-1) ;
      Zn = newZ_RWMH_aCDM(Ytrain,PI,Z,G,mu,R,ssZ,ar);
      if(verbose & (ii % every == 0)) Rcpp::Rcout << "\r CV iteration: " << ii  << ", AR (" << sampler << "): " << std::to_string(ar/(n*ii)) << " ... ";
    }
    piH += PI;
    Z = Zn;
  }
  piH /= nsim;
  double out = 0;
  for(int ii = 0; ii < nmis; ii++){
    arma::uword id = isNA(ii);
    out += Ytest(id)*std::max(arma::datum::log_min, std::log(piH(id))) + (1-Ytest(id))*std::max(arma::datum::log_min, std::log(1-piH(id)));
  }
  out /= -nmis;
  if(verbose){
    if(sampler == "ULA"){
      Rcpp::Rcout << "\r CV iteration: " << nsim << " ... [CV-error: " << std::to_string(out) << "]\n" ;
    } else {
      Rcpp::Rcout << "\r CV iteration: " << nsim  << ", AR (" << sampler << "): " << std::to_string(ar/(n*nsim)) << " ... ";
    }
  }
  arma::vec outYhat = piH(isNA);
  arma::vec outYtst = Ytest(isNA);
  return Rcpp::List::create(Rcpp::Named("Yhat") = Rcpp::wrap(outYhat),
                            Rcpp::Named("Yobs") = Rcpp::wrap(outYtst),
                            Rcpp::Named("CV.error") = out);
}
