#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

double Phi(double x);
double phi(double x);
arma::mat mom(const arma::vec& xinf, const arma::vec& xsup);
arma::rowvec mom5d(double xinf, double xsup);

// [[Rcpp::export]]
Rcpp::List mep(const Rcpp::List& model, const double beta, const double damp, const unsigned int maxIter, const double dlb, const double dub, const double epsil, const bool ff, const unsigned int fIdx, double fMeans, double fVars) {
    
    const arma::mat S = arma::mat(Rcpp::as<arma::sp_mat>(model["S"]));
    arma::vec b = Rcpp::as<arma::vec>(model["b"]);
    arma::vec lb = Rcpp::as<arma::vec>(model["lb"]);
    arma::vec ub = Rcpp::as<arma::vec>(model["ub"]);
    const unsigned int nRxns = S.n_cols;
    const double factor = std::max(arma::abs(lb).max(), arma::abs(ub).max());
    b = b / factor;
    lb = lb / factor;
    ub = ub / factor;
    fMeans = fMeans / factor;
    fVars = fVars / (factor * factor);
    const arma::mat kk = beta * (S.t() * S);
    const arma::vec kb = beta * S.t() * b;
    arma::vec dlbv(nRxns);
    dlbv.fill(dlb);
    arma::vec dubv(nRxns);
    dubv.fill(dub);

    // variables to be returned
    arma::vec a = arma::zeros(nRxns);
    arma::vec d = arma::ones(nRxns);
    arma::vec av = arma::zeros(nRxns);
    arma::vec va = arma::ones(nRxns);
    arma::vec mu(nRxns);
    arma::vec s(nRxns);
    arma::mat D = arma::eye(nRxns, nRxns);

    std::cout << "Begin EP..." << std::endl;

    double err = 100;
    unsigned int iter = 0;

    while (err > epsil && iter < maxIter) {
        iter = iter + 1;
        // fast computation of the means and variances of the truncated Gaussians
        arma::mat I1 = arma::inv(kk + D);
        arma::mat v = I1 * (kb + D * a);
        arma::vec I = I1.diag();
        I = arma::min(I, d);
        arma::vec s1 = arma::min(dlbv, arma::max(dubv, 1/I - 1/d));
        s = 1 / s1;
        mu = (v - a%I/d) / (1 - I/d);
        arma::uvec ied = find(I == d);
        mu(ied) = 0.5 * (lb(ied) + ub(ied));
        // compute means and variances of the tilted distributions
        arma::vec s05 = arma::sqrt(s);
        arma::vec x0 = (lb - mu) / s05;
        arma::vec x1 = (ub - mu) / s05;
        arma::mat tmp = mom(x0, x1);
        arma::vec z = tmp.col(0);
        arma::vec eps = tmp.col(1);
        arma::vec oldva = va;
        arma::vec oldav = av;
        av = mu + z % s05;
        va = arma::max(arma::zeros(nRxns), s % (1 + eps));
        // fix known fluxes
        if (ff) {
            av(fIdx-1) = fMeans;
            va(fIdx-1) = fVars;
        }
        // moments matching (update "a" and "d")
        err = std::max(arma::abs(av-oldav).max(), arma::abs(va-oldva).max());
        arma::vec newd = 1/(1/va - 1/s);
        newd = arma::min(dubv, arma::max(dlbv, newd));
        arma::vec newa = av + (av - mu) % newd % s1;    
        a = damp * a + (1 - damp) * newa;
        d = damp * d + (1 - damp) * newd;
        D = arma::diagmat(1/d);
        if (ff) {
            err = std::max(err, std::abs(a(fIdx-1) - newa(fIdx-1)) + std::abs(d(fIdx-1) - newd(fIdx-1)));
        }
        // print progress
        if (iter % 100 == 0) {
            std::cout << "#iter = " << iter << "; error = " << err << "." << std::endl;
        }
    }

    std::cout << iter << " iterations before converged with error = " << err << "." << std::endl;
    std::cout << "Finished EP." << std::endl;

    mu = mu * factor;
    s = s * factor * factor;
    a = a * factor;
    d = d * factor * factor;
    av = av * factor;
    va = va * factor * factor;
    arma::mat cov = arma::inv(kk + D);

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("means.post") = mu,
                                           Rcpp::Named("vars.post") = s,
                                           Rcpp::Named("means.prior") = a,
                                           Rcpp::Named("vars.prior") = d,
                                           Rcpp::Named("means.trunc") = av,
                                           Rcpp::Named("vars.trunc") = va,
                                           Rcpp::Named("cov") = cov);
    return result;

}


double Phi(double x) {
    return 0.5 * (1 + std::erf(x/std::sqrt(2)));
}

double phi(double x) {
    return std::exp(-x*x/2) / std::sqrt(2*M_PI);
}

arma::mat mom(const arma::vec& xinf, const arma::vec& xsup) {
    const unsigned int n = xinf.n_elem;
    arma::mat res(n, 2);
    for (int i=0; i<n; i++) {
        res.row(i) = mom5d(xinf(i), xsup(i));
    }
    return res;
}

arma::rowvec mom5d(double xinf, double xsup) {
    arma::rowvec res(2);
    double scra1;
    double scra2;

    if (xsup - xinf < 1e-10) { 
        res(0) = 0.5*(xsup + xinf);
        res(1) = -1;
        return res;
    }
    
    if (std::min(std::abs(xinf), std::abs(xsup)) <= 6 || xinf*xsup <= 0) {
        double Phisup = Phi(xsup);
        double phisup = phi(xsup);
        double Phiinf = Phi(xinf);
        double phiinf = phi(xinf);
        scra1 = (phiinf - phisup)/(Phisup - Phiinf);
        scra2 = (xinf*phiinf - xsup*phisup)/(Phisup - Phiinf);
    } else {
        double delta2 = 0.5*(xsup*xsup - xinf*xinf);
        if (delta2 > 40) {
            scra1 = std::pow(xinf,5) / (3 - std::pow(xinf,2) + std::pow(xinf,4));
            scra2 = std::pow(xinf,6) / (3 - std::pow(xinf,2) + std::pow(xinf,4));
        } else {
            scra1 = pow(xinf*xsup,5) * (1 - std::exp(delta2)) / ( -std::exp(delta2)*(3.0 - std::pow(xinf,2) + std::pow(xinf,4))*std::pow(xsup,5) + std::pow(xinf,5)*(3 - std::pow(xsup,2) + std::pow(xsup,4)) );
            scra2 = pow(xinf*xsup,5) * (xsup - xinf*std::exp(delta2)) / ( -std::exp(delta2)*(3.0 - std::pow(xinf,2) + std::pow(xinf,4))*std::pow(xsup,5) + std::pow(xinf,5)*(3 - std::pow(xsup,2) + std::pow(xsup,4)) );
        }
    }

    res(0) = scra1;
    res(1) = scra2 - scra1*scra1;
    return res;
}

