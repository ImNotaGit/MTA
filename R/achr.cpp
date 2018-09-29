#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List achr(const Rcpp::List& model, const Rcpp::List& state, const arma::mat& warmupPoints, const int nPoints, const int stepsPerPoint) {
    
    const double maxMinTol = 1e-09, uTol = 1e-09, dTol = 1e-14;
    arma::vec lb = Rcpp::as<arma::vec>(model["lb"]);
    arma::vec ub = Rcpp::as<arma::vec>(model["ub"]);
    arma::sp_mat S = Rcpp::as<arma::sp_mat>(model["S"]);
    arma::mat* Sd = new arma::mat(S);
    arma::mat N = arma::null(*Sd);
    delete Sd;
    unsigned int nRxns = S.n_cols;
    unsigned int nWarmupPoints = warmupPoints.n_cols;
    arma::vec centerPoint = Rcpp::as<arma::vec>(state["centr.pnt"]);
    arma::vec prevPoint = Rcpp::as<arma::vec>(state["prev.pnt"]);
    unsigned int totalStepCount = Rcpp::as<Rcpp::IntegerVector>(state["n.tot.steps"])[0];

    std::cout << "Start ACHR sampling." << std::endl;

    arma::mat points(nRxns, nPoints);
    arma::vec curPoint;
    arma::vec randVector;
    unsigned int pointCount = 0;
    while (pointCount < nPoints) {
        
        // create the random step size vector
        randVector = arma::randu(stepsPerPoint);
    
        unsigned int stepCount = 0;
        while (stepCount < stepsPerPoint) {
    
            // pick a random warmup point
            int randPointID = std::ceil(nWarmupPoints*arma::randu());
            arma::vec randPoint = warmupPoints.col(randPointID);
            // get a direction from the center point to the warmup point
            arma::vec u = randPoint - centerPoint;
            u = u / norm(u);
            // figure out the distances to upper and lower bounds
            arma::vec distUb = ub - prevPoint;
            arma::vec distLb = prevPoint - lb;
            // figure out if we are too close to a boundary
            arma::uvec validDir = arma::find((distUb > dTol) % (distLb > dTol)); // element-wise multiplication (%) as a workaround for &
            // figure out positive and negative directions
            arma::vec validU = u(validDir);
            arma::uvec posDirn = arma::find(validU > uTol);
            arma::uvec negDirn = arma::find(validU < -uTol);
            // figure out all the possible maximum and minimum step sizes
            arma::vec inverseValidU = 1 / validU;
            arma::vec maxStepTemp = distUb(validDir) % inverseValidU;
            arma::vec minStepTemp = -distLb(validDir) % inverseValidU;
            arma::vec maxStepVec = arma::join_vert(maxStepTemp(posDirn), minStepTemp(negDirn));
            arma::vec minStepVec = arma::join_vert(minStepTemp(posDirn), maxStepTemp(negDirn));
            // figure out the true max & min step sizes
            double maxStep = maxStepVec.min();
            double minStep = minStepVec.max();
            // find new direction if we're getting too close to a constraint
            if ((std::abs(minStep) < maxMinTol && std::abs(maxStep) < maxMinTol) || (minStep > maxStep)) continue;
            // pick a rand out of list_of_rands and use it to get a random step distance
            double stepDist = randVector(stepCount)*(maxStep-minStep)+minStep;
            // advance to the next point
            curPoint = prevPoint + stepDist*u;
            // reproject the current point
            if (totalStepCount % 10 == 0) {
                if (arma::abs(S*curPoint).max() > 1e-9) {
                  curPoint = N*(N.t()*curPoint);
                }
            }
            arma::uvec overInd = ub < curPoint;
            arma::uvec underInd = lb > curPoint;
            if (any(overInd) || any(underInd)) {
                arma::uvec overIndi = arma::find(overInd);
                arma::uvec underIndi = arma::find(underInd);
                curPoint(overIndi) = ub(overIndi);
                curPoint(underIndi) = lb(underIndi);
            }
            // recalculate the center point
            centerPoint = ((nWarmupPoints+totalStepCount)*centerPoint + curPoint)/(nWarmupPoints+totalStepCount+1);
            // next
            prevPoint = curPoint;
            totalStepCount = totalStepCount + 1;
            stepCount = stepCount + 1;
        }
    
        // add the current point to points
        points.col(pointCount) = curPoint;
        pointCount = pointCount + 1;
        // print progress
        std::cout << "Sampling progress: " << pointCount << " / " << nPoints << std::endl;
    }

    std::cout << "Done ACHR sampling." << std::endl;

    Rcpp::List endState = Rcpp::List::create(Rcpp::Named("centr.pnt") = centerPoint,
                                             Rcpp::Named("prev.pnt") = prevPoint,
                                             Rcpp::Named("n.tot.steps") = totalStepCount);
    return Rcpp::List::create(Rcpp::Named("sampl.pnts") = points,
                              Rcpp::Named("stat") = endState);
}