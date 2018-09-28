library(Matrix)
library(MASS)
library(Rcplex)
library(parallel)

sample.model <- function(model, params) {
  warmup.pnts <- sample.warmup.pnts(model, params$n.warmup)
  state <- achr.sampling(model, warmup.pnts, params$n.sampl.state, params$steps.per.pnt)$state
  achr.sampling(model, warmup.pnts, params$n.sampl, params$steps.per.pnt, state)$points
}

sample.warmup.pnts <- function(model, n) {
  orth.pnts <- get.orth.pnts(model, n)
  rand.pnts <- get.rand.pnts(model, n)
  n.rxns <- length(model$rxns)
  r <- rep(runif(n), each=n.rxns)
  dim(r) <- c(n.rxns, n)
  orth.pnts*r + rand.pnts*(1-r)
}

get.orth.pnts <- function(model, n) {
  n.rxns <- length(model$rxns)
  mat <- cbind(Diagonal(n.rxns), Diagonal(n.rxns, x=-1))
  if (n<=2*n.rxns) {
    mat <- mat[, sample(2*n.rxns, n)]
  } else {
    mat <- cbind(mat[, sample(2*n.rxns)], mat[, sample(2*n.rxns, n-2*n.rxns, replace=TRUE)])
  }
  apply(mat, 2, get.opt.pnt, model=model)
}

get.rand.pnts <- function(model, n) {
  n.rxns <- length(model$rxns)
  cs <- runif(n.rxns*n) - 0.5
  dim(cs) <- c(n.rxns, n)
  apply(cs, 2, get.opt.pnt, model=model)
}

get.opt.pnt <- function(model, c) {
  cvec <- c / norm(c,"2")
  objsense <- "max"
  Amat <- rbind(model$S, model$S)
  bvec <- c(model$rowlb, model$rowub)
  sense <- rep(c("G","L"), c(length(model$rowlb), length(model$rowub)))
  lb <- model$lb
  ub <- model$ub
  
  res <- Rcplex(cvec=cvec, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, control=list(maxcalls=Inf))
  if (res$status!=1) warning("Potential problems running LP. Solver status: ", res$status, ".\n")
  res$xopt
}

achr.sampling <- function(model, warmupPts, nPoints, stepsPerPoint, state=NULL) {
    
    N = Null(t(model$S))
    maxMinTol = 1e-09
    uTol = 1e-09
    dTol = 1e-14

    if (is.null(state)) {
      centerPoint = apply(warmupPts, 1, mean)
      totalStepCount = 0
      prevPoint = centerPoint
    } else {
      centerPoint = state$centerPoint
      totalStepCount = state$totalStepCount
      prevPoint = state$prevPoint
    }
    nRxns = ncol(model$S)
    W = ncol(warmupPts)
    nProj = 0
    totalCount = nPoints * stepsPerPoint
    Points = matrix(rep(0, nRxns * nPoints), ncol = nPoints)
    pointCount = 1
    while (pointCount <= nPoints) {
        randVector = runif(stepsPerPoint)
        stepCount = 1
        while (stepCount <= stepsPerPoint) {
            a = ceiling(W * runif(1))
            xa = warmupPts[, a]
            u = (xa - centerPoint)
            u = u/norm(as.matrix(u), "F")
            distUb = (uppbnd(model) - prevPoint)
            distLb = (prevPoint - lowbnd(model))
            validDir = ((distUb > dTol) & (distLb > dTol))
            posDirn = which(u[validDir] > uTol)
            negDirn = which(u[validDir] < -uTol)
            maxStepTemp = distUb[validDir]/u[validDir]
            minStepTemp = -distLb[validDir]/u[validDir]
            maxStepVec = c(maxStepTemp[posDirn], minStepTemp[negDirn])
            minStepVec = c(minStepTemp[posDirn], maxStepTemp[negDirn])
            maxStep = min(maxStepVec)
            minStep = max(minStepVec)
            print(sprintf("step %d: %f  %f", stepCount, minStep, 
                maxStep))
            if ((abs(minStep) < maxMinTol && abs(maxStep) < maxMinTol) || 
                (minStep > maxStep)) {
                print(sprintf("Warning %f %f\n", minStep, maxStep))
                next
            }
            stepDist = randVector[stepCount] * (maxStep - minStep) + 
                minStep
            curPoint = matrix(prevPoint + stepDist * u, ncol = 1)
            if ((totalStepCount%%stepsPerPoint) == 0) {
                print(sprintf("Sampling progress: %5.2f%%\n", 100*(stepCount+pointCount*stepsPerPoint)/(nPoints*stepsPerPoint)))
                if (max(abs(S(model) %*% curPoint)) > 1e-09) {
                  curPoint = N %*% (t(N) %*% curPoint)
                  nProj = nProj + 1
                }
            }
            curPoint[curPoint > uppbnd(model)] = uppbnd(model)[curPoint > 
                uppbnd(model)]
            curPoint[curPoint < lowbnd(model)] = lowbnd(model)[curPoint < 
                lowbnd(model)]
            prevPoint = curPoint
            stepCount = stepCount + 1
            totalStepCount = totalStepCount + 1
            centerPoint = ((W + totalStepCount) * centerPoint + 
                curPoint)/(W + totalStepCount + 1)
        }
        Points[, pointCount] = curPoint
        pointCount = pointCount + 1
    }

    list(state=list(centerPoint=centerPoint, totalStepCount=totalStepCount, prevPoint=prevPoint),
         points=Points)
}

