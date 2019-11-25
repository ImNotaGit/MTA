library(Rcplex.my)

CPX_STAT_CODE = c(
  `1`="CPX_STAT_OPTIMAL",
  `2`="CPX_STAT_UNBOUNDED",
  `3`="CPX_STAT_INFEASIBLE",
  `4`="CPX_STAT_INForUNBD",
  `5`="CPX_STAT_OPTIMAL_INFEAS",
  `6`="CPX_STAT_NUM_BEST",
  `10`="CPX_STAT_ABORT_IT_LIM",
  `11`="CPX_STAT_ABORT_TIME_LIM",
  `12`="CPX_STAT_ABORT_OBJ_LIM",
  `13`="CPX_STAT_ABORT_USER",
  `14`="CPX_STAT_FEASIBLE_RELAXED_SUM",
  `15`="CPX_STAT_OPTIMAL_RELAXED_SUM",
  `16`="CPX_STAT_FEASIBLE_RELAXED_INF",
  `17`="CPX_STAT_OPTIMAL_RELAXED_INF",
  `18`="CPX_STAT_FEASIBLE_RELAXED_QUAD",
  `19`="CPX_STAT_OPTIMAL_RELAXED_QUAD",
  `20`="CPX_STAT_OPTIMAL_FACE_UNBOUNDED",
  `21`="CPX_STAT_ABORT_PRIM_OBJ_LIM",
  `22`="CPX_STAT_ABORT_DUAL_OBJ_LIM",
  `23`="CPX_STAT_FEASIBLE",
  `24`="CPX_STAT_FIRSTORDER",
  `25`="CPX_STAT_ABORT_DETTIME_LIM",
  `30`="CPX_STAT_CONFLICT_FEASIBLE",
  `31`="CPX_STAT_CONFLICT_MINIMAL",
  `32`="CPX_STAT_CONFLICT_ABORT_CONTRADICTION",
  `33`="CPX_STAT_CONFLICT_ABORT_TIME_LIM",
  `34`="CPX_STAT_CONFLICT_ABORT_IT_LIM",
  `35`="CPX_STAT_CONFLICT_ABORT_NODE_LIM",
  `36`="CPX_STAT_CONFLICT_ABORT_OBJ_LIM",
  `37`="CPX_STAT_CONFLICT_ABORT_MEM_LIM",
  `38`="CPX_STAT_CONFLICT_ABORT_USER",
  `39`="CPX_STAT_CONFLICT_ABORT_DETTIME_LIM",
  `101`="CPXMIP_OPTIMAL",
  `102`="CPXMIP_OPTIMAL_TOL",
  `103`="CPXMIP_INFEASIBLE",
  `104`="CPXMIP_SOL_LIM",
  `105`="CPXMIP_NODE_LIM_FEAS",
  `106`="CPXMIP_NODE_LIM_INFEAS",
  `107`="CPXMIP_TIME_LIM_FEAS",
  `108`="CPXMIP_TIME_LIM_INFEAS",
  `109`="CPXMIP_FAIL_FEAS",
  `110`="CPXMIP_FAIL_INFEAS",
  `111`="CPXMIP_MEM_LIM_FEAS",
  `112`="CPXMIP_MEM_LIM_INFEAS",
  `113`="CPXMIP_ABORT_FEAS",
  `114`="CPXMIP_ABORT_INFEAS",
  `115`="CPXMIP_OPTIMAL_INFEAS",
  `116`="CPXMIP_FAIL_FEAS_NO_TREE",
  `117`="CPXMIP_FAIL_INFEAS_NO_TREE",
  `118`="CPXMIP_UNBOUNDED",
  `119`="CPXMIP_INForUNBD",
  `120`="CPXMIP_FEASIBLE_RELAXED_SUM",
  `121`="CPXMIP_OPTIMAL_RELAXED_SUM",
  `122`="CPXMIP_FEASIBLE_RELAXED_INF",
  `123`="CPXMIP_OPTIMAL_RELAXED_INF",
  `124`="CPXMIP_FEASIBLE_RELAXED_QUAD",
  `125`="CPXMIP_OPTIMAL_RELAXED_QUAD",
  `126`="CPXMIP_ABORT_RELAXED",
  `127`="CPXMIP_FEASIBLE",
  `128`="CPXMIP_POPULATESOL_LIM",
  `129`="CPXMIP_OPTIMAL_POPULATED",
  `130`="CPXMIP_OPTIMAL_POPULATED_TOL",
  `131`="CPXMIP_DETTIME_LIM_FEAS",
  `132`="CPXMIP_DETTIME_LIM_INFEAS"
)

rcplex <- function(cvec, Amat, bvec, Qmat=NULL, lb=0, ub=Inf, x0=NULL, control=list(), objsense=c("min", "max"), sense="L", vtype=NULL, n=1) {
  res <- Rcplex(cvec, Amat, bvec, Qmat, lb, ub, x0, control, objsense, sense, vtype, n)
  
  tmpf <- function(x) {
    res <- x[c("xopt", "obj")]
    res$stat <- x$status
    res$stat.str <- CPX_STAT_CODE[as.character(x$status)]
    if (res$stat.str %in% grep("LIM", CPX_STAT_CODE, value=TRUE) && res.stat!=104) {
      warning("in rcplex(): Reached time or other limit.", call.=FALSE)
    }
    if (res$stat %in% c(2,118)) { # unbounded
      warning("in rcplex(): Model unbounded.", call.=FALSE)
      if (objsense=="min") res$obj <- -Inf
      if (objsense=="max") res$obj <- Inf
    } else if (res$stat==20) {
      warning("in rcplex(): Model has an unbounded optimal face.", call.=FALSE)
    } else if (res$stat.str %in% grep("INFEAS", CPX_STAT_CODE, value=TRUE)) {
      warning("in rcplex(): Model may be infeasible, solver status: ", res$stat.str, call.=FALSE)
    } else if (res$stat %in% c(4,119)) { # infeasible or unbounded
      warning("in rcplex(): Model is infeasible or unbounded.", call.=FALSE)
    } else if (!res$stat %in% c(1,101,102,128,129,130)) {
      warning("in rcplex(): Solution may not be optimal or other issues occurred, solver status: ", res$stat.str, call.=FALSE)
    }
    res
  }
  
  if ("status" %in% names(res)) res <- list(res)
  res <- lapply(res, tmpf)
}

