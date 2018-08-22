source("/home/fountain/Documents/Projects/ourtools/MTA/Rwrapper/utils.R")

# specify the metabolic model to use, to initialize things
use.model("recon")

### 1. prepare expression data

# get data from GEO
tmp <- getGEO(GEO="GSE38718", GSEMatrix=TRUE, AnnotGPL=TRUE)
# check number of datasets
length(tmp)
# pre-process data: log-transformation, normalization etc.
m1 <- prep.data(tmp[[1]], norm.method="quantile")
# after examining the pData, create a data.table of the phenotypic variables of interest
m1.phe <- data.table(age=factor(m1[["age group:ch1"]], levels=c("Young","Older")), gender=m1[["gender:ch1"]])
# DE analysis
m1.de <- de(m1, m1.phe, coef="ageOlder")

# create source state discretized gene states, as input for iMAT
m1.expr.ref <- discrt.exprs.for.imat(m1[,m1.phe[,age=="Older"]])
# create direction of changes of reactions based on DE, as input for MTA
m1.dflux <- discrt.dflux.for.mta(m1.de, topn=100)

### 2. run iMAT in MATLAB

# start MATLAB server and setup environment
start.matlab()
init.matlab.for.mta()

# load model
load.model()
# run iMAT
m1.imat.res <- imat(m1.expr.ref, vout1="m1_v_ref", vout2="m1_sampl_pnts_ref", milp.solver="cplex", lp.solver="cplex")
# iMAT results:
# fluxes
m1.imat.res$v
# sample points
m1.imat.res$sampl.pnts

### 3. run MTA in MATLAB

# run MTA, continuing from the previous step
# del=NULL is default: screen all reaction deletions as well as the control (no deletion)
m1.mta.res <- mta(v.ref="m1_v_ref", dflux=m1.dflux, del=NULL, vout1="m1_mta_scores", vout2="m1_mta_stats", solver="cplex")
# or equivalently, pass v.ref also from R
m1.v.ref <- m1.imat.res$v
m1.mta.res <- mta(v.ref=m1.v.ref, dflux=m1.dflux, del=NULL, vout1="m1_mta_scores", vout2="m1_mta_stats", solver="cplex")
# MTA results: a data.table, del.rxn is the index of deleted reaction, 0 is the control; score is MTA score; stat is the returned status of the solver
m1.mta.res

# map the reactions to genes
m1.mta.res[, genes:=rxns2genes(del.rxn)]
# may need to check the returned solver status and proceed as appropriate. For now (2018.08.22), the value is dependent on which solver is used
m1.mta.res[, table(stat)]
# plot the distribution of MTA scores, and label the position of the control score
hist(m1.mta.res$score)
score.ctrl <- m1.mta.res[del.rxn==0, score]
abline(v=score.ctrl, col="red")
# those with scores > the control score
m1.mta.res.sub <- m1.mta.res[score>score.ctrl]
# *melt the above result to "long" form so that each gene mapped to the same reaction is assigned the score for that reaction, note that now it's possible that one gene has multiple scores
m1.mta.res.sub[, .(genes=unlist(genes), score), by=del.rxn][order(-score)]



