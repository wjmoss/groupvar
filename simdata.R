# Title     : simdata
# Objective : script for running simulations with given inputing parameters
# Created by: wj
# Created on: 2022/3/2

rm(list = ls())

## source files and required packages
source("generate.R")
source("greedysearch.R")

# "C:/Users/wj/Documents/R/win-library/4.0"
# "C:/Program Files/R/R-4.0.5/library"
.libPaths("C:/Users/wj/Documents/R/win-library/4.0")
.libPaths()
library(pcalg)
library(doParallel)


## on windows: Rscript simdata.R p=5 n=100 sp="s"
## on linux: R CMD BATCH --no-save '--args p=5 n=100 sp="s"' test.R test.out &
args <- commandArgs(trailingOnly = TRUE)
#p = args[1]
#n = args[2]
#sp = args[3]
#print(p)

if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
    p = 5
    n = 100
    sp = 's'
} else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

if (sp == 's'){
  prob = 3/(2*p-2)
} else if (sp == 'd'){
  prob = 0.3
} else{
  stop("Invalid arguments!")
}

#equivalent.eps = 1e-10
#verbose = FALSE
Oscale=1
faithful.eps=0
bic="std_bic"
#pop.version=FALSE (for constraint based algorithm, tbd)

## parameters for greedy search
replicate = 100 # number of replicate, for computing statistical properties
n.restarts = 5  # number of restarts in the greedy search algorithms, in each single replicate
mc.cores = 1   # number of cores in the parallelization
max.steps = Inf # maximum greedy search steps
max.in.degree = Inf # maximum in-degree (number of parents)
neighbourhood.size = 300 # the maximum size of searching neighborhood


## T: random ground truth, F: fix
randomgraph = T
## T: randomstart, F: true graph start
randomstart = T


## different seed can be fixed
seed <- 1
set.seed(seed)


## run simulations
res <- list()
for (r in 1:replicate){
  # >= p/3 blocks
  while (T){
    part <- sample(ceiling(p/3), p, replace = T)
    if (length(unique(part)) == ceiling(p/3))
      break()
  }

  ## generate ground truth
  if (randomgraph){
    gt <- GenerateGT(p, prob=prob, part, max.in.degree=2, Oscale, faithful.eps=faithful.eps)
  } else{
    # tbd: load some 0-1 matrix data?
    G <- matrix(0, p, p)
    tmp <- GenerateParams(mg=G, part=part, Oscale=1)
    gt <- list()
    gt$mg <- G
    gt$params$L <- tmp$L
    gt$params$O <- tmp$O
    gt$covmat <- GetSigma(tmp)
  }

  ## generate data, n*p, compute empirical covariance matrix
  data <- GenerateData(n, gt$params)
  covmat <- cov(data)

  res[[r]] <- list()
  res[[r]]$part <- part
  res[[r]]$graph <- gt
  res[[r]]$n.edges <- sum(gt$mg)
  res[[r]]$cpdag <- computeCPDAG_bk(gt$mg, part=part)
  if (randomstart){
    mg.start <- NULL
  } else {
    mg.start <- gt$mg
  }

  ## ges
  score <- new("GaussL0penObsScore", data, lambda =  0.5*log(nrow(data)))
  ges.fit <- ges(score)
  res[[r]]$ges.fit <- as(ges.fit$essgraph, "matrix") * 1


  ## ges-ev
  res[[r]]$ev_path <- greedySearch(
    cov.mat = covmat,
    mg.start = mg.start,
    part = rep(1,p),
    n = n,
    n.restarts = n.restarts,
    max.steps = max.steps,
    max.in.degree = max.in.degree,
    neighbourhood.size = neighbourhood.size,
    mc.cores = mc.cores
  )
  res[[r]]$gesev.fit <- res[[r]]$ev_path$final.mg

  ## ges-gev, save the search path
  res[[r]]$gev_path <- greedySearch(
    cov.mat = covmat,
    mg.start = mg.start,
    part = part,
    n = n,
    n.restarts = n.restarts,
    max.steps = max.steps,
    max.in.degree = max.in.degree,
    mc.cores = mc.cores
  )
  res[[r]]$gesgev.fit <- computeCPDAG_bk(res[[r]]$gev_path$final.mg, part)
  print(r)

  ## save loglikelihoods
  final.mg <- res[[r]]$gev_path$final.mg
  score <- res[[r]]$gev_path$final.score
  n.edges <- sum(final.mg)
  res[[r]]$loglike_gt <- sum( mvtnorm::dmvnorm(data, sigma=GetSigma(gt$params), log=TRUE) ) / n
  res[[r]]$loglike_est <- score + set.edge.penalty * (log(n)/2/n * (n.edges + p))

  # record best graph in each replicate?
  #if (i==50) {save(file = paste("tmp.p=", p, "n=", n, "seed1", ".RData"), res)}

}


## compute structural Hamming distances
err <- matrix(0,3, r, dimnames=list(c('gev','ev','nv')))
for (r in 1:replicate){
  err['gev', r] <- computeSHD(res[[r]]$gesgev.fit, res[[r]]$cpdag)
  err['ev', r] <- computeSHD(res[[r]]$gesev.fit, res[[r]]$cpdag)
  err['nv', r] <- computeSHD(res[[r]]$ges.fit, res[[r]]$cpdag)
  #err[4, r] <- shd(as(res[[r]]$gesgev.fit, "graphNEL"), as(res[[r]]$cpdag, "graphNEL"))
  #err[5, r] <- shd(as(res[[r]]$gesev.fit,"graphNEL"), as(res[[r]]$cpdag, "graphNEL"))
  #err[6, r] <- shd(as(res[[r]]$ges.fit,"graphNEL"), as(res[[r]]$cpdag, "graphNEL"))
}
err.avg <- apply(err, 1, mean)
#sapply(res, function(x) x$n.edges)

save(file = paste("p=", p, "n=", n, "sp=", sp, "seed=", seed, ".RData"), res, err, err.avg)