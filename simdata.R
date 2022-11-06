# Title     : simdata
# Objective : script for running simulations with given input parameters
# Created by: wj
# Created on: 2022/11/6

rm(list = ls())

## source files and required packages
source("generate.R")
source("greedysearch.R")
source("utils.R")

.libPaths("/home/wuju/linux/R/x86_64-pc-linux-gnu-library/4.1")
library(pcalg)
library(doParallel)


## on windows: Rscript simdata.R p=5 n=100 sp="s"
## on linux: R CMD BATCH --no-save '--args p=5 n=100 sp="s"' simdata.R simdata.out &
args <- commandArgs(trailingOnly = TRUE)
p = as.numeric(args[1])
n = as.numeric(args[2])
sp = as.character(args[3])


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
randomstart = F


## different seed can be fixed
seed <- 1
set.seed(seed)


## run simulations
res <- list()
for (r in 1:replicate){
  # >= p/3+1 blocks
  while (T){
    part <- sample(ceiling(p/3)+1, p, replace = T)
    if (length(unique(part)) == ceiling(p/3)+1)
      break()
  }

  ## generate ground truth
  if (randomgraph){
    gt <- GenerateGT(p, prob=prob, part, Oscale, faithful.eps=faithful.eps)
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

  score <- new("GaussL0penObsScore", data, lambda =  0.5*log(nrow(data)) )

  ## ges
  start.time <- Sys.time()
  ges.fit <- ges(score)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  res[[r]]$ges.fit <- as(ges.fit$essgraph, "matrix") * 1
  res[[r]]$ges.time <- time.taken

  ## pc
  alpha <- 1e-5 # Compute the path, starting with this value of alpha
  alpha.base <- 1.1 # The next value of alpha is alpha * 1.1
  path <- list()
  counter <- 1
  while(T) {
    #suffstat <- list(C = cov2cor(covmat), n = n)
    suffstat <- list(C = cor(data), n = n)
    start.time <- Sys.time()
    pc.fit <- pc(suffStat = suffstat, indepTest = gaussCItest, alpha = alpha, p = p,
                 u2pd='retry', verbose=F, skel.method="stable")
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    #if (!ggm::isAcyclic(as(pc.fit, 'amat')) ){
      #Rgraphviz::plot(pc.fit@graph)
      #stop("Cyclic graph!")
    #}


    u.amat <- as(pc.fit@graph, "matrix") * 1
    if(sum(u.amat) == 0) {
      alpha <- alpha * alpha.base
      next()
    }

    est.cpdag <- as(pc.fit@graph, "matrix") * 1
    onedag <- cpdag_repr(est.cpdag)

    # Record these results
    path[[counter]] <- list()
    path[[counter]]$alpha <- alpha
    path[[counter]]$graph <- est.cpdag
    path[[counter]]$score <- score$global.score(as(onedag, 'GaussParDAG'))
    path[[counter]]$fitting.time <- time.taken

    counter <- counter + 1
    alpha <- alpha * alpha.base
    if (alpha > 0.8) {
      break()
    }

  }

  best <- which.max(sapply(1:length(path), function(i) path[[i]]$score))
  best.alpha <- path[[best]]$alpha
  best.score <- path[[best]]$score
  res[[r]]$pc.fit <- path[[best]]$graph
  res[[r]]$pc.time <- path[[best]]$fitting.time
  res[[r]]$pc.cumtime <- sum(sapply(1:length(path), function(i) path[[i]]$fitting.time))

  ## ges-ev
  #res[[r]]$ev_path <- greedySearch(
    #cov.mat = covmat,
    #mg.start = mg.start,
    #part = rep(1,p),
    #n = n,
    #n.restarts = n.restarts,
    #max.steps = max.steps,
    #max.in.degree = max.in.degree,
    #neighbourhood.size = neighbourhood.size,
    #mc.cores = mc.cores
  #)
  #res[[r]]$gesev.fit <- res[[r]]$ev_path$final.mg

  ## ges-gev, save the search path
  start.time <- Sys.time()
  #suffstat <- list(C = cor(data), n = n)
  #pc.fit <- pc(suffStat = suffstat, indepTest = gaussCItest, alpha = 0.1, p = p,
               #u2pd='retry', verbose=F, skel.method="stable")
  #est.cpdag <- as(pc.fit@graph, "matrix") * 1
  #onedag <- cpdag_repr(est.cpdag)
  res[[r]]$gev_path <- greedySearch(
    cov.mat = covmat,
    mg.start = list(matrix(0,p,p)),
    #mg.start = list(onedag),
    part = part,
    n = n,
    n.restarts = n.restarts,
    max.steps = max.steps,
    max.in.degree = max.in.degree,
    neighbourhood.size = neighbourhood.size,
    mc.cores = mc.cores
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  res[[r]]$gev.fit <- computeCPDAG_bk(res[[r]]$gev_path$final.mg, part)
  res[[r]]$gev.time <- time.taken
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
err1 <- matrix(0,3, replicate, dimnames=list(c('ges_gev','pc','ges_nv')))
for (r in 1:replicate){
  err1['ges_gev', r] <- computeSHD(res[[r]]$gev.fit, res[[r]]$cpdag)
  err1['pc', r] <- computeSHD(res[[r]]$pc.fit, res[[r]]$cpdag)
  err1['ges_nv', r] <- computeSHD(res[[r]]$ges.fit, res[[r]]$cpdag)
  #err[4, r] <- shd(as(res[[r]]$gesgev.fit, "graphNEL"), as(res[[r]]$cpdag, "graphNEL"))
  #err[5, r] <- shd(as(res[[r]]$gesev.fit,"graphNEL"), as(res[[r]]$cpdag, "graphNEL"))
  #err[6, r] <- shd(as(res[[r]]$ges.fit,"graphNEL"), as(res[[r]]$cpdag, "graphNEL"))
}
err1.avg <- apply(err1, 1, mean)
#sapply(res, function(x) x$n.edges)

err2 <- matrix(0,3, replicate, dimnames=list(c('ges_gev','pc','ges_nv')))
for (r in 1:replicate){
  err2['ges_gev', r] <- computeSHD(res[[r]]$gev.fit, res[[r]]$cpdag, reversal=2)
  err2['pc', r] <- computeSHD(res[[r]]$pc.fit, res[[r]]$cpdag, reversal=2)
  err2['ges_nv', r] <- computeSHD(res[[r]]$ges.fit, res[[r]]$cpdag, reversal=2)
  #err[4, r] <- shd(as(res[[r]]$gesgev.fit, "graphNEL"), as(res[[r]]$cpdag, "graphNEL"))
  #err[5, r] <- shd(as(res[[r]]$gesev.fit,"graphNEL"), as(res[[r]]$cpdag, "graphNEL"))
  #err[6, r] <- shd(as(res[[r]]$ges.fit,"graphNEL"), as(res[[r]]$cpdag, "graphNEL"))
}
err2.avg <- apply(err2, 1, mean)

#mean(sapply(1:replicate, function(r) res[[r]]$ges.time))
#mean(sapply(1:replicate, function(r) res[[r]]$pc.cumtime))
#mean(sapply(1:replicate, function(r) res[[r]]$gev.time))

# time
time.taken <- matrix(0, 3, replicate, dimnames=list(c('ges_gev','pc','ges_nv')))
for (r in 1:replicate){
  time.taken['ges_gev', ] <- sapply(1:replicate, function(r) res[[r]]$gev.time)
  time.taken['pc', ] <- sapply(1:replicate, function(r) res[[r]]$pc.cumtime)
  time.taken['ges_nv', ] <- sapply(1:replicate, function(r) res[[r]]$ges.time)
}
time.avg <- apply(time.taken, 1, mean)
time.max <- apply(time.taken, 1, max)


Sys.time()
fnm <- paste0("./empty_stdbic/p=", p, "n=", n, "sp=", sp, "seed=", seed, blocks, ".RData")
save(file=fnm,res, err1, err1.avg, err2, err2.avg)

