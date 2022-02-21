source("generate.R")
source("greedysearch.R")
source("plotting.R")


## parameters for randomly generating
p=6
part=c(1,1,2,2,3,3)
n=1000
mc.cores = 1
max.steps = 100
equivalent.eps = 1e-10
verbose = FALSE
Bdist="sunif"
Oscale=1
faithful.eps=0
pop.version=FALSE
bic="std_bic"


## parameters for greedy search
replicate = 10 # number of replicate, for computing statistical properties
R = 99  # number of restarts in the greedy search algorithms, in each single replicate
M = 1   # number of cores in the parallelization
maxSteps = Inf # maximum greedy search steps
max.in.degree = 2


## T: random, F: fix
randomgraph = T
## T: randomstart, F: true graph start
randomstart = T

res <- list()

#different seed can be set
seed <- 1
set.seed(seed)


for (r in 1:replicate){

  ## generate graph
  if (randomgraph){
    gt <- GenerateGT(p, part, max.in.degree=1, Oscale, faithful.eps=faithful.eps)
  } else{
    G <- matrix(0, p, p)
    tmp <- GenerateParams(mg=G, part=part, Oscale=1)
    gt <- list()
    gt$mg <- G
    gt$params$B <- tmp$B
    gt$params$Omega <- tmp$Omega
    gt$covMat <- GetSigma(tmp)
  }

  ## generate data
  data <- GenerateData(n, gt$params)
  covMat <- cov(data)

  ## greedy search, save the search path
  res[[r]] <- list()
  res[[r]]$graph <- gt
  if (randomstart){
    mg.start <- NULL
  } else{
    mg.start <- gt$mg
  }
  res[[r]]$path <- greedySearch(
    cov.mat = covMat,
    mg.start = mg.start,
    part = part,
    #gt$mg
    n = n,
    n.restarts = R,
    max.steps = maxSteps,
    max.in.degree = max.in.degree,
    mc.cores = M
  )
  print(r)

  ## save loglikelihoods
  final.mg <- res[[r]]$path$final.mg
  score <- res[[r]]$path$final.score
  n.edges <- sum(final.mg)
  res[[r]]$loglike_gt <- sum( mvtnorm::dmvnorm(data, sigma=GetSigma(gt$params), log=TRUE) ) / n
  if (bic == "ext_bic"){
    res[[r]]$loglike_est <- score + set.edge.penalty * (log(n)/2/n * (n.edges + p)) - (n.edges + p) * (2*log(p) + log(3)) / n
  }
  else{
    res[[r]]$loglike_est <- score + set.edge.penalty * (log(n)/2/n * (n.edges + p))
  }

  # record best graph in each replicate?
  #if (i==50) {save(file = paste("tmp.p=", p, "n=", n, "seed1", ".RData"), res)}

}

#save(file = paste("p=", p, "n=", n, "seed=", seed, ".RData"), res)
