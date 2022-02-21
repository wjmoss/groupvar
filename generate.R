#generate a simple mixed graph uniformly, by MCMC sampling
GenerateGraph <- function(p, N=1, iter=p^4, p1=1, max.in.degree=Inf, names=paste("V", 1:p, sep="")){
  # p1: P(directed | empty) = P(empty | directed)
  # 1-p1: P(bidirected | empty) = P(empty | bidirected)
  # p2: P(bidirected | directed) = P(directed | bidirected)
  # p3: P(direction switch | directed)

  # p2 must be < min(p1, 1-p1)
  # p3 must be < 1-p1-p2

  # only directed edges: p1=1
  mg <- matrix(0, p, p)
  rownames(mg) <- names
  colnames(mg) <- names
  index <- 1
  res <- list()

  for (k in 1:(N * iter)){

    #sample a position
    pos <- sample(p*p, 1)
    j <- ceiling(pos / p)
    i <- (pos - 1) %% p + 1
    #indegree_i <- (length(which(mg[,i]>0)) < max.in.degree)
    indegree_j <- (length(which(mg[,j]>0)) < max.in.degree)

    # the condition of if {in the bracket "if ()"} must be True of False,
    # NA and other value lead to returning value NA, and hence it raises error
    if (i != j && mg[i,j] == 0 && mg[j,i] == 0){
      #no edge in [i,j]
      if (indegree_j){
        mg[i,j] <- 1
      }
    }
      # an edge
    else if (mg[i,j] != 0){
      # remove
      mg[i,j] <- 0
    }
    #print(sum(mg))

    #take the sample when k%%iter == 0
    if (k >= iter && k %% iter == 0){
      res[[index]] <- mg
      index <- index + 1
    }

  }


  return(res)
}


GenerateParams <- function(Linit=NULL, mg=NULL, part=NULL, L.lb=0.5, L.ub=0.9,
                           O.lb=0.3, O.ub=1, Oscale=1, paramsign='posneg'){
  if (is.null(Linit)){
    Linit <- mg
    Linit[Linit==100] <- 0
  }

  indL <- which(Linit != 0)
  if(paramsign == "pos") {Lvals <- runif(length(indL), L.lb, L.ub)}
  else {Lvals <- (2*rbinom(length(indL),1,0.5)-1)*runif(length(indL), -L.ub, -L.lb)}
  L <- Linit
  L[indL] <- Lvals


  O <- rep(0, p)
  for (k in 1:max(part)){
    ind <- which(part == k)
    O[ind] <- runif(1, O.lb, O.ub)
  }

  #rescale source nodes
  #sourcenodes <- all(Linit[,j]==0)
  #sourcenodes <- which(sapply(1:ncol(L), function(j) sum(Linit[,j]==0)))
  #O[,sourcenodes] <- O[,sourcenodes] * sqrt(Oscale)
  #O[sourcenodes,] <- O[sourcenodes,] * sqrt(Oscale)

  return (list(L=L, O=O))
}


GetSigma <- function(params){
  p <- ncol(params$L)
  # b %*% a = b * a
  return (solve(diag(p)-t(params$L)) %*% (params$O * solve(diag(p)-params$L)))
}


#isFaithful??
isFaithful <- function(mg, Lhat, Ohat, Shat, faithful.eps) {
  L.ind <- which((abs(Lhat) < faithful.eps) & (mg == 1))
  O.ind <- which((abs(Ohat) < faithful.eps) & (mg == 100))
  S.ind <- which((abs(Shat) < faithful.eps) & (mg > 0))
  if ((length(L.ind > 0)) || (length(O.ind > 0)) || (length(S.ind > 0))) return(list(flag=FALSE, L.ind=L.ind, O.ind=O.ind, S.ind=S.ind))
  return(list(flag=TRUE))
}


# generate ground truth
GenerateGT <- function(p, part, max.in.degree=Inf, Oscale=1, faithful.eps=0, paramsign = "posneg"){
  res <- list()
  res$mg <- GenerateGraph(p, N=1, max.in.degree=max.in.degree)[[1]]
  res$params <- GenerateParams(mg=res$mg, part=part, Oscale=Oscale, paramsign = paramsign)
  while (! isFaithful(res$mg, t(res$params$L), res$params$O, GetSigma(res$params), 10*faithful.eps)$flag) {
    print("Ground Truth not faithgful - regenerating...")
    res$params <- GenerateParams(mg=res$mg, Oscale=Oscale, paramsign = paramsign)
  }
  res$covmat <- GetSigma(res$params)
  return (res)
}


GenerateData <- function(n, params){
  err <- mvtnorm::rmvnorm(n, sigma=diag(params$O))
  data <- err %*% solve(diag(ncol(params$L))-params$L)
  #colnames(data) <- if (is.null(colnames(params$L))) 1:ncol(params$L) else colnames(params$L)
  return (data)
}

