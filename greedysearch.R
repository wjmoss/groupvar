# default setups
set.edge.penalty = 1
bic = "std_bic"


#not needed?
computeOmega <- function(mg, part, k, covMat){
  # mg: zero pattern of directed graph
  # part: a vector indicating the partition of nodes
  # k: the interesting partition block
  # covMat: covariance matrix from data

  ind <- which(part == k)
  omega <- 0
  for (i in ind){
    pa <- which(mg[, i] != 0)
    alpha <- solve(covMat[pa, pa], covMat[pa, i])
    err <- err + covMat[i,i] - t(alpha) %*% covMat[pa, i]
  }
  omega <- omega / length(ind)
  return (omega)
}


partScore <- function(mg, part, k, n, covMat, edge.penalty=set.edge.penalty, bic = "std_bic"){
  # mg: zero pattern of directed graph
  # part: a vector indicating the partition of nodes
  # k: the interesting partition block
  # n: sample size
  # covMat: sample covariance matrix
  # edge.penalty: hyperparameter lambda
  # bic: std_bic or ext_bic

  p <- ncol(mg)
  ind <- which(part == k)

  omega <- 0
  for (i in ind){
    pa <- which(mg[, i] != 0)
    if (length(pa) > 0){
      alpha <- solve(covMat[pa, pa], covMat[pa, i])
      omega <- omega + covMat[i,i] - t(alpha) %*% covMat[pa, i]
    }
    else{
      omega <- omega + covMat[i,i]
    }
  }
  omega <- omega / length(ind)

  n.edges <- sum(mg[, ind] != 0)
  score <- -0.5 * (length(ind) * omega + length(ind)) # second part is a cst, omit?

  if (bic == 'std_bic'){
    score <- score - edge.penalty * (log(n)/2/n * (n.edges + p))
  }
  else {
    stop("Other penalty is not implemented yet!")
  }

  return (score)
}


getSkeleton <- function(mg){
  sk <- (mg+t(mg)>0)*1
  return (sk)
}


#need to be modified
getVStructures <- function(mg){
  res <- list()
  index <- 1

  for (j in 1:ncol(mg)){
    for (i in 1:(ncol(mg)-1)){
      for (k in (i+1):ncol(mg)){
        #print(c(i,j,k))
        if (i != j && k!= j && mg[i,j] != 0 && mg[k,j] != 0){
          res[[index]] <- c(i,j,k)
            index <- index + 1
        }
      }
    }
  }
}


fastGreedySearch <- function(mg.start, part, n, covMat, maxSteps=Inf, max.in.degree=Inf, direction=3,
                             edge.penalty=1, verbose=TRUE, faithful.eps=0, max.pos=Inf, dags.only=FALSE,
                             eps.conv=1e-12, bic = "std_bic")
{
  # mg.start: simple mixed graph, with entries 1/100
  # part: a vector indicating the partition of nodes
  # n: sample size
  # covMat: sample covariance matrix
  # maxSteps: maximal number of greedy search steps
  # max.in.degree: maximum in-degree (size of parents)
  # direction: control forward/backward search, edge reversal
  # edge.penalty: hyperparameter lambda
  # verbose: print info
  # faithful.eps: faithful eps in ricf
  # max.pos: maximal number of positions for adding/deleting/reversing
  # eps.conv: eps for greedy search
  #W: restarts of ricf
  #bic: currently only "std_bic"

  t <- proc.time()[3]
  p <- ncol(mg.start)
  if (is.null(colnames(mg.start))) {
    nnames <- 1:ncol(mg.start)
  } else {
    nnames <- colnames(mg.start)
  }

  #initial score
  # Compute initial score
  #score <- computeScore(mg.start, covMat, n, maxIter, edge.penalty, faithful.eps=faithful.eps, W=W, bic = bic)
  scores <- c()
  for (k in 1:max(part)){
    scores[k] <- partScore(mg.start, part, k, n, covMat, edge.penalty, bic = bic)
  }
  score <- sum(scores)

  #greedy search
  iter <- 1

  #state
  state <- list()
  state$mg <- mg.start
  state$scores <- scores
  state$score <- score
  states <- list()

  # Take time for first "step"
  state$t <- proc.time()[3] - t
  state$ct <- state$t  # Cumulative time
  t <- proc.time()[3]


  while(iter < maxSteps){
    # 1 -- only forward
    # 2 -- only backward
    # 3 -- both
    # <=3 do edge reversal

    #forward
    cand.add <- list()
    i.a <- 0
    if (direction != 2){
      poslist <- which(state$mg + t(state$mg) == 0 & lower.tri(state$mg))
      # if length == 1, resample command sample over 1:pos
      if (length(poslist) > 1) poslist <- sample(poslist, min(length(poslist), max.pos))
      for (pos in poslist){
        #(i,j) entry
        j <- ceiling(pos / p)
        i <- (pos - 1) %% p + 1
        #if (i == j) {
          #print(poslist)
          #print(pos)
          #stop("Error!")
        #}

        #add directed edge, two directions
        newstate <- state
        newstate$mg[i,j] <- 1
        if (sum(newstate$mg[,j]) <= max.in.degree && ggm::isAcyclic(newstate$mg)) {
          k <- part[j]
          newstate$scores[k] <- partScore(newstate$mg, part, k, n, covMat, edge.penalty, bic = bic)
          newstate$score <- sum(newstate$scores)
          i.a <- i.a + 1
          cand.add[[i.a]] <- newstate
        }


        newstate <- state
        newstate$mg[j,i] <- 1
        if (sum(newstate$mg[,i]) < max.in.degree && ggm::isAcyclic(newstate$mg)) {
          k <- part[i]
          newstate$scores[k] <- partScore(newstate$mg, part, k, n, covMat, edge.penalty, bic = bic)
          newstate$score <- sum(newstate$scores)
          i.a <- i.a + 1
          cand.add[[i.a]] <- newstate
        }

      }
    }

    #backward
    cand.del <- list()
    i.d <- 0
    if (direction != 1){
      poslist <- which(state$mg + t(state$mg) != 0 & lower.tri(state$mg))
      if (length(poslist)>1) poslist <- sample(poslist, min(length(poslist), max.pos))
      for (pos in poslist){
        #(i,j) entry
        j <- ceiling(pos / p)
        i <- (pos - 1) %% p + 1

        #delete edge
        newstate <- state
        newstate$mg[i,j] <- newstate$mg[j,i] <- 0
        k <- part[j]
        newstate$scores[k] <- partScore(newstate$mg, part, k, n, covMat, edge.penalty, bic = bic)
        newstate$score <- sum(newstate$scores)
        i.d <- i.d + 1
        cand.del[[i.d]] <- newstate
      }
    }

    #edge type change
    cand.cha <- list()
    i.c <- 0
    if (direction <= 3){
      # traverse directed edges
      poslist <- which(state$mg == 1)
      if (length(poslist)>0) poslist <- sample(poslist, min(length(poslist), max.pos))
      for (pos in poslist){
        #(i,j) entry
        j <- ceiling(pos / p)
        i <- (pos - 1) %% p + 1
        # check covered edge and prune equivalent DAGs
        if ((sum(part == part[i]) == 1) & (sum(part == part[j]) == 1) &
          (all(state$mg[, i] == state$mg[, j])))
          next()

        #reverse direction
        newstate <- state
        newstate$mg[i,j] <- 0
        newstate$mg[j,i] <- 1
        if (sum(newstate$mg[,i]) < max.in.degree && ggm::isAcyclic(newstate$mg)) {
          if (part[i] == part[j]){
            k <- part[i]
            newstate$scores[k] <- partScore(newstate$mg, part, k, n, covMat, edge.penalty, bic = bic)
          }
          else{
            k1 <- part[i]
            k2 <- part[j]
            newstate$scores[k1] <- partScore(newstate$mg, part, k1, n, covMat, edge.penalty, bic = bic)
            newstate$scores[k2] <- partScore(newstate$mg, part, k2, n, covMat, edge.penalty, bic = bic)
          }
          newstate$score <- sum(newstate$scores)
          i.c <- i.c + 1
          cand.cha[[i.c]] <- newstate
        }
      }
    }

    # Compare candidates and select best one
    if (i.a > 0){
      best.a <- which.max(sapply(1:i.a, function(j) cand.add[[j]]$score))
      best.a.score <- cand.add[[best.a]]$score
    } else{
      best.a <- NULL
      best.a.score <- -Inf
    }

    # if some block has no nodes, the score is NaN and best.d is integer(0)
    if (i.d > 0){
      best.d <- which.max(sapply(1:i.d, function(j) cand.del[[j]]$score))
      #if (isTRUE(best.d==0)) stop("!")
      #print(c(i.d, best.d))
      #print(cand.del[[2]]$scores)
      #print(part)
      #print(best.d)
      best.d.score <- cand.del[[best.d]]$score
    } else{
      best.d <- NULL
      best.d.score <- -Inf
    }

    if (i.c > 0){
      best.c <- which.max(sapply(1:i.c, function(j) cand.cha[[j]]$score))
      best.c.score <- cand.cha[[best.c]]$score
    } else{
      best.c <- NULL
      best.c.score <- -Inf
    }

    states[[iter]] <- state

    # Break out of loop, if no improvement of score is possible
    if (max(best.a.score, best.d.score, best.c.score) <= state$score + eps.conv) break

    #Otherwise change state to best candidate
    ind <- which.max(c(best.a.score, best.d.score, best.c.score))
    if (ind == 1){
      state <- cand.add[[best.a]]
    } else {
      if (ind == 2){
        state <- cand.del[[best.d]]
      } else{
        state <- cand.cha[[best.c]]
      }
    }

    # Take time for step
    state$t <- proc.time()[3] - t
    state$ct <- state$ct + state$t
    t <- proc.time()[3]

    iter <- iter + 1

    # Info output
    action <- c("ADD", "DEL", "CHA")[ind]
    lengths <- paste0("(", length(cand.add), ", ", length(cand.del), ", ", length(cand.cha), ")", sep="")
    if (verbose) print( paste(i, action, "; (adds, dels, changes) =", lengths, ";", "steptime =", round(state$t, 1),
                             "; score =", round(state$score, 4) ) )

  }

  if (iter >= maxSteps) print("MAX STEPS ACHIEVED")

  return (list(mg=state$mg, score=state$score, it=iter, states=states, cand.add=cand.add,
               cand.del=cand.del, cand.cha=cand.cha))
}


greedySearch <- function(
  cov.mat,
  n,
  mg.start = NULL,
  part,
  n.restarts = 1,
  mc.cores = 1,
  max.steps = 10,
  max.in.degree = 5,
  neighbourhood.size = Inf,  # max.pos
  eps.conv = 1e-12,
  dags.only = FALSE,
  direction = "both",
  verbose = FALSE,
  bic="std_bic"
)
{
  p <- ncol(cov.mat)

  if (is.null(colnames(cov.mat))) colnames(cov.mat) <- 1:ncol(cov.mat)
  if (is.null(rownames(cov.mat))) rownames(cov.mat) <- 1:nrow(cov.mat)

  # list of start mixed graphs
  mg.list <- if(is.null(mg.start)){
    GenerateGraph(p, N=n.restarts, p1=1, max.in.degree=2, names=rownames(cov.mat))
  } else{
    if (class(mg.start) == "list"){
      if (length(mg.start) == n.restarts) mg.start else sample(mg.start, n.restarts, replace = TRUE)
    } else{
      stop("Invalid starting graphs")
    }
  }

  directionMap <- c(
    "forward" = 1,
    "backward" = 2,
    "both" = 3
  )

  #one replication of greedy search, at start graph i
  oneRep <- function(i)
    fastGreedySearch(
      mg.list[[i]],
      part = part,
      n = n,
      maxSteps = max.steps,
      max.in.degree = max.in.degree,
      direction = directionMap[direction],
      covMat = cov.mat,
      verbose = verbose,
      max.pos = neighbourhood.size,
      eps.conv = eps.conv,
      bic = bic
    )

    res <- if(mc.cores > 1){
      doParallel::registerDoParallel(cores = mc.cores)
      foreach::foreach(i = 1:(n.restarts)) %dopar% { oneRep(i) }
    } else{
      lapply(1:(n.restarts), oneRep)
    }

    #return (res)
    i.best <- which.max(sapply(1:(n.restarts), function(j) res[[j]]$score))
    list(
      # must use "=", otherwise the attribute names are not saved
      final.mg = res[[i.best]]$mg,
      final.score = res[[i.best]]$score,
      #all.mgs <- lapply(1:n.restarts, function(obj){lapply(1:length(obj$states), function(state){state$mg})} ),
      all.mgs = lapply(res, function(obj) lapply(obj$states, function(state) state$mg) ),
      all.scores = lapply(res, function(obj) sapply(obj$states, function(state) state$score) ),
      all.times = lapply(res, function(obj) sapply(obj$states, function(state) state$ct))
    )
}
