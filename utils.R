## some util functions

# modified from pcalg.R, line 590
computeCPDAG <- function(mg, phase=1){
  if (phase == 1){
    # without backgound knowledge, mg is a dag and need to be changed to a cpdag
    mg[mg != 0] <- 1
    skel <- mg + t(mg)
    skel[skel == 2] <- 1
    cpdag <- skel
  }
  else {
    # with background knowledge, mg is a pdag
    cpdag <- mg
  }


  ## search the v-structures in the DAG
  ind <- which((mg == 1 & t(mg) == 0), arr.ind = TRUE)
  tripleMatrix <- matrix(,0,3)
  ## Go through all edges
  for (i in seq_len(nrow(ind))) { ## MM(FIXME): growth of tripleMatrix
    x <- ind[i,1]
    y <- ind[i,2]
    indY <- setdiff(which((mg[,y] == 1 & mg[y,] == 0), arr.ind = TRUE),x) ## x-> y <- z
    if(length(newZ <- indY[mg[x,indY] == 0])) ## deparse.l.=0: no colnames
      tripleMatrix <- rbind(tripleMatrix, cbind(x, y, newZ, deparse.level=0), deparse.level=0)
  }
  if ((m <- nrow(tripleMatrix)) > 0) {
    deleteDupl <- logical(m)# all FALSE
    for (i in seq_len(m))
        if (tripleMatrix[i,1] > tripleMatrix[i,3])
          deleteDupl[i] <- TRUE
    if(any(deleteDupl))
      tripleMatrix <- tripleMatrix[!deleteDupl,, drop=FALSE]

    ## orient the v-structures in the CPDAG
    for (i in seq_len(nrow(tripleMatrix))) {
      x <- tripleMatrix[i,1]
      y <- tripleMatrix[i,2]
      z <- tripleMatrix[i,3]
      cpdag[x,y] <- cpdag[z,y] <- 1
      cpdag[y,x] <- cpdag[y,z] <- 0
    }
  }

  ## orient the edges with the 3 orientation rules
  repeat{
    old_cpdag <- cpdag
    ## rule 1
    ind <- which((cpdag == 1 & t(cpdag) == 0), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))){
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- ((cpdag[b, ] == 1 & cpdag[, b] == 1) &
              (cpdag[a, ] == 0 & cpdag[, a] == 0))
      if (any(isC)){
        indC <- which(isC)
        cpdag[b, indC] <- 1
        cpdag[indC, b] <- 0
      }
    }
    ## rule 2
    ind <- which((cpdag == 1 & t(cpdag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))){
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- ((cpdag[a, ] == 1 & cpdag[, a] == 0) &
              (cpdag[b, ] == 0 & cpdag[, b] == 1))
      if (any(isC)){
        cpdag[a, b] <- 1
        cpdag[b, a] <- 0
      }
    }
    ## rule 3
    ind <- which((cpdag == 1 & t(cpdag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))){
      a <- ind[i, 1]
      b <- ind[i, 2]
      indC <- which((cpdag[a, ] == 1 & cpdag[, a] == 1) &
              (cpdag[b, ] == 0 & cpdag[, b] == 1))
      if (length(indC) >= 2) {
        cmb.C <- combn(indC, 2) # all 2-pairs in indC
        cC1 <- cmb.C[1, ]
        cC2 <- cmb.C[2, ]
        for (j in seq_along(cC1)){
          c1 <- cC1[j]
          c2 <- cC2[j]
          if (c1 != c2 && cpdag[c1, c2] == 0 && cpdag[c2,c1] == 0) {
            cpdag[a, b] <- 1
            cpdag[b, a] <- 0
            break
          }
        }
      }
    }
    ## if phase == 2, rule 4
    if (phase == 2){
      ind <- which((cpdag == 1 & t(cpdag) == 1), arr.ind = TRUE)
      for (i in seq_len(nrow(ind))){
        a <- ind[i, 1]
        b <- ind[i, 2]
        indC <- which((cpdag[a, ] == 1 & cpdag[, a] == 1) &
          (cpdag[b, ] == 0 & cpdag[, b] == 1))
        for (c in indC){
          indD <- ((cpdag[c, ] == 0 & cpdag[, c] == 1) &
            (cpdag[a, ] == 1 & cpdag[, a] == 1))
          if (any(indD)){
            cpdag[a, b] <- 1
            cpdag[b, a] <- 0
            break
          }
        }
      }
    }
    if (all(cpdag == old_cpdag))
      break
  }
  return (cpdag)
}


# cpdag with background knowledge
computeCPDAG_bk <- function(mg, part=NULL){
  cpdag <- computeCPDAG(mg, phase=1)
  if (is.null(part))
    return (cpdag)
  else{
    for (i in unique(part)){
      ind <- which(part == i)
      if (length(ind) > 1){
        for (j in ind){
          cpdag[j, ] <- mg[j, ]
          cpdag[, j] <- mg[, j]
        }
      }
    }
    #return (cpdag)
    return (computeCPDAG(cpdag, phase=2))
  }
}


# shd with dist += 2 for an edge reversal, modified from pcalg
computeSHD <- function(mg1, mg2, reversal=1){
  if (nrow(mg1) != ncol(mg1) | nrow(mg2)!=ncol(mg2))
    stop("DAG or CPDAG must be a square matrix!")
  if (nrow(mg1) != nrow(mg2))
    stop("The two graph must have the same number of nodes!")
  shd <- 0
  s1 <- mg1 + t(mg1)
  s2 <- mg2 + t(mg2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  mg1[ind] <- 0
  shd <- shd + length(ind) / 2
  ind <- which(ds < 0)
  mg1[ind] <- mg2[ind]
  shd <- shd + length(ind) / 2
  d <- abs(mg1 - mg2)
  if (reversal == 1){
    shd + sum((d + t(d)) > 0)/2 #this line is for classical shd, same as that in pcalg
  } else {
    shd + sum((d + t(d)) > 0) # dist+=2 for each pair of reversed edges
  }
}


## cpdag to a representative dag:
# helper
validation <- function(mg, i, j){
  if (mg[i, j] == 0){
    # if not i->j, exchange i and j
    tmp <- j
    j <- i
    i <- tmp
  }
  pa <- which(mg[, j] == 1)
  pa <- pa[pa != i]
  #pa <- which(mg[-i, j] == 1) # indices may change!
  #pa <- setdiff(pa, i)

  # new unshielded collider triples -- return false
  if (length(pa) > 0 & any( (mg+t(mg))[pa, i] == 0) )
    return (F)

  tmp <- mg
  tmp[tmp == 10] <- 0
  if (ggm::isAcyclic(tmp)){
    return (T)
  }
  return (F)
}


# generate a DAG of the given CPDAG
cpdag_repr <- function(cpdag){
  mg <- cpdag
  mg[mg > 0 & mg == t(mg)] <- 10
  poslist <- which(mg == 10 & t(mg) == 10 & lower.tri(mg))
  k <- 1

  while (k > 0){
    if (k == length(poslist) + 1)
      return (mg)

    j <- ceiling(poslist[k] / ncol(mg))
    i <- (poslist[k] - 1) %% ncol(mg) + 1
    if (mg[i, j] == 10 & mg[j, i] == 10){
      mg[i, j] <- 1
      mg[j, i] <- 0
    } else {
      if (mg[i, j] == 1){
        mg[i, j] <- 0
        mg[j, i] <- 1
      } else {
        mg[i, j] <- mg[j, i] <- 10
        k <- k - 1
        next()
      }

    }
    if (validation(mg, i, j)){
      k <- k + 1
    }
  }
}