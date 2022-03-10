# set vc (vertex color) and vfc (vertex boudary color) as vector of string for different partition blocks
# if the input is a cpdag, it will plot a cpdag
plotG <- function(
  dag,
  part=NULL,
  #part=1:ncol(dag),
  tcltk=FALSE,
  vfc="black",
  ...
){
  mg <- dag
  poslist <- which(lower.tri(dag) & dag == 1 & dag == t(dag))
  for (pos in poslist){
    #(i,j) entry
    j <- ceiling(pos / ncol(mg))
    i <- (pos - 1) %% ncol(mg) + 1
    mg[i, j] <- mg[j, i] <- 10
  }
  ggm::plotGraph(mg, layout=igraph::layout.circle, tcltk=tcltk, vc=part, vfc=vfc, directed=F, ...)
}
