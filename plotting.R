# set vc (vertex color) and vfc (vertex boudary color) as vector of string for different partition blocks
plotDAG <- function(
  dag,
  part=1:ncol(dag),
  tcltk=FALSE,
  vfc="black",
  ...
) ggm::plotGraph(dag, layout=igraph::layout.circle, tcltk=tcltk, vc=part, vfc=vfc, ...)

