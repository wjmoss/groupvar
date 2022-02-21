source("generate.R")
source("greedysearch.R")
source("plotting.R")


#GenerateGraph(p=5,N=2,iter=50)

r=0
r=r+1
plotDAG(res[[r]]$graph$mg, part=part)
plotDAG(res[[r]]$path$final.mg, part=part)
