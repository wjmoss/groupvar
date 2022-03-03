source("generate.R")
source("greedysearch.R")
source("plotting.R")


#GenerateGraph(p=5,N=2,iter=50)

#r=0
#r=r+1
#plotDAG(res[[r]]$graph$mg, part=part)
#plotDAG(res[[r]]$path$final.mg, part=part)


##First read in the arguments listed at the command line
args <- commandArgs(trailingOnly = TRUE)

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default value
  a = 1
  b = c(1,1,1)
  c='ss'
}else{
  for(i in 1:3){
    print(args[[i]])
    eval(parse(text=args[[i]]))
  }
}

#print(a*2)
#print(b*3)
#print(c)
#d = args[[4]]
#print(eval(parse(text=d)))
#.libPaths()
#version

#print(args[1])
#print(args[2])