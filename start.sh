#! /bin/bash

pvals=(5 10 20) #40)
nvals=(100 500 1000)
spvals=("s" "d")

for ((i1=0; i1<${#pvals[@]} ;i1++)) 
do
for ((i2=0; i2<${#nvals[@]} ;i2++))
do
for ((i3=0; i3<${#spvals[@]} ;i3++))
do

  p=${pvals[$i1]}
  n=${nvals[$i2]}
  sp=${spvals[$i3]}
  seed=1
  
  fnm="./logging/$p-$n-$sp-$seed-empty-stdbic.out"
  OMP_NUM_THREADS=1 nohup nice Rscript simdata.R $p $n $sp 2>&1 | tee $fnm &
  echo "OMP_NUM_THREADS=1 nohup nice Rscript simdata.R $p $n $sp 2>&1 | tee $fnm &"

done
done  
done
