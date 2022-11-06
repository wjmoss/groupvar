# groupvar

Greedy search algorithm for Gaussian linear structural equation models with groupwise equal variances.


# files
generate.R: graph and data generation
greedysearch.R: search algorithm
simdata.R: currently used main function for simulations
start.sh: start the simulations

make_plot.R, theme_for_plots.R: result visualisation

plotting.R: to plot CPDAG, could be put in utils?
utils.R: some util functions


# simulations
To run simulations with [p/3]+1 blocks (comparing with standard GES and PC algorithm):

```R
bash start.sh
```
