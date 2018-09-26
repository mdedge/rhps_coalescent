This directory contains code to run a small working example of simulations like those in figures 3-5 of the paper. You can run the example if you have the repository downloaded (not necessary to include the maintext_sims or height_analyses directories), you have R installed (with the ape package also installed), and you have a version of java compatible with RENT+. To run the example, navigate to the small_example/example folder and open R. Then run

```
source("loop_pheno_sims_reps.R")
```

The script will, for several polygenic scores, simulate allele frequency histories and trees, estimate trees using RENT+ at each locus, and run the estimators and tests proposed in the paper on both the true and RENT+-estimated trees. Plots analogous to figures 3-5 will be produced with the trajectories for a single trial, as well as bias/mean squared error and confidence interval coverage. (The bias/MSE and confidence interval plots will look noisy if few trials are run, as is the default.) By default, the parameters  

```
N <- 10000
sel.intenses <- .005
n.locis <- 20
n_chromss <- 30
ts <- 0.04
t.offs <- 0.02
phen_nums <- 1:5 
```

specify a constant effective population size of 10000, a selection intensity of .005, 20 loci per polygenic score/trait, 30 chromosomes per simulated sample, selection that is "on" between .02 and .04 coalescent units in the past, and to run 5 simulated trials. The number of simulations, loci per trait, and chromosomes per sample are all smaller than in the paper so that this example will run relatively quickly. On my machine, it runs in a few (~10) minutes.
