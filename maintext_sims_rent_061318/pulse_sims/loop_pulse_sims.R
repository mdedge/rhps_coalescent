
#complete phenotype simulations by looping through parameter values and calling 
#pheno_sim_1iter.R for each set of parameters.

N <- 10000
herit <- 1
out_dir <- "out_061318/"

#nsamps <- 20
sel.intenses <- .01
n.locis <- c(100)
n_chromss <- c(20, 50, 100, 200, 500, 1000)
t.offs <- .08
ts <- t.offs + .005
phen_nums <- c(1:100) #1 number for each rep we want to do at each combination of parameters

traj.fn <- "temp_out10/temp.txt"
msout.fn <- "temp_out10/ms_out.txt"
rent_in_fn <- "temp_out10/rent_in.txt"


helper_fn <- "../../helper_functions_coal_sel.R"
source(helper_fn) #read in helper functions and load ape package

ms_dir <- "../../msseldir/"
rentplus_fn <- "../../RentPlus.jar"
len_hap <- 200000 #the length (in base pairs) of the haplotype -- short because 
#we are not worrying about recombination--just need the sel site.
sel_site <- 100000 #the position of the selected site in the haplotype
u <- 2e-8 #the neutral mutation rate per base pair/generation
r <- 2.5e-8 #the recombination rate per base pair/generation
options(scipen = 999) #disable scientific notation so that parameters passed to ms are given as numbers
sd.trait <- 1

time <- c(seq(0, 4, by = 0.001))
pars <- expand.grid(n_chromss, sel.intenses, n.locis, ts, t.offs, phen_nums)
qxtest_mat <- matrix(nrow = nrow(pars), ncol = 4)


for(k in 1:dim(pars)[1]){
	n_chroms <- pars[k,1]
	sel.intens <- pars[k,2]
	n.loci <- pars[k,3]	
	t <- pars[k,4]
	t.off <- pars[k,5]
	phen_num <- pars[k,6]
	if((k - 1) %% length(n_chromss) == 0){
		print(k)
		print(paste("simulating polygenic score trajectory", as.character(t.off)))
		source("pheno_sim_aftrajs.R")
		fn_str <- paste("trajs_loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")
	save(list = c("trajs", "time", "eff_sizes"),file = paste(out_dir, fn_str, ".RData", sep = ""))
	}	
	source("pheno_sim_trees.R")

	print(paste("trial", as.character(k), "complete."))
}

fn_str <- paste("qxmat_ton", as.character(t), "_toff", as.character(t.off), "phen_nums_", as.character(min(phen_nums)), "_to_", max(phen_nums), sep = "")
save(list = c("pars", "qxtest_mat"), file =  paste(fn_str, ".RData", sep = ""))





