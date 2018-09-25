
#Simulate a bunch of trees under neutrality and compute SDS, saving mean and sd.

N <- 10000
herit <- 1
out_dir <- "out_073118/"

#nsamps <- 20
sel.intense <- .01
#n.locis <- c(100)
n.loci <- 500
n_chromss <- c(20, 50, 100, 200, 500, 1000)
t <- 0
t.off <- 0
#t.offs <- .08
#ts <- t.offs + .005
dafs <- (1:199)/200

traj.fn <- "temp_out/temp.txt"
msout.fn <- "temp_out/ms_out.txt"

helper_fn <- "../../helper_functions_coal_sel.R"
source(helper_fn) #read in helper functions and load ape package

ms_dir <- "../../msseldir/"
rentplus_fn <- "../../RentPlus.jar"
len_hap <- 2000 #the length (in base pairs) of the haplotype -- short because 
#we are not worrying about recombination--just need the sel site.
sel_site <- 1000 #the position of the selected site in the haplotype
u <- 2e-8 #the neutral mutation rate per base pair/generation
r <- 0 #the recombination rate per base pair/generation
options(scipen = 999) #disable scientific notation so that parameters passed to ms are given as numbers
sd.trait <- 1

time <- c(seq(0, 4, by = 0.001))
pars <- expand.grid(n_chromss, dafs)
tbl.list <- list()

mean.sds <- numeric(dim(pars)[1])
sd.sds <- numeric(dim(pars)[1])


for(k in 1:dim(pars)[1]){
	n_chroms <- pars[k,1]
	daf <- pars[k,2]
	if((k - 1) %% length(n_chromss) == 0){
		print(paste("simulating allele frequency trajectories", as.character(daf)))
		source("pheno_sim_aftrajs.R")
		fn_str <- paste("trajs_daf", as.character(daf), sep = "")
		save(list = c("trajs"),file = paste(out_dir, fn_str, ".RData", sep = ""))
	}	
	source("pheno_sim_trees.R")
	tbl.list[[k]] <- mat.tbl
	mean.sds[k] <- mean(mat.tbl[,5])
	sd.sds[k] <- sd(mat.tbl[,5])
	print(paste("daf", as.character(daf), "sample size", n_chroms, "complete."))
}

save(list = c("n.loci", "pars", "tbl.list", "mean.sds", "sd.sds"), file = "neutral_SDS_sims_summary.RData")

#fn_str <- paste("qxmat_ton", as.character(t), "_toff", as.character(t.off), "phen_nums_", as.character(min(phen_nums)), "_to_", max(phen_nums), sep = "")
#save(list = c("pars", "qxtest_mat"), file =  paste(fn_str, ".RData", sep = ""))





