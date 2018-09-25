
#complete phenotype simulations by looping through parameter values and calling 
#pheno_sim_1iter.R for each set of parameters.

N <- 10000
herit <- 1
#out_dir <- "out_061318/"
out_dir <- "out_09test/"

#nsamps <- 20
#sel.intenses <- c(.001,.005, .01)
sel.intenses <- .01
n.locis <- c(100)
n_chromss <- c(1000)
#ts <- c(.001, 0.0025, .005)
ts <- .095
#t.offs <- c(0.0)
t.offs <- .09
#phen_nums <- 1:40 #1 number for each rep we want to do at each combination of parameters
phen_nums <- 1:20

traj.fn <- "temp_out9/temp.txt"
msout.fn <- "temp_out9/ms_out.txt"
rent_in_fn <- "temp_out9/rent_in.txt"


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
pars <- expand.grid(sel.intenses, n.locis, n_chromss, ts, t.offs, phen_nums)

qxtest_mat <- matrix(nrow = nrow(pars), ncol = 4)

for(k in 1:dim(pars)[1]){
	sel.intens <- pars[k,1]
	n.loci <- pars[k,2]
	n_chroms <- pars[k,3]
	t <- pars[k,4]
	t.off <- pars[k,5]
	phen_num <- pars[k,6]
	#source("../pheno_sim_1iter.R")
	source("../pheno_sim_1iter_norent.R")

	trajs_neut <- matrix(nrow = length(time), ncol = length(ms_trees_list))
	for(i in 1:length(ms_trees_list)){	
		trajs_neut[,i] <- est_af_traj_neut(lins.list[[i]])
	}
	traj.phen.neut <- 2 * trajs_neut %*%  eff_sizes 
	
	qxtest_mat[k,1:3] <- Qx_test(trajs_neut[time %in% ((0:10)/100),], eff_sizes, perms = 0)
	qxtest_mat[k,4] <- Qx_test(trajs_neut[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
	#print(qxtest_mat)
	print(paste("trial", as.character(k), "complete."))
}


save.image(paste("test_grads", ".RData", sep = ""))





