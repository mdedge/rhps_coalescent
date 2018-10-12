
#complete phenotype simulations by looping through parameter values and calling 
#pheno_sim_1iter.R for each set of parameters.

N <- 10000
herit <- 1
out_dir <- "out_092618/"

sel.intenses <- .005
n.locis <- c(1000)
n_chromss <- c(200)
ts <- c(0.02)
t.offs <- c(0.0)
phen_nums <- c(1:1000) #1 number for each rep we want to do at each combination of parameters


helper_fn <- "../../helper_functions_coal_sel.R"
source(helper_fn) #read in helper functions and load ape package
helper_fn_stable <- helper_fn

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
pars.mat <- expand.grid(sel.intenses, n.locis, n_chromss, ts, t.offs, phen_nums)

mat.true.phentrajs <- matrix(-1000, nrow =length(time), ncol = dim(pars.mat)[1])

err.array <- array(dim = c(length(time),5,dim(pars.mat)[1]))   
err.std.array <- array(dim = c(length(time),3,dim(pars.mat)[1]))  
qxtest_mat <- matrix(-1, nrow = dim(pars.mat)[1], ncol = 12)

#err.array.rent <- array(dim = c(length(time),3,dim(pars.mat)[1]))   
#err.std.array.rent <- array(dim = c(length(time),3,dim(pars.mat)[1]))  
#qxtest_mat_rent <- matrix(-1, nrow = dim(pars.mat)[1], ncol = 12)


for(iter in 1:dim(pars.mat)[1]){
	sel.intens <- pars.mat[iter,1]
	n.loci <- pars.mat[iter,2]
	n_chroms <- pars.mat[iter,3]
	t <- pars.mat[iter,4]
	t.off <- pars.mat[iter,5]
	phen_num <- pars.mat[iter,6]
	fn <- paste(out_dir, "loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), ".RData", sep = "")
	load(fn)
	source(helper_fn_stable)
	source("../analyze_sim_true.R")
	#source("../analyze_sim_rent.R")
	print(paste("trial", as.character(iter), "complete."))
}


save.image(paste("analyzed_true_trees_100518", ".RData", sep = ""))
#save.image(paste("analyzed_rent_trees_072718", ".RData", sep = ""))





