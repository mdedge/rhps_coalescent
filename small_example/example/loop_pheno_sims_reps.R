
#complete phenotype simulations by looping through parameter values and calling 
#pheno_sim_1iter.R for each set of parameters.
N <- 10000
herit <- 1
out_dir <- "out/"

sel.intenses <- .005
n.locis <- 20
n_chromss <- 30
ts <- 0.04
t.offs <- 0.02
phen_nums <- 1:5 #1 number for each rep we want to do at each combination of parameters

traj.fn <- "temp/temp.txt"
msout.fn <- "temp/ms_out.txt"
rent_in_fn <- "temp/rent_in.txt"


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


for(k in 1:dim(pars)[1]){
	sel.intens <- pars[k,1]
	n.loci <- pars[k,2]
	n_chroms <- pars[k,3]
	t <- pars[k,4]
	t.off <- pars[k,5]
	phen_num <- pars[k,6]
	source("../pheno_sim_1iter.R")
	#source("../pheno_sim_1iter_norent.R")
	print(paste("trial", as.character(phen_num), "complete."))
}


save.image(paste("sim_trees", ".RData", sep = ""))




err.array <- array(dim = c(length(time),5,dim(pars)[1]))   
err.std.array <- array(dim = c(length(time),3,dim(pars)[1]))  
qxtest_mat <- matrix(-1, nrow = dim(pars)[1], ncol = 12)

err.array.rent <- array(dim = c(length(time),3,dim(pars)[1]))   
err.std.array.rent <- array(dim = c(length(time),3,dim(pars)[1]))  
qxtest_mat_rent <- matrix(-1, nrow = dim(pars)[1], ncol = 12)
mat.true.phentrajs <- matrix(nrow = length(time), ncol = dim(pars)[1])

for(iter in 1:dim(pars)[1]){
	sel.intens <- pars[iter,1]
	n.loci <- pars[iter,2]
	n_chroms <- pars[iter,3]
	t <- pars[iter,4]
	t.off <- pars[iter,5]
	phen_num <- pars[iter,6]
	fn <- paste(out_dir, "loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), ".RData", sep = "")
	load(fn)
	source("../analyze_sim_true.R")
	source("../analyze_sim_rent.R")
	print(paste("trial", as.character(iter), "complete."))
}


save.image(paste("analyzed_trees", ".RData", sep = ""))

source("make_figures.R")

