

#Run through phenotype simulations and fix rent+

#complete phenotype simulations by looping through parameter values and calling 
#pheno_sim_1iter.R for each set of parameters.

N <- 10000
herit <- 1
out_dir <- "out_061318/"

#nsamps <- 20
sel.intenses <- .005
n.locis <- c(100)
n_chromss <- c(200)
ts <- c(0.0)
t.offs <- c(0.0)
phen_nums <- 1:100 #1 number for each rep we want to do at each combination of parameters

traj.fn.rf <- "temp_out9/temp.txt"
msout.fn.rf <- "temp_out9/ms_out.txt"
rent_in_fn.rf <- "temp_out9/rent_in.txt"


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
pars.fr <- expand.grid(sel.intenses, n.locis, n_chromss, ts, t.offs, phen_nums)



#Feed in a list of trees (and optionally the last locus tree of the previous trial)
#and get back a vector of the ones that failed.
find_fails <- function(trees, last_tree = NULL){
	fails <- numeric(0)
	if(is.null(last_tree)){
		fails[1] <- FALSE	
	}
	if(!is.null(last_tree)){
		fails[1] <- all.equal(last_tree, trees[[1]])	
	}	
	for(i in 2:length(trees)){
		fails[i] <- all.equal(trees[[i-1]], trees[[i]])
	}
	fails
}


lasttree <- NULL

for(flop in 1:dim(pars.fr)[1]){
	sel.intens <- pars.fr[flop,1]
	n.loci <- pars.fr[flop,2]
	n_chroms <- pars.fr[flop,3]
	t <- pars.fr[flop,4]
	t.off <- pars.fr[flop,5]
	phen_num <- pars.fr[flop,6]
	#source("../pheno_sim_1iter.R")
	
	fn_str <- paste("loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")
	load(paste(out_dir, fn_str, ".RData", sep = ""))
	
	failed_loci <- which(find_fails(rent_trees_list, lasttree) == 1)
	print(failed_loci)
	
	for(locus in failed_loci){
		i <- locus
		source("../simwrent_onelocus.R")
	}
	if(length(failed_loci) > 0){
		save.image(paste(out_dir, fn_str, ".RData", sep = ""))
	}
	#lasttree <- rent_trees_list[[n.loci]]	
	print(paste("trial", as.character(flop), "complete."))
}


