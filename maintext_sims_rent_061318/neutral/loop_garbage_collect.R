
#Loop through saved environments and delete extraneous objects.

N <- 10000
herit <- 1
out_dir <- "out_061318/"

#nsamps <- 20
sel.intenses <- .005
n.locis <- c(100)
n_chromss <- c(200)
ts <- c(0.0)
t.offs <- c(0.0)
#phen_nums <- 1:1000 #1 number for each rep we want to do at each combination of parameters
phen_nums <- 1:100 #1 number for each rep we want to do at each combination of parameters


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
pars.gc <- expand.grid(sel.intenses, n.locis, n_chromss, ts, t.offs, phen_nums)

#objects not to delete
essent.obs <- c("anc_trees_ms","anc_trees_rent","ct_ms_all","ct_ms_anc","ct_ms_der",
"ct_rent_all","ct_rent_anc", "ct_rent_der", "curr.freqs",             
"der_trees_ms", "der_trees_rent",      
"eff_sizes", "eff_sizes_unscaled", "herit", "len_hap", "lins.list",   
"lins.list.rent", "ms_haplotypes_list", "ms_trees_list", "N", "n_chroms", "n_ders",            
"n.loci", "out_dir","phen_num","phen.traj", "pt.time", "r",                   
"rent_trees_list", "sd.trait", "sel_site", "sel.intens", "sel.shift", "t", "t.off", "target.shift", "time", "times.c", "times.c.rent", "trait.sd", "trajs", "u", "lasttree")
   

for(iter in 1:dim(pars.gc)[1]){
	sel.intens <- pars.gc[iter,1]
	n.loci <- pars.gc[iter,2]
	n_chroms <- pars.gc[iter,3]
	t <- pars.gc[iter,4]
	t.off <- pars.gc[iter,5]
	phen_num <- pars.gc[iter,6]
	fn <- paste(out_dir, "loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), ".RData", sep = "")
	e <- new.env()	
	load(fn, envir = e)
	if(length(setdiff(ls(e), essent.obs)) >  0){	
		rm(list = setdiff(ls(e), essent.obs), envir = e)
		save(list = ls(e), file = fn, envir = e)
	}
	print(paste("gc for trait", as.character(iter), "complete."))
}


