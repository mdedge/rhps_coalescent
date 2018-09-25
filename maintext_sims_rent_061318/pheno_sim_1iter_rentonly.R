

#Doc Edge, 6/13/18
#Goal: Open up a file with saved trees and just run rent.


fn_str <- paste("loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")


load(paste(out_dir, fn_str, ".RData", sep = ""))

rent_in_fn <- rent_in_fn.ro

rent_trees_list <- list()


#for each locus, write trajectory, run ms to simulate a sample and tree,
#and run rent+ to infer tree. Save both the ms and the rent+ trees, as well
#as the number of derived alleles for each locus.
for(i in 1:n.loci){
	write(ms_haplotypes_list[[i]], rent_in_fn)
	
	#run rent+
	rent.string <- paste("java -jar ", rentplus_fn, " -t -l ", as.character(len_hap), " ", rent_in_fn, sep = "")
	system(rent.string)

	#Read in trees at selected site in rent and rent+
	rent_trees_fn <- paste(rent_in_fn, ".trees", sep = "")
	rent_trees <- readLines(rent_trees_fn)
	#which tree is for the selected site?
	rent_sel_ind <- grep(paste(as.character(sel_site), "\t", sep = ""), rent_trees)

	#get tree at selected site in newick format, both true and rent-estimated.
	rent_sel_tree_newick <- sub(paste(as.character(sel_site), "\t", sep = ""), "", rent_trees[rent_sel_ind])
	rent_trees_list[[i]] <- read.tree(text = paste(rent_sel_tree_newick, ";", sep = ""))
	print(paste("simulation locus", as.character(i), "of trait", as.character(phen_num), "complete"))
}


#Run through list and check for sites where 2+ variants had same coordinates as
#selected site. At those places, check which tree(s) are monophyletic for derived tips.
#if more than one, select randomly from among them.
for(i in 1:n.loci){
	if(is.null(names(rent_trees_list[[i]]))){
		cand.trees <- rent_trees_list[[i]]
		is.mono <- numeric(0)		
		for(j in 1:length(cand.trees)){
			der_tips <- which(cand.trees[[j]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
			is.mono[j] <- is.monophyletic(cand.trees[[j]], der_tips)
		}
		if(mean(is.mono) > 0){
			ind.to.take <- 	sample(which(is.mono == 1), 1)
		}		
		if(mean(is.mono) == 0){ #if none are monophyletic, choose at random
			ind.to.take <- sample(1:length(cand.trees), 1)
		}
		rent_trees_list[[i]] <- rent_trees_list[[i]][[ind.to.take]]
	}
}



#split into ancestral and derived trees and retrieve coalescence times.

anc_trees_rent <- list()
der_trees_rent <- list()
ct_rent_all <- list()
ct_rent_anc <- list()
ct_rent_der <- list()

for(i in 1:n.loci){
	anc_tips_rent <- which(rent_trees_list[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
	der_tips_rent <- which(rent_trees_list[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
	anc_trees_rent[[i]] <- drop.tip(rent_trees_list[[i]], der_tips_rent)
	der_trees_rent[[i]] <- drop.tip(rent_trees_list[[i]], anc_tips_rent)
	ct_rent_all[[i]] <- coal.times(rent_trees_list[[i]])
	ct_rent_anc[[i]] <- coal.times(anc_trees_rent[[i]])
	ct_rent_der[[i]] <- coal.times(der_trees_rent[[i]])
	print(paste("round", as.character(i), "of tree processing complete"))
}



times.c.rent <- list()
for(i in 1:length(rent_trees_list)){
	times.c.rent[[i]] <- trees_to_times(rent_trees_list[[i]], anc_trees_rent[[i]], der_trees_rent[[i]], time, sure.alt.is.derived = FALSE, units_in = 2)
}

lins.list.rent <- list()
for(i in 1:length(rent_trees_list)){
	lins.list.rent[[i]] <- times_to_lins(times.c.rent[[i]], time)
}




save.image(paste(out_dir, fn_str, ".RData", sep = ""))




