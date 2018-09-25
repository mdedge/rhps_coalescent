

#Doc Edge, 6/13/18
#given a locus, simulate and run rent.
#This is called by the script that attempts to rerun rent+ loci that had problems.



#for each locus, write trajectory, run ms to simulate a sample and tree,
#and run rent+ to infer tree. Save both the ms and the rent+ trees, as well
#as the number of derived alleles for each locus.

	write.traj(traj.fn.rf, trajs[[i]])
	curr.freq.der <- curr.freqs[i]
	n_der <- rbinom(1, n_chroms, curr.freq.der) #number of chroms with derived allele.
	if(n_der == 0){n_der <- n_der + 1}
	if(n_der == n_chroms){n_der <- n_der - 1}
	n_ders[i] <- n_der
	#run mssel
	ms.string <- paste(ms_dir, "mssel ", as.character(n_chroms), " 1 ", as.character(n_chroms - n_der), " ", as.character(n_der), " ", traj.fn.rf, " ", sel_site, " -r ", as.character(4*N*(1 - dbinom(0,len_hap, r))), " ", as.character(len_hap), " -t ", as.character(4*N* (1 - dbinom(0, len_hap, u))), " -T -L > ", msout.fn.rf, sep = "")
	system(ms.string)

	#cut down mssel output to just the part rent+ needs
	process.string <- paste("sed '/positions: /,$!d' ", msout.fn.rf, " | sed 's/.*positions: //' > ", rent_in_fn.rf, sep = "")
	system(process.string)

	if(file.exists(paste(rent_in_fn.rf, ".trees", sep = ""))){
		file.remove(paste(rent_in_fn.rf, ".trees", sep = ""))
	}

	#run rent+
	rent.string <- paste("java -jar ", rentplus_fn, " -t -l ", as.character(len_hap), " ", rent_in_fn.rf, sep = "")
	system(rent.string)

	ms_haplotypes_list[[i]] <- readLines(rent_in_fn.rf) #save the haplotypes produced by ms (and fed to rent+) in a list entry.
	#the result will be a character vector. The first entry in the vector is relative (0-1) positions, with the selected
	#site at 0.05. The remaining entries are string haplotypes of 0s and 1s.
	

	#Read in trees at selected site in rent+
	rent_trees_fn <- paste(rent_in_fn.rf, ".trees", sep = "")

	if(file.exists(rent_trees_fn)){
		rent_trees <- readLines(rent_trees_fn)
	#which tree is for the selected site?
		rent_sel_ind <- grep(paste(as.character(sel_site), "\t", sep = ""), rent_trees)

	#get tree at selected site in newick format, both true and rent-estimated.
		rent_sel_tree_newick <- sub(paste(as.character(sel_site), "\t", sep = ""), "", rent_trees[rent_sel_ind])
		rent_trees_list[[i]] <- read.tree(text = paste(rent_sel_tree_newick, ";", sep = ""))

	}

		
	#pull trees from ms
	msoutlines <- readLines(msout.fn.rf)
	dists <- numeric(0)
	for(m in 1:(length(msoutlines) - 4)){
		suppressWarnings(dists[m] <- as.numeric(substring(strsplit(msoutlines[m + 4], "]")[[1]][1],2)))
	}
	dists <- dists[!is.na(dists) & dists != Inf]
	posits.rs <- cumsum(dists)
	ms_to_extract <- which.max(posits.rs[posits.rs < sel_site]) + 1
	if(length(ms_to_extract)  == 0){ms_to_extract <- 1}
	ms_trees_list[[i]] <- read.tree(text = msoutlines[4 + ms_to_extract])		

	print(paste("simulation locus", as.character(i), "of trait", as.character(flop), "complete"))



#Run through list and check for sites where 2+ variants had same coordinates as
#selected site. At those places, check which tree(s) are monophyletic for derived tips.
#if more than one, select randomly from among them.

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
		ms_trees_list[[i]] <- cand.trees[[ind.to.take]]
		rent_trees_list[[i]] <- rent_trees_list[[i]][[ind.to.take]]
	}




#split into ancestral and derived trees and retrieve coalescence times.


	anc_tips_rent <- which(rent_trees_list[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
	anc_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
	der_tips_rent <- which(rent_trees_list[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
	der_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
	anc_trees_rent[[i]] <- drop.tip(rent_trees_list[[i]], der_tips_rent)
	der_trees_rent[[i]] <- drop.tip(rent_trees_list[[i]], anc_tips_rent)
	anc_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], der_tips_ms)
	der_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], anc_tips_ms)
	ct_rent_all[[i]] <- coal.times(rent_trees_list[[i]])
	ct_rent_anc[[i]] <- coal.times(anc_trees_rent[[i]])
	ct_rent_der[[i]] <- coal.times(der_trees_rent[[i]])
	ct_ms_all[[i]] <- coal.times(ms_trees_list[[i]])
	ct_ms_anc[[i]] <- coal.times(anc_trees_ms[[i]])
	ct_ms_der[[i]] <- coal.times(der_trees_ms[[i]])
	print(paste("round", as.character(i), "of tree processing complete"))







times.c[[i]] <- trees_to_times(ms_trees_list[[i]], anc_trees_ms[[i]], der_trees_ms[[i]], time, sure.alt.is.derived = FALSE, units_in = 4)



lins.list[[i]] <- times_to_lins(times.c[[i]], time)




times.c.rent[[i]] <- trees_to_times(rent_trees_list[[i]], anc_trees_rent[[i]], der_trees_rent[[i]], time, sure.alt.is.derived = FALSE, units_in = 2)


lins.list.rent[[i]] <- times_to_lins(times.c.rent[[i]], time)







