

#Doc Edge, 6/13/18
#Goal: Run through ms and save trees to estimate phenotype trajectory using coalescent times.




n_ders <- numeric(n.loci)
ms_trees_list <- list()
#ms_haplotypes_list <- list()


#for each locus, write trajectory, run ms to simulate a sample and tree,
#and run rent+ to infer tree. Save both the ms and the rent+ trees, as well
#as the number of derived alleles for each locus.
for(i in 1:n.loci){
	write.traj(traj.fn, trajs[[i]])
	curr.freq.der <- curr.freqs[i]
	n_der <- rbinom(1, n_chroms, curr.freq.der) #number of chroms with derived allele.
	if(n_der == 0){n_der <- n_der + 1}
	if(n_der == n_chroms){n_der <- n_der - 1}
	n_ders[i] <- n_der
	#run mssel
	ms.string <- paste(ms_dir, "mssel ", as.character(n_chroms), " 1 ", as.character(n_chroms - n_der), " ", as.character(n_der), " ", traj.fn, " ", sel_site, " -r ", as.character(4*N*(1 - dbinom(0,len_hap, r))), " ", as.character(len_hap), " -t ", as.character(4*N* (1 - dbinom(0, len_hap, u))), " -T -L > ", msout.fn, sep = "")
	system(ms.string)

	#cut down mssel output to just the part rent+ needs
	#process.string <- paste("sed '/positions: /,$!d' ", msout.fn, " | sed 's/.*positions: //' > ", rent_in_fn, sep = "")
	#system(process.string)

	#ms_haplotypes_list[[i]] <- readLines(rent_in_fn) #save the haplotypes produced by ms (and fed to rent+) in a list entry.
	#the result will be a character vector. The first entry in the vector is relative (0-1) positions, with the selected
	#site at 0.05. The remaining entries are string haplotypes of 0s and 1s.

	#pull trees from ms since rent+ hasn't been run
	msoutlines <- readLines(msout.fn)
	dists <- numeric(0)
	for(m in 1:(length(msoutlines) - 4)){
		suppressWarnings(dists[m] <- as.numeric(substring(strsplit(msoutlines[m + 4], "]")[[1]][1],2)))
	}
	dists <- dists[!is.na(dists) & dists != Inf]
	posits.rs <- cumsum(dists)
	ms_to_extract <- which.max(posits.rs[posits.rs < sel_site]) + 1
	if(length(ms_to_extract)  == 0){ms_to_extract <- 1}
	ms_trees_list[[i]] <- read.tree(text = msoutlines[4 + ms_to_extract])		
	if(i %% 10 == 0){
		print(paste("simulation locus", as.character(i), "of trait", as.character(k), "complete"))
	}
}


#Run through list and check for sites where 2+ variants had same coordinates as
#selected site. At those places, check which tree(s) are monophyletic for derived tips.
#if more than one, select randomly from among them.
for(i in 1:n.loci){
	if(is.null(names(ms_trees_list[[i]]))){
		cand.trees <- ms_trees_list[[i]]
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
	}
}



#split into ancestral and derived trees and retrieve coalescence times.
anc_trees_ms <- list()
der_trees_ms <- list()
ct_ms_all <- list()
ct_ms_anc <- list()
ct_ms_der <- list()


for(i in 1:n.loci){
	anc_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
	der_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
	anc_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], der_tips_ms)
	der_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], anc_tips_ms)
	ct_ms_all[[i]] <- coal.times(ms_trees_list[[i]])
	ct_ms_anc[[i]] <- coal.times(anc_trees_ms[[i]])
	ct_ms_der[[i]] <- coal.times(der_trees_ms[[i]])
	
}

print(paste("tree processing complete"))



times.c <- list()
for(i in 1:length(ms_trees_list)){
	times.c[[i]] <- trees_to_times(ms_trees_list[[i]], anc_trees_ms[[i]], der_trees_ms[[i]], time, sure.alt.is.derived = FALSE, units_in = 4)
}

lins.list <- list()
for(i in 1:length(ms_trees_list)){
	lins.list[[i]] <- times_to_lins(times.c[[i]], time)
}



trajs_neut <- matrix(nrow = length(time), ncol = length(ms_trees_list))
for(i in 1:length(ms_trees_list)){	
	trajs_neut[,i] <- est_af_traj_neut(lins.list[[i]])
}
trajs_neut[time == 0,] <- n_ders / n_chroms
traj.phen.neut <- 2 * trajs_neut %*%  eff_sizes 

qxtest_mat[k,1:3] <- Qx_test(trajs_neut[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat[k,4] <- Qx_test(trajs_neut[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
print(pars[k,])
print(qxtest_mat[k,])
print(aggregate(qxtest_mat[,4], by = pars[,c(1,5)], FUN = function(x){mean(x < .05, na.rm = TRUE)}  ))


fn_str <- paste("trees_loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")


save(list = c("time", "eff_sizes", "n_ders", "ms_trees_list"),file = paste(out_dir, fn_str, ".RData", sep = ""))




