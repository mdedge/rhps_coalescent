

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

	msoutlines <- readLines(msout.fn)
	ms_trees_list[[i]] <- read.tree(text = msoutlines[5])		
	if(i %% 10 == 0){
		print(paste("simulation locus", as.character(i), "of daf", as.character(k), "complete"))
	}
}





#split into ancestral and derived trees.
anc_trees_ms <- list()
der_trees_ms <- list()

for(i in 1:n.loci){
	anc_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
	der_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
	anc_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], der_tips_ms)
	der_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], anc_tips_ms)
	
}

extract.tbl.sum <- function(x){sing.bl(branch.length.by.n.descendants(x))}


tbl.anc <- sapply(anc_trees_ms, extract.tbl.sum)
tbl.der <- sapply(der_trees_ms, extract.tbl.sum)

mat.tbl <- cbind(n_chroms - n_der, tbl.anc, n_der, tbl.der, tbl.anc / (n_chroms - n_der) - tbl.der / n_der)




print(paste("tree processing complete"))




fn_str <- paste("trees_nchr", as.character(n_chroms), "daf_", as.character(daf), sep = "")


save(list = c("n_ders", "ms_trees_list", "anc_trees_ms", "der_trees_ms", "mat.tbl"),file = paste(out_dir, fn_str, ".RData", sep = ""))




