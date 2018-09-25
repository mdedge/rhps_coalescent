

#Doc Edge, 6/13/18
#Goal: Simulate phenotypic scores (possibly under selection) that are additive functions
# of genotypes at many loci. 


#Plan: simulate n_loci allele frequency trajectories.
#to do this, select an effect size from normal dist with mean 0 and sd
#equal to heritability over n_loci (from Martin et al). locus undergoes selection 
#according to effect size and selection gradient of trait (which changes through time). 
#Generate an allele frequency trajectory obeying the selection coefficient.
#Do this for all loci and save "phenotype" trajectory.
#Run through ms and save trees to estimate phenotype trajectory using coalescent times.
#Run haplotypes through RENT+ to see whether estimated trees recover phenotype trajectory.


trajs <- list()
sel.parts <- list()
eff_sizes <- numeric(0)
ss <- numeric(0)
curr.freqs <- numeric(0)


#need a constant that reflects a harmonic series from 1 to 2N-1, but only the
#terms where i/(2N) \in [.01, .99]. for N large this will be very close to the 
#analogous integral, which is ln(.99) - ln(.01) ~= 4.595
ser <- 1:(2*N-1)
harm.const <- sum(1/ser[ser >= .01*2*N & ser <= .99*2*N])



#Simulate and write n.loci allele-frequency trajectories. None of these should have fixed,
#so reject and do it again if fixed. Reject if minor allele frequency is < 0.01.
shift.achieved <- 0
while(shift.achieved == 0){
	pre.sel.freqs <- numeric(0)
	post.sel.freqs	<- numeric(0)
	for(i in 1:n.loci){
		fix <- 1
		while(fix == 1){
			eff_size <- rnorm(1,0,sqrt(herit * sd.trait^2 * harm.const /n.loci))
			s_locus <- eff_size * sel.intens / sd.trait #Charlesworth B3.7.7
			ss[i] <- s_locus
			p0 <- gen.neut.sfs.condit(1, 2*N, .01) #select from neutral sfs conditional on common.
			pre.sel.freqs[i] <- p0 
			sel.part <- sel.traj(p0, s = s_locus, N = N, t = t-t.off)
			post.sel.freqs[i] <- sel.part[nrow(sel.part),2]
			if(t.off > 0){
				post.drift <- neut.traj.time(p0 = sel.part[nrow(sel.part),2], N, t.off)
				sel.part <- rbind(sel.part, cbind(post.drift[,1] + t - t.off, post.drift[,2] ))			
			}
			if(sel.part[nrow(sel.part),2] >= 0.01 & sel.part[nrow(sel.part),2] <= 0.99){fix = 0}
		}
		eff_sizes[i] <- eff_size 
		curr.freqs[i] <- sel.part[nrow(sel.part),2]
		sel.parts[[i]] <- sel.part
	}	
	#Check whether the achieved shift is close to the target.
	trait.sd <- sqrt(sum(2 * curr.freqs * (1 - curr.freqs) *eff_sizes^2))
	sel.shift <- sum(2 * eff_sizes * (post.sel.freqs - pre.sel.freqs)) / trait.sd
	target.shift <- (sel.intens * sd.trait) * (t - t.off) * 2 * N
	if(target.shift*.95 <= sel.shift & sel.shift <= target.shift*1.05){shift.achieved <- 1}
#	print("trait.attempted")
}


for(i in 1:n.loci){
	sel.part <- sel.parts[[i]]
	driftup <- neut.traj(pre.sel.freqs[i], N, loss = TRUE)
	traj.fwd <- rbind(cbind(driftup[,1], rev(driftup[,2])), cbind(sel.part[-1,1] + max(driftup[,1]), sel.part[-1,2] ))
	traj <- traj.fwd
	traj[,2] <- rev(traj.fwd[,2])
	trajs[[i]] <- traj
}

eff_sizes_unscaled <- eff_sizes
eff_sizes <- eff_sizes / trait.sd #The SD of the polygenic score is set
#to be 1 in the present


rm(ser)
rm(harm.const)
rm(eff_size)
rm(s_locus)
rm(p0)
rm(sel.part)
#rm(post.drift)
rm(post.sel.freqs)
rm(pre.sel.freqs)


mat.trajs <- matrixify.list.of.trajs(trajs)

phen.traj <- as.numeric(2 * eff_sizes %*% t(mat.trajs))  
pt.time <- seq(0, by = 1/(2*N), length.out = max(sapply(trajs, length)/2) )




n_ders <- numeric(n.loci)
ms_trees_list <- list()
rent_trees_list <- list()
ms_haplotypes_list <- list()


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
	process.string <- paste("sed '/positions: /,$!d' ", msout.fn, " | sed 's/.*positions: //' > ", rent_in_fn, sep = "")
	system(process.string)

	#run rent+
	rent.string <- paste("java -jar ", rentplus_fn, " -t -l ", as.character(len_hap), " ", rent_in_fn, " ", msout.fn, sep = "")
	system(rent.string)

	ms_haplotypes_list[[i]] <- readLines(rent_in_fn) #save the haplotypes produced by ms (and fed to rent+) in a list entry.
	#the result will be a character vector. The first entry in the vector is relative (0-1) positions, with the selected
	#site at 0.05. The remaining entries are string haplotypes of 0s and 1s.
	

	#Read in trees at selected site in rent+
	rent_trees_fn <- paste(rent_in_fn, ".trees", sep = "")
	rent_trees <- readLines(rent_trees_fn)
	#which tree is for the selected site?
	rent_sel_ind <- grep(paste(as.character(sel_site), "\t", sep = ""), rent_trees)

	#get tree at selected site in newick format, both true and rent-estimated.
	rent_sel_tree_newick <- sub(paste(as.character(sel_site), "\t", sep = ""), "", rent_trees[rent_sel_ind])
	rent_trees_list[[i]] <- read.tree(text = paste(rent_sel_tree_newick, ";", sep = ""))
	
	#pull trees from ms
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

	print(paste("simulation locus", as.character(i), "of trait", as.character(k), "complete"))
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
		#rent_trees_list[[i]] <- rent_trees_list[[i]][[ind.to.take]]
	}
}

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
		#ms_trees_list[[i]] <- cand.trees[[ind.to.take]]
		rent_trees_list[[i]] <- rent_trees_list[[i]][[ind.to.take]]
	}
}



#split into ancestral and derived trees and retrieve coalescence times.
anc_trees_ms <- list()
der_trees_ms <- list()
ct_ms_all <- list()
ct_ms_anc <- list()
ct_ms_der <- list()

anc_trees_rent <- list()
der_trees_rent <- list()
ct_rent_all <- list()
ct_rent_anc <- list()
ct_rent_der <- list()

for(i in 1:n.loci){
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
}




times.c <- list()
for(i in 1:length(ms_trees_list)){
	times.c[[i]] <- trees_to_times(ms_trees_list[[i]], anc_trees_ms[[i]], der_trees_ms[[i]], time, sure.alt.is.derived = FALSE, units_in = 4)
}

lins.list <- list()
for(i in 1:length(ms_trees_list)){
	lins.list[[i]] <- times_to_lins(times.c[[i]], time)
}


times.c.rent <- list()
for(i in 1:length(rent_trees_list)){
	times.c.rent[[i]] <- trees_to_times(rent_trees_list[[i]], anc_trees_rent[[i]], der_trees_rent[[i]], time, sure.alt.is.derived = FALSE, units_in = 2)
}

lins.list.rent <- list()
for(i in 1:length(rent_trees_list)){
	lins.list.rent[[i]] <- times_to_lins(times.c.rent[[i]], time)
}


fn_str <- paste("loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")


save.image(paste(out_dir, fn_str, ".RData", sep = ""))




