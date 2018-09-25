

#Doc Edge, 4/10/18
#Goal: Simulate phenotypes (possibly under selection) that are additive functions
# of genotypes at many loci. This version bypasses fitting the ARG in rent+.
#Because we do not need to infer recombination, no recombination is specified.
#Because we only need


#Plan: simulate n_loci allele frequency trajectories.
#to do this, select an effect size from normal dist with mean 0 and sd
#equal to heritability over n_loci (from Martin et al). locus undergoes selection 
#according to effect size and selection gradient of trait (which changes through time). 
#Generate an allele frequency trajectory obeying the selection coefficient.
#Do this for all loci and save "phenotype" trajectory.
#Run through ms and save trees to estimate phenotype trajectory using coalescent times.
#Run haplotypes through RENT+ to see whether estimated trees recover phenotype trajectory.





#simulate one "true" phenotypic trajectory at set parameters.
#randomness across reps will be in sampling given moden freqs
#and in coalescent process.
trajs <- list()
eff_sizes <- numeric(0)
ss <- numeric(0)


#Simulate and write n.loci allele-frequency trajectories. None of these should have fixed,
#so reject and do it again if fixed. Optionally, reject if minor allele frequency is < 0.01.
for(i in 1:n.loci){
	fix <- 1
	while(fix == 1){
		eff_size <- rnorm(1,0,sqrt(herit/n.loci))
		s_locus <- eff_size * sel.intens / sd.trait #Charlesworth B3.7.7
		ss[i] <- s_locus
		p0 <- gen.neut.sfs(1, 2*N)
		sel.part <- sel.traj(p0, s = s_locus, N = N, t = t-t.off)
		if(t.off > 0){
			post.drift <- neut.traj.time(p0 = sel.part[nrow(sel.part),2], N, t.off)
			sel.part <- rbind(sel.part, cbind(post.drift[,1] + t - t.off, post.drift[,2] ))			
		}
		if(sel.part[nrow(sel.part),2] >= 0.01 & sel.part[nrow(sel.part),2] <= 0.99){fix = 0}
	}
	driftup <- neut.traj(p0, N, loss = TRUE)
	traj.fwd <- rbind(cbind(driftup[,1], rev(driftup[,2])), cbind(sel.part[-1,1] + max(driftup[,1]), sel.part[-1,2] ))
	traj <- traj.fwd
	traj[,2] <- rev(traj.fwd[,2])
	trajs[[i]] <- traj
	eff_sizes[i] <- eff_size 
}

#current allele frequencies.
curr.freqs <- sapply(trajs, function(x){x[1,2]})
trait.sd <- sqrt(sum(2 * curr.freqs * (1 - curr.freqs) *eff_sizes^2))
eff_sizes_unscaled <- eff_sizes
eff_sizes <- eff_sizes / trait.sd #The SD of the "genetic" part of the trait is set
#to be 1 in the present

#recover trajectory of phenotype over time
matrixify.list.of.trajs <- function(trajs){
	lns <- sapply(trajs, length)/2
	max.length <- max(lns)
	mat.ret <- matrix(0, nrow = max.length, ncol = length(trajs))
	for(i in 1:length(trajs)){
		mat.ret[1:lns[i],i] <- trajs[[i]][,2]
	}
	mat.ret
}

mat.trajs <- matrixify.list.of.trajs(trajs)

phen.traj <- as.numeric(2 * eff_sizes %*% t(mat.trajs))  
pt.time <- seq(0, by = 1/(2*N), length.out = max(sapply(trajs, length)/2) )




for(l in 1:length(n_chromss)){

n_chroms <- n_chromss[l]

n_ders <- numeric(n.loci)
ms_trees_list <- list()
#rent_trees_list <- list()

fn_str <- paste("loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")



#for each locus, write trajectory, run ms to simulate a sample and tree,
#and save the tree, as well
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

	#There's only one tree without recombination, and it's on the fifth line.
	true_sel_tree_newick <- readLines("ms_out.txt", n = 5)[5]

	ms_trees_list[[i]] <- read.tree(text = true_sel_tree_newick)
	print(paste("simulation round", as.character(i), "complete"))
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
	print(paste("round", as.character(i), "of tree processing complete"))
}

#time <- c(seq(0, .1, by = 0.001), seq(0.11, 1, by = 0.01), seq(1.02, 2, by = 0.02), seq(2.05, 4, by = 0.05))


times.c <- list()
for(i in 1:length(ms_trees_list)){
	times.c[[i]] <- trees_to_times(ms_trees_list[[i]], anc_trees_ms[[i]], der_trees_ms[[i]], time, sure.alt.is.derived = TRUE, units_in = 4)
}

lins.list <- list()
for(i in 1:length(ms_trees_list)){
	lins.list[[i]] <- times_to_lins(times.c[[i]], time)
}


save.image(paste(out_dir, fn_str, ".RData", sep = ""))

}


