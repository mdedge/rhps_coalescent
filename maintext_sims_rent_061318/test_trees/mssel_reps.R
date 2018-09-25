
N <- 10000
herit <- 1

sel.intenses <- 0
n.loci <- 5000
n_chroms <- 200
t <- 0
t.off <- 0
phen_num <- 1 #1 number for each rep we want to do at each combination of parameters

helper_fn <- "../../helper_functions_coal_sel.R"
source(helper_fn) #read in helper functions and load ape package
traj.fn <- "temp.txt"

ms_dir <- "../../msseldir/"

len_hap <- 200000 #the length (in base pairs) of the haplotype -- short because 
#we are not worrying about recombination--just need the sel site.
sel_site <- 100000 #the position of the selected site in the haplotype
u <- 2e-8 #the neutral mutation rate per base pair/generation
r <- 2.5e-8 #the recombination rate per base pair/generation
options(scipen = 999) #disable scientific notation so that parameters passed to ms are given as numbers

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


pre.sel.freqs <- numeric(0)
post.sel.freqs	<- numeric(0)
for(i in 1:n.loci){
	fix <- 1
	while(fix == 1){
		#eff_size <- rnorm(1,0,sqrt(herit * 1^2 * harm.const /n.loci))
		#s_locus <- eff_size * sel.intens / sd.trait #Charlesworth B3.7.7			
		p0 <- gen.neut.sfs.condit(1, 2*N, .01) #select from neutral sfs conditional on common.
		pre.sel.freqs[i] <- p0 
		sel.part <- sel.traj(p0, s = 0, N = N, t = t-t.off)
		post.sel.freqs[i] <- sel.part[nrow(sel.part),2]
		if(t.off > 0){
			post.drift <- neut.traj.time(p0 = sel.part[nrow(sel.part),2], N, t.off)
			sel.part <- rbind(sel.part, cbind(post.drift[,1] + t - t.off, post.drift[,2] ))			
		}
		if(sel.part[nrow(sel.part),2] >= 0.01 & sel.part[nrow(sel.part),2] <= 0.99){fix = 0}
			}
		#eff_sizes[i] <- eff_size 
		curr.freqs[i] <- sel.part[nrow(sel.part),2]
		sel.parts[[i]] <- sel.part
}


for(i in 1:n.loci){
	sel.part <- sel.parts[[i]]
	driftup <- neut.traj(pre.sel.freqs[i], N, loss = TRUE)
	traj.fwd <- rbind(cbind(driftup[,1], rev(driftup[,2])), cbind(sel.part[-1,1] + max(driftup[,1]), sel.part[-1,2] ))
	traj <- traj.fwd
	traj[,2] <- rev(traj.fwd[,2])
	trajs[[i]] <- traj
	print(paste("locustraj", as.character(i)))
}



mat.trajs <- matrixify.list.of.trajs(trajs)
pt.time <- seq(0, by = 1/(2*N), length.out = max(sapply(trajs, length)/2) )



n_ders <- numeric(n.loci)
ms_trees_list <- list()


#for each locus, write trajectory, run ms to simulate a sample and tree,
#and run rent+ to infer tree. Save both the ms and the rent+ trees, as well
#as the number of derived alleles for each locus.
for(i in 1:n.loci){
	msout.fn <- paste("out2/msout", as.character(i), ".txt", sep = "")	
	write.traj(traj.fn, trajs[[i]])
	curr.freq.der <- curr.freqs[i]
	n_der <- rbinom(1, n_chroms, curr.freq.der) #number of chroms with derived allele.
	if(n_der == 0){n_der <- n_der + 1}
	if(n_der == n_chroms){n_der <- n_der - 1}
	n_ders[i] <- n_der
	#run mssel
	ms.string <- paste(ms_dir, "mssel ", as.character(n_chroms), " 1 ", as.character(n_chroms - n_der), " ", as.character(n_der), " ", traj.fn, " ", sel_site, " -r ", as.character(4*N*(1 - dbinom(0,len_hap, r))), " ", as.character(len_hap), " -t ", as.character(4*N* (1 - dbinom(0, len_hap, u))), " -T -L > ", msout.fn, sep = "")
	system(ms.string)
	
	#pull trees from ms
	msoutlines <- readLines(msout.fn)
	dists <- numeric(0)
	for(m in 1:(length(msoutlines) - 4)){
		suppressWarnings(dists[m] <- as.numeric(substring(strsplit(msoutlines[m + 4], "]")[[1]][1],2)))
	}
	dists <- dists[!is.na(dists) & dists != Inf]
	posits.rs <- cumsum(dists)
	ms_to_extract <- which.max(posits.rs[posits.rs < sel_site]) + 1
	ms_trees_list[[i]] <- read.tree(text = msoutlines[4 + ms_to_extract])		

	print(paste("simulation locus", as.character(i), "complete"))
}


#Run through list and check for sites where 2+ variants had same coordinates as
#selected site. At those places, check which tree(s) are monophyletic for derived tips.
#if more than one, select randomly from among them.
for(i in 1:n.loci){
#for(i in 1:1600){
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

#check whether derived tips are monophyletic.
mono.ders <- numeric(0)
mono.ancs <- numeric(0)
#for(i in 1:1600){
for(i in 1:n.loci){
	der_tips <- which(ms_trees_list[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
	mono.ders[i] <- is.monophyletic(ms_trees_list[[i]], der_tips)
	anc_tips <- which(ms_trees_list[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
	mono.ancs[i] <- is.monophyletic(ms_trees_list[[i]], anc_tips)
}



