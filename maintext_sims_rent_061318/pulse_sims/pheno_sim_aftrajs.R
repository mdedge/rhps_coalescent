

#Doc Edge, 6/13/18
#Goal:  simulate n_loci allele frequency trajectories.
#to do this, select an effect size from normal dist with mean 0 and sd
#equal to heritability over n_loci  locus undergoes selection 
#according to effect size and selection gradient of trait (which changes through time). 
#Generate an allele frequency trajectory obeying the selection coefficient.
#Do this for all loci and save "phenotype" trajectory.

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
	print("trait attempted")
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




