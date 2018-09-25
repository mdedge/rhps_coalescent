#3/20/18 by Doc Edge
#Goal: Assemble functions useful for project on coalescent approcahes to polygenic selection.
#Many of these functions relate to working with trees stored as package ape phylo-type objects.

#package ape is needed for representing sample trees.
if(!("ape" %in% installed.packages())){install.packages("ape")}
library(ape)

#generate derived allele frequences for n.loci loci under neutral sfs.
#n can either be size of sample of chromosomes (for sample SFS)
#or size of population (in # of chromosomes) for population sfs
gen.neut.sfs <- function(n.loci, n){
  probs.uw <- 1/(1:(n-1))
  probs <- probs.uw/sum(probs.uw)
  cdf <- cumsum(probs)
  percs <- runif(n.loci, 0 , 1)
  get.cdfind <- function(perc){
    (sum(cdf <= perc) + 1)/n 
  }
  sapply(percs, get.cdfind)
}

#generate derived allele frequences for n.loci loci under neutral sfs.
#n can either be size of sample of chromosomes (for sample SFS)
#or size of population (in # of chromosomes) for population sfs
#condition.freq is a decimal such that all returned allele frequencies
#should be between condition.freq and 1 - condition.freq
gen.neut.sfs.condit <- function(n.loci, n, condition.freq = 0){
  ser <- 1:(n-1)
  ser <- ser[ser/n >= condition.freq & ser/n <= 1 - condition.freq]
  probs.uw <- 1/ser
  probs <- probs.uw/sum(probs.uw)
  cdf <- cumsum(probs)
  percs <- runif(n.loci, 0 , 1)
  get.cdfind <- function(perc){
    ser[(sum(cdf <= perc) + 1)]/n 
  }
  sapply(percs, get.cdfind)
}

#recover trajectory of polygenic score over time
matrixify.list.of.trajs <- function(trajs){
	lns <- sapply(trajs, length)/2
	max.length <- max(lns)
	mat.ret <- matrix(0, nrow = max.length, ncol = length(trajs))
	for(i in 1:length(trajs)){
		mat.ret[1:lns[i],i] <- trajs[[i]][,2]
	}
	mat.ret
}


#simulates allele frequency trajectory forward in time conditional on loss.
#returns a two-column matrix; left column is time, right is derived allele frequency.
#Uses Lee & Coop eq. A.14-15 (biorxiv version)
neut.traj.loss <- function(p0, N = 10000, delta.t = 1/(2*N)){
	traj <- p0
	time <- 0
	time.ct <- time
	while(traj[length(traj)] != 1 & traj[length(traj)] != 0){
    		curr <- traj[length(traj)]
    		nex.freq <- rnorm(1, curr-curr*delta.t, sqrt(delta.t) * sqrt(curr*(1-curr)) )
    		if(nex.freq > 1){nex.freq <- 1}
    		if(nex.freq < 0){nex.freq <- 0}
    		traj <- c(traj, nex.freq)
    		time <- time + delta.t
    		time.ct <- c(time.ct, time)
  	}
	if(traj[length(traj)] == 0){
		return(cbind(time.ct, traj))	
	}
	neut.traj.loss(p0, N, delta.t)
}



#Function to generate neutral allele frequency trajectory forward in time.
#Trajectory can be reversed to give trajectory up to a given frequency.
#Starts at frequency p0. Frequency at next time step is drawn from a
#normal with expectation equal to current frequency (curr) and variance equal
#to delta.t*curr*(1-curr), where delta.t is the size of the time step.
#(The only use of N is to scale the default time step size to 1/(2*N) coalescent
#units, or one generation under standard model). The default of loss = NULL
#returns the first trajectory generated, regardless of whether it ends in
#loss or fixation. If loss == TRUE, then the function runs conditional on loss.
# If loss == FALSE, then the function
#runs diffusion conditional on fixation.
neut.traj <- function(p0, N = 10000, loss = NULL, delta.t = 1/(2*N)){
  if(is.null(loss)){
    traj <- p0
    time <- 0
    time.ct <- time
    while(traj[length(traj)] != 1 & traj[length(traj)] != 0){
      curr <- traj[length(traj)]
      nex.freq <- rnorm(1, curr, sqrt(delta.t) * sqrt(curr*(1-curr)) )
      if(nex.freq > 1){nex.freq <- 1}
      if(nex.freq < 0){nex.freq <- 0}
      traj <- c(traj, nex.freq)
      time <- time + delta.t
      time.ct <- c(time.ct, time)
    }
    return(cbind(time.ct, traj))
  }
  if(loss == TRUE){
    return(neut.traj.loss(p0, N, delta.t))
  }
  if(loss == FALSE){
    #to simulate an allele that fixes, simulate an allele that's lost
    #and flip the allele frequencies.
    traj.anc <-  neut.traj.loss(1-p0, N, delta.t)
    return(cbind(traj.anc[,1], 1 - traj.anc[,2]))
  }
}


#function to create drifting allele frequency for a set amount of time, starting
#at frequency p0. This differs from neut.traj because neut.traj simulates
#to either loss or fixation. Simulates forward in time, so reverse time if a 
#back-in-time trajectory is desired. p0 is starting frequency, t is time (in units
#of 2N generations), N is population size, and delta.t is the time interval to use.
neut.traj.time <- function(p0, N = 10000, t = 0.1, delta.t = 1/(2*N)){
  traj <- p0
  time <- 0
  time.ct <- time
  while(time <= t & !isTRUE(all.equal(time, t))){
    curr <- traj[length(traj)]
    if(curr > 0 & curr < 1){
      nex.freq <- rnorm(1, curr, sqrt(delta.t) * sqrt(curr*(1-curr)) ) 
      if(nex.freq > 1){nex.freq <- 1}
      if(nex.freq < 0){nex.freq <- 0}
    }else{nex.freq <- curr}
    traj <- c(traj, nex.freq)
    time <- time + delta.t
    time.ct <- c(time.ct, time)
  }
  cbind(time.ct, traj)
}


#Generate an allele frequency trajectory forward in time under selection.
#p0 is the frequency of the allele when selection starts.
#s is the selection coefficient, measured as the difference between the heterozygote 
#and one of the homozygotes, which is
#why the infinitesimal mean here is Nsp(1-p) rather than 2Nsp(1-p).
#N is the (constant, diploid) population size.
#t is the time (in units of 2N gnerations) over which the allele frequency trajectory should be computed.
#delta.t is the interval between time steps.
sel.traj <- function(p0, s = 0.001, N = 10000, t = 0.1, delta.t = 1/(2*N)){
  traj <- p0
  time <- 0
  time.ct <- time
  while(time <= t & !isTRUE(all.equal(time, t))){
    curr <- traj[length(traj)]
    if(curr > 0 & curr < 1){
      next.freq <- rnorm(1,  curr + delta.t*2*N*s*curr*(1-curr), sqrt(delta.t) * sqrt(curr*(1-curr)) ) 
      if(next.freq > 1){next.freq <- 1}
      if(next.freq < 0){next.freq <- 0}
    }else{next.freq <- curr}
    traj <- c(traj, next.freq)
    time <- time + delta.t
    time.ct <- c(time.ct, time)
  }
  cbind(time.ct, traj)
}

#tests
#st <- sel.traj(0.05, s = 0.001, t = 0.5)
#plot(st[,1], st[,2], ylim = c(0,1), type = "l")
#summary(st)


#Generate trajectory involving drift to specified frequency
#followed by selection for fixed time period t
#if t = 0, then only a neutral trajectory generated.
#neutral trajectory is reversed so that selection starts at p0.
#derived = TRUE or FALSE allows considitioning on whether selected allele is derived.
#startat0 FALSE makes the neutral times negative; selection starts at time 0.
#startat0 TRUE makes time start when the locus first becomes dimorphic.
#Trajectory is meant to run BACKWARD in time.
neut.then.sel.traj <- function(p0, s = 0.001, t = 0.1, N = 10000, derived = NULL, delta.t = 1/(2*N)){
  neut.part <- neut.traj(p0, N, loss = derived, delta.t)
  neut.part <- neut.part[-1,] #remove first row because will be duplicated in selected part
  neut.part <- neut.part[rev(1:nrow(neut.part)),]
  neut.part[,1] <- -neut.part[,1]
  traj <- neut.part
  if(t > 0){
    sel.part <- sel.traj(p0, s, N, t, delta.t)
    traj <- rbind(neut.part, sel.part)
  }
  traj[,1] <- traj[,1] - min(traj[,1])
  traj[,2] <- rev(traj[,2]) #makes it run backward in time.
  traj
}


#Write an allele-frequency trajectory to a file that can be read
#by mssel.
#mssel file tracks the derived allele frequency. This function writes down the derived allele frequency
#and returns both the starting frequency and whether the derived or ancestral allele experienced
#selection. (Probability the ancestral allele is the one to experience selection is p0 if derived
#is set to NULL.)
write.trajfile <- function(filename, p0 = 0.1, s = 0.001, t = 0.1, N = 10000, derived = NULL, delta.t = 1/(2*N)){
  traj <- neut.then.sel.traj(p0, s, t, N, derived, delta.t)
  s.allele <- "derived"
  if(traj[nrow(traj),2] == 1){
    s.allele = "ancestral"
    traj[,2] <- 1 - traj[,2]
  }
  traj[,1] <- traj[,1]/2 #mssel measures time in units of 4N generations, not 2N.
  comment.line <- paste("#parameters: p0=", as.character(p0), ", s=", as.character(s), ", t=",
                        as.character(t), ", N=", as.character(N), ", derived=", as.character(derived),
                        ", delta.t=", as.character(delta.t), sep = "")
  len.traj <- nrow(traj)
  head.lines <- paste("ntraj: 1\nnpop: 1\nn: ", as.character(len.traj), sep = "")
  fileConnection <- file(filename)
  writeLines(c(comment.line, head.lines), fileConnection)
  close(fileConnection)
  write(t(traj), filename, ncolumns = 2, append = TRUE)
  c(traj[1,2] , s.allele, max(traj[,1])) 
}


#Somewhat redundant with write.trajfile()--should merge these
#This function takes an allele frequency trajectory and writes it to a file,
#whereas write.trajfile simulates the trajectory *and* writes to a file.
#the default in.time = 2 assumes the input trajectory is in units of 2N generations.
#The default out.time says to convert to units of 4N generations (which ms uses).
#the first column of the traj object is assumed to contain times, and the second
#contains allele frequencies for the derived allele.
write.traj <- function(filename, traj, in.time = 2, out.time = 4){
	traj[,1] <- traj[,1] * (in.time / out.time)
	len.traj <- nrow(traj)
  	head.lines <- paste("ntraj: 1\nnpop: 1\nn: ", as.character(len.traj), sep = "")
  	fileConnection <- file(filename)
	writeLines(head.lines, fileConnection)
	close(fileConnection)
  	write(t(traj), filename, ncolumns = 2, append = TRUE)
}



#Given a filename of an mssel output file (fn_msout), and the length of the
#haplotype simulated by ms, writes a "sites" format file that can be read
#into argweaver to filename fn_sites.
convert_ms_out_to_sites <- function(fn_msout, fn_sites, length.hap = 200000){
	#read in the file and identify the haplotypes	
	mslines <- readLines(fn_msout)
	posline <- pmatch("positions", mslines)	
	haps <- mslines[(posline + 1):length(mslines)]
	n.haps <- length(haps)
	n.sites <- nchar(haps[1])	
	#the (arbitrary) haplotype names are the first line.	
	names.str <- "NAMES"
	for(i in 1:n.haps){
		names.str <- paste(names.str, "\tn", as.character(i), sep = "") 	
	}
	#Region string is the second line of the output	
	region.str <- paste("REGION\tchr\t1\t", as.character(length.hap), sep = "")
	#To get the haplotypes in the right format, turn them into a matrix of characters,
	#change 0 to A and 1 to T, and then put the characters together so that each
	#entry in the resulting vector is data for a site (rather than a haplotype).
	#Could speed this up, but it's not close to being the limiting factor in an 
	#argweaver simulation.	
	hap.mat <- matrix("", nrow = n.sites, ncol = n.haps)	
	for(i in 1:n.haps){
		hap.mat[,i] <- strsplit(haps[i], "")[[1]]
	}
	hap.mat[hap.mat == "0"] <- "A"
	hap.mat[hap.mat == "1"] <- "T"
	per.locus.data <- character(n.sites)	
	for(i in 1:n.sites){
		per.locus.data[i] <- paste(hap.mat[i,], sep = "", collapse = "")
	}
	pos <- as.numeric(strsplit(mslines[posline], " ")[[1]][-1])
	pos <- round(length.hap * pos)
	#Dealing with duplicated positions by shifting one of them to a nearby spot.
	#This is somewhat ad hoc and could be improved.
	l <- 1
	while(any(duplicated(pos)) & l <= 10){
		pos[duplicated(pos)] <- pos[duplicated(pos)] + round(min(diff(pos[!duplicated(pos)]))/2 + 0.01)
		l <- l + 1
	}
	#Assemble all the lines in the output file and write them to a connection.	
	lines <- c(names.str, region.str)
	for(i in 1:n.sites){
		lines[i + 2] <- paste(pos[i], per.locus.data[i], sep = "\t")
	}
	file.create(fn_sites)	
	fileCon <- file(fn_sites)
	writeLines(lines, fileCon)
	close(fileCon)
}


#Takes in an argweaver .smc format output file and reads in the
#tree for the target site (which is at a specified location, sel.pos).
#The resulting tree is returned as a phylo object.
read_smc <- function(fn_smc, sel.pos){
	smclines <- readLines(fn_smc)
	treelines <- smclines[seq(3, length(smclines), 2)]
	#start position of each tree
	starts <- as.numeric(sapply(strsplit(treelines, "\t"), "[[", 2))
	sel_tree_ind <- which.max(starts[starts <= sel.pos])
	key_tree <- strsplit(treelines[sel_tree_ind], "\t")[[1]][4]	
	names.perm <- as.numeric(strsplit(smclines[1], "\tn")[[1]][-1])
	key_phylo <- read.tree(text = key_tree)
	#ARGw messes with the ordering of the tips---this step changes the labels so they correspond
	#to the original 1:n.chrom labels.
	key_phylo$tip.label <- as.character(names.perm[as.numeric(key_phylo$tip.label) + 1])
	key_phylo
}


#Wrapper for read_smc() that reads in a bunch of argweaver outputs and returns them in a list.
#prefix is the prefix of the argweaver output .smc files, inds are the samples for which
#we want to extract the key tree and store in a list, and sel.pos is the position (in base pairs)
#of the key site.
read_specified_smcs <- function(prefix, inds = c(800,840,880,920,960,1000), sel.pos){
	toreturn <- list()
	for(i in 1:length(inds)){
		fn <- paste(prefix, ".", as.character(inds[i]), ".smc.gz", sep = "")
		toreturn[[i]] <- read_smc(fn, sel.pos)
	}
	toreturn
}


#Function to obtain descendant tips from internal node of a phylo object named tree,
#from Liam Revell's blog.
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters > length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

#Get the number of descendant tips from a given node. If node is a tip, return 1. tree is a phylo object.
getNDescendantTips <- function(tree, node){
   ntips <- sum(getDescendants(tree, node) <= length(tree$edge.length) - tree$Nnode + 1)
   if(ntips == 0){ntips <- 1}
   ntips
}


#Sum of branch lengths according to number of descendant tips per branch. Conditional on tree,
#this is related by a constant to expected unfolded SFS.
#tree is an ape-style phylo object.
branch.length.by.n.descendants <- function(tree){
   ntips <- length(tree$tip.label)
   ndescs.nodes <- sapply(1:(ntips + tree$Nnode), getNDescendantTips, tree = tree)
   ndescs.edges <- ndescs.nodes[tree$edge[,2]] 
   bl <- aggregate(tree$edge.length, list(ndescs.edges), FUN = sum)
   empty.fs <- (1:(ntips-1))[!((1:(ntips-1)) %in% bl[,1])]
   bl <- rbind(cbind(bl$Group.1, bl$x), cbind(empty.fs, rep(0, length(empty.fs)) ))
   if(nrow(bl) >= 2){bl <- bl[order(bl[,1]),]}
   colnames(bl) <- c("n.desc", "br.ln")
   bl
}


#total branch length of a tree (phylo object)
total.bl <- function(tree){
   sum(tree$edge.length)
}


#function to compute the Watterson constant, or the (k-1)th harmonic number
Watterson_cons <- function(k){
   vec <- 1:(k-1)
   sum(1/vec)
}

#given a matrix with branch lengths by # of descendants, return sum length of terminal branches
sing.bl <- function(eSFS){
   eSFS[eSFS[,1] == 1,2]
}

#given a matrix with branch lengths by # of descendants, return sum length of 
#branches with num descendants.
num.bl <- function(eSFS, num = 1){
   eSFS[eSFS[,1] == num,2]
}


#Function to obtain total time from present to an
#internal node. tree is an ape-style phylo object
coal.time <- function(tree, node, time = 0){
	desc <- getDescendants(tree, node)
	if(length(desc) < 2){
		print(paste("Warning: node provided (", as.character(node),") is external.", sep = ""))
		return(NA)
	}
	brnch <- which(tree$edge[,1] == node)[1]
	nodenext <- tree$edge[brnch,2]
	bl <- tree$edge.length[brnch]
	time <- time + bl
	if(getNDescendantTips(tree, nodenext) == 1){return(time)}
	coal.time(tree, nodenext, time)		
}


#function to return the coalescent times (always starting from 0)
#input is an ape-style phylo object.
coal.times <- function(tree){
	if(nrow(tree$edge) == 1){return(numeric(0))}
	tree <- collapse.singles(tree) #singleton nodes get deleted.	
	int.nodes <- unique(tree$edge[,1])
	sort(sapply(int.nodes, coal.time, tree = tree, time = 0))
}


#Return the time at which two subtrees are joined.
#input is coalescence time vectors for the full tree (allt) and two subtrees (anct and dert).
#the mut.time argument is an option to override this function
#with a user-provided number.
tree_join_time <- function(allt, anct, dert, mut.time = NULL){
	if(!is.null(mut.time)){return(mut.time)}	
	subs <- sort(c(anct,dert))
	bools <- allt[-length(allt)] == subs
	if(mean(bools) == 1){
		return(allt[length(allt)])
	}
	allt[!bools][1]
}


#takes a set of coalescence times from the whole tree and from a reference
#and alternative tree. Adds the mutation to the topmost branch of either tree.
#(By "adds a mutation", I mean that it adds a coalescence time to the vector, as if the derived type
#is coalescing from one lineage down to zero.)
#Assumes that the derived allele appeared between the coalescence of the derived tree and the
#joining of the two trees. If the two subtrees coalesce before they join, then the alternative
#tree is assumed to be derived.
#If sure.alt.is.derived is set to true, then the alt is forced to be the derived allele
#even if the ref tree seems to coalesce into the alt tree. This is done by changing the 
#overall tree so that the trees join just earlier in time than the alt tree coalesces.
#place determines where the mutation is placed on the topmost branch, relative to its length
#place = 0.5 (the default) puts it in the middle, lower numbers are closer to the present
add_mutation <- function(allt, reft, altt, tj, times, sure.alt.is.derived = FALSE, place = 0.5){
	#If both trees have at least two leaves (and  at least one coalescence time)	
	if(length(reft) > 0 & length(altt) > 0){
		#If we are sure that the alternative tree is derived but the
		#tree-join time is more recent than the final
		#coalescence in the alternative tree, then change the
		#tree-join time to be just earlier than the last coalescence
		#in the alternative tree.		
		if(sure.alt.is.derived & max(altt) > tj){
			tj <- max(altt) + 0.001
		}
		if(max(altt) > tj){ #if the alt tree hasn't coalesced when the trees join
			reft <- c(reft, (max(reft)*(1-place) + tj*place)) 
			altt <- c(altt, max(times) + 1)
		}
		if(max(altt) <= tj){
			reft <- c(reft, max(times) + 1)
			altt <- c(altt, (max(altt)*(1-place) + tj*place))
		}
	}
	if(length(altt) == 0){
		reft <- c(reft, max(times) + 1)		
		altt <- tj*place
	}
	if(length(reft) == 0 & !sure.alt.is.derived){
		reft <- tj*place		
		altt <- c(altt, max(times) + 1)	
	}
	if(length(reft) == 0 & sure.alt.is.derived){
		reft <- max(times) + 1		
		altt <- altt #We are working with one fewer coalescence than in the other
			     #cases, but we don't know when that coalescence was.
	}
	list(allt, reft, altt)
}






#Function to compute the value of population-size integral for adjusting
#coalescent times. Given a starting population size N0, a population size at time
#t Nt, a time t over which to compute the integral, a ploidy, and a functional
#form for the population-size change from time 0 to time t, computes the integral
#\int_0^t 1/(ploidy*N(u)) du, which answers the question, "if the population
#changes from N0 to Nt over t generations, how many coalescent-units have passed?".
#The form argument has three possible values, "linear", in which population size
#changes linearly from N0 to Nt; "piecewise.constant" in which population size
#switches from constant N0 to constant Nt at time t*frac.swit, and "exponential",
#in which the population size grows/decays exponentially from N0 to Nt in time t.
#frac.swit only applies to the piecewise constant functional form, which switches
#from N0 to Nt after t*frac.swit of the time has passed.
#func.user allows the user to input his/her own function giving population size
#at time u. If it is non-null, func.user overrides the form argument.
time.scaling <- function(N0, Nt, t, form = "linear", ploidy = 1, frac.swit = 0.5, func.user = NULL){
	if(form == "linear"){
		Nu <- function(u){N0 + u*(Nt-N0)/t}
	}
	if(form == "piecewise.constant"){
		Nu <- function(u){
			N0 + (Nt - N0)*(u >= frac.swit*t)
		}
	}
	if(form == "exponential"){
		Nu <- function(u){N0*((Nt/N0)^(1/t))^u}
	}
	if(!is.null(func.user)){Nu <- func.user}
	integrand <- function(u){1/(ploidy * Nu(u))}
	integrate(integrand, 0, t)[[1]]
}

#define the rising and falling factorial functions, which are used to compute the exact
#expectation for the number of lineages.
rising.factorial <- function(n, i){
	factorial(n+i-1)/factorial(n-1)
}
falling.factorial <- function(n, i){
	factorial(n)/factorial(n-i)
}



#compute the expected number of lineages ancestral to a
#haploid sample of size m assuming that the haploid population size
#changes from N0 to Nt over time t (in generations) into the past. (N0 is the present).
#There are three options for the expression argument, "exact", which gives
#the exact result from Tavare (1984) and given in Chen & Chen, "Maruvka",
#which gives the approximation of Maruvka et al., and "Slatkin", which
#gives the approximation of Slatkin & Rannala. 
#The form, frac.swit, and func.user arguments are passed directly to 
#time.scaling().
expec.ancs <- function(m, N0, Nt, t, expression = "Maruvka", form = "linear", frac.swit = 0.5, func.user = NULL){
	t.eff <- time.scaling(N0, Nt, t, form, ploidy = 1, frac.swit, func.user)
	if(expression == "Maruvka"){
		return(m/(m + (1-m)*exp(-t.eff/2)))
	}	
	if(expression == "Slatkin"){
		return(m/(1 + m*t.eff/2))
	}
	if(expression == "exact"){
		l <- 1:m		
		return(sum( (2*l-1) * (falling.factorial(m,l)/rising.factorial(m,l)) * exp(-(l*(l-1)/2) * t.eff) ) )
	}
}



#Computes equations 19 and 20 of Chen and Chen., which give the expec and variance
#of a normal approximation to the distribution of the number of lineages
#ancestral to a sample of m in the present at time t. Assumes haploid N0 and Nt.
#m is the sample size in the present (should be large), t is in generations.
#N0 is pop size now, Nt is pop time t in the past. Other arguments are passed
#directly to time.scaling() and control the population-size change from N0 to Nt.
ancs.approx.norm.mv <- function(m, N0, Nt, t, form = "linear", frac.swit = 0.5, func.user = NULL){
	gt <- time.scaling(N0, Nt, t, form, ploidy = 1, frac.swit, func.user)
	alpha <- m*gt/2
	beta <- gt/2
	eta <- (alpha * beta)/(alpha*(exp(beta) - 1) + beta*exp(beta))
	norm.mean <- 2*eta/gt
	norm.var <- (2*eta/gt)*((eta + beta)^2)*(1 + eta/(eta + beta) - eta/alpha - eta/(alpha + beta) - 2*eta)/beta^2
	c(norm.mean, norm.var)
}

#Computes the expectation of the approximately poisson number of coalescences
#that occur when m is large and t is very short. This appears in Chen & Chen just after 
#Equation 23.
ncoals.approx.poisson.m <- function(m, N0, Nt, t, form = "linear", frac.swit = 0.5, func.user = NULL){
	gt <- time.scaling(N0, Nt, t, form, ploidy = 1, frac.swit, func.user)
	0.5*m*(m-1)*gt
}

#gives the exact probability that there are mt lineages ancestral to a 
#sample of m0 chromosomes after t generations in the past. N0 is the 
#haploid pop size in the present, Nt at time t, and population size
#changed according to form argument.
ancs.exact.dist <- function(m0, mt, N0, Nt, t, form = "linear", frac.swit = 0.5, func.user = NULL){
	gt <- time.scaling(N0, Nt, t, form, ploidy = 1, frac.swit, func.user)
	i <- mt:m0
	sum( (((-1)^(i-mt)) * (2*i-1) * rising.factorial(mt, i-1) * falling.factorial(m0, i)/(factorial(mt) * factorial(i-mt) * rising.factorial(m0, i) ) ) * exp( -(i*(i-1))*gt/(2) )   )
}


#return the probability that m0 lineages coalesced down to mt lineages given
#that the haploid population size changed from N0 to Nt in t generations according
#to a specified trajectory. If the number of lineages at the start is small
#(<=20), then the exact probability is computed. Otherwise, if the normal
#approximation looks good (i.e. late enough in coal that right-truncation
#is not severe, meaning here that the maximum possible # of lineages is
#more than 3 sds above the expectation), then use it. If it is too early
#in the coalescent for the normal approximation, then use the poisson approximation.
p.ancs <- function(m0, mt, N0, Nt, t, form = "linear", frac.swit = 0.5, func.user = NULL, logp = FALSE){
	if(m0 <= 10){
		p.ex <- ancs.exact.dist(m0, mt, N0, Nt, t, form, frac.swit, func.user)
		if(logp == TRUE){return(log(p.ex))}
		return(p.ex)
	}
	mv <- ancs.approx.norm.mv(m0, N0, Nt, t, form, frac.swit, func.user)	
	expec.n <- mv[1]
	sd.n <- sqrt(mv[2])
	if((m0 - expec.n)/sd.n > 3){
		return( dnorm(mt, expec.n, sd.n, logp) )
	}
	expec.p <- ncoals.approx.poisson.m(m0, N0, Nt, t, form, frac.swit, func.user)
	dpois(m0 - mt, expec.p, logp)
}







#The scaled time expected before each of the n-1 coalescent events from n lineages
#down to 1. 
coal.exp.scaledtime <- function(n){
	m <- (n-1):1
	2*(1/m - 1/n)
}

#Get a method-of-moments estimate of the coalescent-units scaled time
#that has elapsed as a sample of n0 coalesces to nt lineages.
#Solves the Maruvka approximation of the expected # of ancestral lineages.
#If the sample has coalesced to 1 lineage, assume that it's been 2*(1-1/n0) 
#coalescent units, the expectation.
mom.scaledtime <- function(nt, n0){
	if(nt == 1){return( 2*(1-1/n0) )}	
	-2*log(n0*(1/nt-1)/(1-n0))
}

#Given a starting y value, ending y value, and vector of x coordinates,
#returns a vector of y coordinates along a line at requested x coordinates.
smooth.rise <- function(start, end, x){
	run <- x[length(x)] - x[1]
	slope = (end - start)/(run)
	start + slope*(x - x[1])
}

#"smooths out" a vector of step-like function values y of an independent
#variable x. If y[i]-y[i-1]==0, then the next time that a rise occurs, it
#is split up at an even slope, so there are no flat "step-like" portions.
smooth.steps <- function(x, y){
	chg_pts <- y[1:(length(y)-1)] != y[2:length(y)]
	smooth.y <- numeric(length(y))
	last_chg <- 0
	for(i in 1:length(chg_pts)){
		if(chg_pts[i] == TRUE){
			smooth.y[(last_chg+1):(i+1)] <- smooth.rise(y[last_chg+1], y[i+1], x[(last_chg+1):(i+1)])
			last_chg <- i
		}
	}
	chg <- which(y[1:(length(y)-1)] != y[2:length(y)])
	last.start.ind <- chg[length(chg)-1] + 1
	if(length(last.start.ind) == 0){last.start.ind <- 1} #deals with a bug that arises when lineage starts with only 2 copies
	last.end.ind <- chg[length(chg)] + 1
	last.slope <- (y[last.end.ind] - y[last.start.ind])/(x[last.end.ind] - x[last.start.ind])
	if(last.end.ind<length(y)){
		smooth.y[(last.end.ind + 1):length(y)] <- y[last.end.ind] + last.slope*(x[(last.end.ind + 1):length(y)] - x[last.end.ind])
	}
	smooth.y
}




#computes approximate variance of p estimate based on number of active lineages
#after t generations
#parameters are
#n: # of lineages of type 1 sampled
#m: # of lineages of type 2 sampled
#N: population size of type 1 alleles
#M: population size of type 2 alleles (true value of p is N/(N+M))
#t: time elapsed in generations
#NOTE: The function gives the same answers if N and M are decimals adding to 1 (p and 1-p),
#and t is given in coalescent units.
#var_phat_lineages <- function(n, m, N, M, t){
#	num1 <- exp(-t/(2*M))*(-((m-1)^2)*M + exp(t/M)*M*m^2 + exp(t/(2*M))*(M-2*m*M-(m-1)*m*t))
#	den1 <- (m-1)*m*M*(( log(m/(m-1)) - log(exp(t/(2*M))*m/(m-1)))^4)
#	num2 <- exp(-t/(2*N))*(-((n-1)^2)*N + exp(t/N)*N*n^2 + exp(t/(2*N))*(N-2*n*N-(n-1)*n*t))
#	den2 <- (n-1)*n*N*(( log(n/(n-1)) - log(exp(t/(2*N))*n/(n-1)))^4)
#	num3 <- num2
#	num4 <- 2*num2
#	den3 <- den2*N^2
#	den4 <- den2*N*(M+N)
#	(N*N*t*t/(4*(M+N)^2))*((num1/den1 + num2/den2)/((M+N)^2) + num3/den3 - num4/den4)
#}

#Idea: choose tau (based on current estimated p and N) that will minimize variance of estimate.
#could penalize for longer lengths of time to penalize the unrealistic assumption that
#everything stays the same for a long time.


#computes approximate variance of p estimate based on number of active lineages
#after tau coalescent time (in the whole sample)
#parameters are
#n: # of lineages of type 1 sampled
#m: # of lineages of type 2 sampled
#t: time elapsed (generations)
#N: population size of type 1 alleles
#M: population size of type 2 alleles (p is N/(N+M))
#Calculates:
#nt: # of lineages of type 1 remaining at the ancient end of interval
#mt: # of lineages of type 2 remaining at the ancient end of interval
#p: frequency of type 1
#NOTE: The function gives the same answers if N and M are 
#replaced by decimals adding to 1 (p and 1-p),
#and t is given in coalescent units.
#NOTE: This function is equivalent to the other
# var_phat_lineages, but it's easier to relate to the expression shown in
#the text. (previous expression came from mathematica)
var_phat_lineages <- function(n, m, N, M, t){
	p <- N/(N+M)
	tau <- t / (N+M)
	nt <- n/(n + (1-n) * exp(-t/(2*N)))
	mt <- m/(m + (1-m) * exp(-t/(2*M)))	
	if(n == 1 & m == 1){return(1/12)} #if there is one lineage left of either type, return 1/12 
					#variance of a uniform(0,1) RV.	
	if(N == 0 | M == 0 | is.na(N) | is.na(M)){return(0)}	#If allele has fixed or lost, return 0.
	fun.n <- 1/nt + 1/(nt-1) - 1/n - 1/(n-1) - tau/p
	fun.m <- 1/mt + 1/(mt-1) - 1/m - 1/(m-1) - tau/(1-p)
	var.est <- (4*(p^2)*((1-p)^2)/(tau^2))*((p^2)*fun.n + ((1-p)^2)*fun.m)
	#If the estimated variance is larger than the largest possible variance
	#for waiting-time estimator, return that instead.	
	if(var.est > 2*p^2*(1-p)^2 | is.nan(var.est)){
		var.est <- 2*p^2*(1-p)^2
	}
	#if(var.est > 1/12){var.est <- 1/12} #Don't exceed the variance of a uniform(0,1) RV.
	var.est
}


#This estimated variance can be very large
#var.grid <- numeric(191*191*99)
#l <- 1
#for(i in 20:200){
#	for(j in 20:200){
#		for(k in 1:99){
#			var.grid[l] <- var_phat_lineages(i, j, k/100, 1-k/100, .001)
#			l <- l + 1
#		}	
#	}
#}


#variance/N^2 for estimate of N^2 based on waiting-time estimator
var.mult <- function(m, l){
	k <- (m+1):(m+l)
	sum(1/((k^2)*(k-1)^2))/(1/m - 1/(m+l))^2
}


#Takes a vector of the number of lineages ancestral to the sample
#(lins) at each of a vector of times. Estimates a factor proportional
#to the population size between each pair of timepoints in the times vector.
mom.smoothtime <- function(lins, time){
	n0 <- lins[1]
	lins.g0 <- lins[lins>0]
	time.g0 <- time[lins>0]
	#This if statement effectively rounds up if all coalescences happen before
	#the first timepoint 
	if(length(lins.g0) == 1 & lins.g0[1] > 1){
		lins.g0[2] <- 1
		time.g0 <- time[1:2]
	}	
	coal.times <- sapply(lins.g0, mom.scaledtime, n0 = n0)
	coal.smooth <- smooth.steps(time.g0, coal.times)
	rise <- coal.smooth[2:length(coal.smooth)] - coal.smooth[1:(length(coal.smooth) - 1)]
	run <- time.g0[2:length(time.g0)] - time.g0[1:(length(time.g0) - 1)]
	c(run/rise, rep(0, sum(lins==0))) #population size is proportional to reciprocal of slope relating generations to coal time.
}


#Takes a matrix of lineage counts (ref in 1st col, alt in second) and a vector of times
#at which the lineages are counted. Times should be in units proportional to generations.
#Returns an estimated trajectory of alternative allele frequency at times.
est_af_traj_mom.smoothtime <- function(lins, times = seq(0.005, 8.005, by = 0.01)){
	traj <- rep(-1, length(lins[,1]))
	traj[1] <- lins[1,2]/(lins[1,2] + lins[1,1])
	Nref <- mom.smoothtime(lins[,1], times)
	Nalt <- mom.smoothtime(lins[,2], times)
	traj[traj == -1] <- Nalt / (Nref + Nalt)
	traj
}

#Takes a matrix of lineage counts (ref in 1st col, alt in second) and a vector of times
#at which the lineages are counted. Times should be in ***generations***.
#Returns approximate variance of estimated trajectory of alternative allele frequency at times.
# variance of uniform(0,1) is 1/12, so if asymptotic var estimate is larger than
#that, we bring it down to that.
est_af_var_mom.smoothtime <- function(lins, times = 10000*seq(0.005, 8.005, by = 0.01)){
	vars.traj <- rep(-1, length(lins[,1]))
	vars.traj[1] <- (lins[1,2]*lins[1,1]/(lins[1,2] + lins[1,1])^3)
	Nref <- mom.smoothtime(lins[,1], times)
	Nalt <- mom.smoothtime(lins[,2], times)
	for(i in which(vars.traj == -1)){
		vars.traj[i] <- var_phat_lineages(lins[i-1,2], lins[i-1,1], Nalt[i], Nref[i], times[i] - times[i-1])			
	}
	vars.traj
}


#Takes a matrix of lineage counts (ref in 1st col, alt in second) and a vector of times
#at which the lineages are counted. Times should be in ***generations***.
#Returns approximate variance of estimated trajectory of alternative allele frequency at times.
# variance of uniform(0,1) is 1/12, so if asymptotic var estimate is larger than
#that, we bring it down to that.
#This version uses p estimate and coalescent time.
est_af_var_mom.smoothtime_coaltimes <- function(lins, times = seq(0.005, 8.005, by = 0.01)){
	vars.traj <- rep(-1, length(lins[,1]))
	vars.traj[1] <- (lins[1,2]*lins[1,1]/(lins[1,2] + lins[1,1])^3)
	Nref <- mom.smoothtime(lins[,1], times)
	Nalt <- mom.smoothtime(lins[,2], times)
	Nsum <- Nref + Nalt
	for(i in which(vars.traj == -1)){
		vars.traj[i] <- var_phat_lineages(lins[i-1,2], lins[i-1,1], Nalt[i]/Nsum[i], Nref[i]/Nsum[i], times[i] - times[i-1])			
	}
	vars.traj
}




#Take in vector of coalescence times and a value of ell 
#(number of coalescences to wait for before making
#an estimate). Return a matrix with column 1 equal to the 
#max time at which N estimate applies, column 2 equal to N estimate,
#and column 3 variance of the N estimate. If there are less than 
#ell coalescences left after last estimate, just use the ones that
#are left.
estN_waittimes <- function(ctimevec, ell){
	ctimevec <- sort(ctimevec)
	if(length(ctimevec) < ell){inds <- length(ctimevec)}
	if(length(ctimevec) == ell){inds <- ell}	
	if(length(ctimevec) > ell){
		inds <- seq(ell, ell*floor(length(ctimevec)/ell), by = ell)
		if(length(ctimevec)%%ell != 0){
			inds <- c(inds, length(ctimevec))
		}
	}
	ctimes <- ctimevec[inds]
	N.ests <- numeric(length(ctimes))
	N.vars <- N.ests
	wt <- ctimes[1]
	n <- length(ctimevec) + 1
	l <- inds[1]
	N.ests[1] <- wt/(2*(1/(n-l) - 1/n))
	N.vars[1] <- (N.ests[1]^2)*var.mult(n-l, l)
	if(length(ctimes) > 1){
		for(i in 2:length(ctimes)){
			wt <- ctimes[i] - ctimes[i-1]
			n <- length(ctimevec) + 1 - ell*(i-1)
			l <- inds[i] - inds[i-1]		
			N.ests[i] <- wt/(2*(1/(n-l) - 1/n))
			N.vars[i] <- (N.ests[i]^2)*var.mult(n-l, l)
		}
	}	
	cbind(ctimes, N.ests, N.vars)
}



#Given a matrix of the form produced by estN_waittimes() and a target time,
#pulls out the Nestimate and estimated variance of N for the target time.
getN_estNmat <- function(estNmat, targ.time, est.only = FALSE){
	ind <- nrow(estNmat)
	if(estNmat[ind, 1] > targ.time){
		ind <- which(estNmat[,1] == min(estNmat[estNmat[,1] > targ.time,1]))
	}
	if(est.only == TRUE){return(as.numeric(estNmat[ind,2]))}	
	as.numeric(estNmat[ind,2:3])
}


#compute the first-order Taylor-series approximation of the variance
#of a quotient of two random variables.
#mu_n is the expectation of the numerator; var_n is the variance of the numerator;
#mu_d is the expectation of the denominator; var_d is the variance of the denominator;
#cov_nd is the covariance of the numerator and denominator.
#returns 0 if mu_n or var_n is 0.
ts_var_quotient <- function(mu_n, var_n, mu_d, var_d, cov_nd){
	if(mu_n == 0 | var_n == 0){return(0)}	
	((mu_n^2)/(mu_d^2))*(var_n/mu_n^2 - 2*cov_nd/(mu_n*mu_d) + var_d/mu_d^2)
}



#Take in a list of three vectors of coalescence times,
#one for the whole tree (element [[1]]), 
#one for the "ref" subtree (element [[2]]), 
#and one for the "alt" subtree (element [[3]]).
#Also take in a list of times (in the same units as the vector of 
#coalescent times), finally a parameter ell for the ref tree
#and the alt tree. (ell controls how many coalescent events we wait
#for before making an estimate).
#returns estimates of alt allele frequency and estimated variance of 
#frequency estimate. 
#This version assumes that the alt allele is derived and assigns alt
#frequency to 0 before place proportion on the branch on which the mutation
#must have occurred. 
p_ests_wait <- function(ctime.list, time.eval, ell.ref = 5, ell.alt = 5, place = 0.5, ord2adj = FALSE){
	Ns_ref <- estN_waittimes(ctime.list[[2]], ell.ref)
	Ns_alt <- estN_waittimes(ctime.list[[3]], ell.alt)
	#tree join time is the time in full tree that doesn't appear in either subtree.
	tj <- tree_join_time(ctime.list[[1]], ctime.list[[2]][-length(ctime.list[[2]])], ctime.list[[3]][-length(ctime.list[[3]])], NULL)
	lca <- Ns_alt[nrow(Ns_alt), 1]	
	if(tj >= lca){
		Ns_alt <- rbind(Ns_alt, c(lca + (tj - lca)*place, 0, 0))
	}
	if(tj < lca){
		Ns_alt <- rbind(Ns_alt, c(lca + 0.001, 0, 0))
	}
	p.ests <- numeric(length(time.eval))
	var.ests <- p.ests
	for(i in 1:length(time.eval)){
		Ns.r.t <- getN_estNmat(Ns_ref, time.eval[i])
		Ns.a.t <- getN_estNmat(Ns_alt, time.eval[i])
		p.ests[i] <- Ns.a.t[1]/(Ns.a.t[1] + Ns.r.t[1])
		if(ord2adj == TRUE){
			p.ests[i] <- p.ests[i] + Ns.a.t[2]/(Ns.a.t[1] + Ns.r.t[1])^2 - Ns.a.t[1]*(Ns.a.t[2] + Ns.r.t[2])/(Ns.a.t[1] + Ns.r.t[1])^3
		}
		var.ests[i] <- ts_var_quotient(Ns.a.t[1], Ns.a.t[2], Ns.a.t[1] + Ns.r.t[1], Ns.a.t[2] + Ns.r.t[2], Ns.a.t[2])
	}
	cbind(p.ests, var.ests)
}

#Function to compute harmonic mean.
harm.mean <- function(x,...){
	1/mean(1/x,...)
}

#Return the actual harmonic mean Ns for a type between coalescences.
mean_waittimes <- function(pt.time, af.traj, N, ctimevec, ell, mean.fn = "harmonic"){
	ctimevec <- sort(ctimevec)
	if(length(ctimevec) < ell){inds <- length(ctimevec)}
	if(length(ctimevec) == ell){inds <- ell}	
	if(length(ctimevec) > ell){
		inds <- seq(ell, ell*floor(length(ctimevec)/ell), by = ell)
		if(length(ctimevec)%%ell != 0){
			inds <- c(inds, length(ctimevec))
		}
	}
	ctimes <- ctimevec[inds]
	Ns <- numeric(length(ctimes))
	wt <- ctimes[1]
	Ns[1] <- harm.mean(N*af.traj[pt.time < wt])	
	if(length(ctimes) > 1){
		for(i in 2:length(ctimes)){
			if(mean.fn == "harmonic"){
				Ns[i] <- harm.mean(N*af.traj[pt.time < ctimes[i] & pt.time >= ctimes[i-1]])			
			}
			if(mean.fn == "arithmetic"){
				Ns[i] <- mean(N*af.traj[pt.time < ctimes[i] & pt.time >= ctimes[i-1]])			
			}
			
		}
	}	
	cbind(ctimes, Ns, rep(0, length(ctimes)))
}

#ctimevec <- times.c[[200]][[3]]
#af.traj <- mat.trajs[,200]
#harmN_waittimes(pt.time, af.traj, N, ctimevec, 5)

#Function that uses the actual harmonic mean Ns in intervals between wait times.
#Goal is to see whether waiting-time estimator would be biased if estimation
#of harmonic-mean Ns of each type was perfect.
#af.traj is true allele freq trajectory measured at times pt.time.
#N is the population size (assumed to be constant)
p_ests_wait_cheat.mean <- function(pt.time, af.traj, N, ctime.list, time.eval, ell.ref = 5, ell.alt = 5, place = 0.5, mean.fn = "harmonic"){
	Ns_ref <- mean_waittimes(pt.time, 1-af.traj, N, ctime.list[[2]], ell.ref, mean.fn)
	Ns_alt <- mean_waittimes(pt.time, af.traj, N, ctime.list[[3]], ell.alt, mean.fn)
	tj <- tree_join_time(ctime.list[[1]], ctime.list[[2]][-length(ctime.list[[2]])], ctime.list[[3]][-length(ctime.list[[3]])], NULL)
	lca <- Ns_alt[nrow(Ns_alt), 1]	
	if(tj >= lca){
		Ns_alt <- rbind(Ns_alt, c(lca + (tj - lca)*place, 0, 0))
	}
	if(tj < lca){
		Ns_alt <- rbind(Ns_alt, c(lca + 0.001, 0, 0))
	}
	p.ests <- numeric(length(time.eval))
	for(i in 1:length(time.eval)){
		Ns.r.t <- getN_estNmat(Ns_ref, time.eval[i])
		Ns.a.t <- getN_estNmat(Ns_alt, time.eval[i])
		p.ests[i] <- Ns.a.t[1]/(Ns.a.t[1] + Ns.r.t[1])
	}
	p.ests
}


#p.wt <- p_ests_wait(times.c[[50]], time, ell.ref = 1, ell.alt = 1)
#p.hm <- p_ests_wait_cheat.mean(pt.time, af.traj, N, times.c[[50]], time, ell.ref = 1, ell.alt = 1)
#p.am <- p_ests_wait_cheat.mean(pt.time, af.traj, N, times.c[[50]], time, ell.ref = 1, ell.alt = 1, mean.fn = "arithmetic")
#plot(pt.time, af.traj, type = "l")
#lines(time, p.wt[,1], col = "red")
#lines(time, p.am, col = "purple")
#lines(time, p.hm, col = "green")

##Need to refactor waiting-time estimator functions.
#Goal is to allow search across different partitions of the # of coalescent
#events to maximize AIC.

#Take in vector of coalescence times and a vector of ell  values
#(number of coalescences to wait for before making
#an estimate). The entries in the l vector should be a partition
#of the length of ctimevec (i.e. all positive integers, with 
#sum(l.vec) == length(ctimevec)). Return a matrix with column 1 equal to the 
#max time at which N estimate applies, column 2 equal to N estimate,
#and column 3 variance of the N estimate. 
estN_waittimes_partition <- function(ctimevec, l.vec){
	ctimevec <- sort(ctimevec)
	if(sum(l.vec) != length(ctimevec) | mean(l.vec >= 1) != 1 |  !(all.equal(l.vec, as.integer(l.vec)))){
		stop("l.vec must be an ordered partition of the number of coalescences.")
	}
	inds <- cumsum(l.vec)
	ctimes <- ctimevec[inds]
	N.ests <- numeric(length(ctimes))
	N.vars <- N.ests
	wt <- ctimes[1]
	n <- length(ctimevec) + 1
	l <- inds[1]
	N.ests[1] <- wt/(2*(1/(n-l) - 1/n))
	N.vars[1] <- (N.ests[1]^2)*var.mult(n-l, l)
	if(length(ctimes) > 1){
		for(i in 2:length(ctimes)){
			wt <- ctimes[i] - ctimes[i-1]
			n <- length(ctimevec) + 1 - inds[i-1]
			l <- l.vec[i]		
			N.ests[i] <- wt/(2*(1/(n-l) - 1/n))
			N.vars[i] <- (N.ests[i]^2)*var.mult(n-l, l)
		}
	}	
	cbind(ctimes, N.ests, N.vars)
}


#Compute the log likelihood associated with a set of N estimates
#in an matrix like the ones produced by estN_waittimes().
#Requires all the coalescent times (in ctimevec).
loglike.Nests.wait <- function(ctimevec, Nestmat){
	inds.vec <- which(ctimevec %in% Nestmat[,1]) #the events at which N is assessed.
	l.vec <- diff(c(0, inds.vec))
	if(length(l.vec) == length(ctimevec)){
		N.ests <- Nestmat[,2]
	}
	if(length(l.vec) < length(ctimevec)){
		N.ests <- rep(Nestmat[,2], times = l.vec)
	}	
	lin.choice <- choose((length(ctimevec)+1):2, 2)
	lls <- log(lin.choice) - log(N.ests) - c(ctimevec[1], diff(ctimevec)) * lin.choice / N.ests
	sum(lls)
}


#takes in a vector and an index, ind. Returns a vector
#1 shorter than vec, where the indth entry is equal to the sum
#of the indth and ind+1th entries of vec, and all other entries
#are equal to entries in vec.
lumpvec <- function(vec, ind){
	if(ind >= length(vec)){
		stop("The index must be less than the length of the vector.")	
	}
	sum.ent <- vec[ind] + vec[ind + 1]
	if(ind == 1){return(c(sum.ent, vec[-(1:2)]))}
	if(ind == (length(vec) - 1)){return(c(vec[-(ind:(ind+1))], sum.ent ) )}
	c(vec[1:(ind-1)], sum.ent, vec[(ind+2):length(vec)])
}


#Take a vector that represents an ordered partition of an integer.
#Return a list of every vector that represents an ordered partition
#of the same integer and is formed by "lumping" two entries in the 
#original vector.
lumps.lvec <- function(lvec){
	outlist <- list()
	for(i in 1:(length(lvec) - 1)){
		outlist[[i]] <- lumpvec(lvec, i)	
	}
	outlist
}


#Do a greedy search for smoothing with best AIC, where smoothing is determined
#by an ordered partition of #s of coalescences to wait for.
#Start with every interval having its own N estimate, then lump
#neighboring intervals in the most favorable way (AIC-wise) until
#AIC no longer improves or there is just one lump 
#weight.K = 1 gives AIC penalty. weight.K = log(length(ctimevec))/2 gives BIC penalty
find.lvec.AIC.lump <- function(ctimevec, weight.K = 1){
	ctimevec <- sort(ctimevec)
	ctimevec <- breakties.ctimes(ctimevec)	
	S <- length(ctimevec)
	if(S == 1){return(1)}
	K <- S
	best.lvec <- rep(1, length(ctimevec))
	Nmat <- estN_waittimes_partition(ctimevec, best.lvec)
	ll <- loglike.Nests.wait(ctimevec, Nmat)
	best.AIC <- ll - K 
	best.contender <- best.AIC + 1
	while(best.AIC <= best.contender & K > 1){
		contender.AICs <- numeric()
		contenders <- lumps.lvec(best.lvec)
		K <- length(contenders[[1]])		
		for(i in 1:length(contenders)){
			Nmat <- estN_waittimes_partition(ctimevec, contenders[[i]])
			contender.AICs[i] <- loglike.Nests.wait(ctimevec, Nmat) - K * weight.K
		}
		best.contender <- max(contender.AICs)
		if(best.contender > best.AIC){
			best.AIC <- best.contender
			best.lvec <- contenders[[which(contender.AICs == best.contender)[1]]]
		}
	}
	best.lvec
}

#Given an integer to partition (int.to.part), and an integer l,
#returns (as a vector) an ordered partition of int.to.part 
#where all but the final entry
#are equal to l.
get.fixedl.part <- function(int.to.part, l = 1){
	n.ls <- int.to.part %/% l
	rem <- int.to.part %% l
	part <- rep(l, n.ls)	
	if(rem >= 1){part <- c(part, rem)}
	part
}


#Do a search for smoothing parameter l (the number of coalescent events to wait for) 
#with best AIC. If there are fewer than l remaining coalescent events,
#we just wait for all of them to occur.
#weight.K = 1 gives AIC penalty. weight.K = log(length(ctimevec))/2 gives BIC penalty
find.l.AIC <- function(ctimevec, weight.K = 1){
	ctimevec <- sort(ctimevec)
	ctimevec <- breakties.ctimes(ctimevec)	
	S <- length(ctimevec)
	if(S == 1){return(1)}
	parts <- list()
	AICs <- numeric(length(ctimevec))
	for(i in 1:length(ctimevec)){
		part <- get.fixedl.part(length(ctimevec), i)
		Nmat <- estN_waittimes_partition(ctimevec, part)
		AICs[i] <- loglike.Nests.wait(ctimevec, Nmat) - nrow(Nmat) * weight.K
	}
	which.max(AICs)[1]
}



#breaks ties in a vector of coalescent times by adding a small amount to the latter of any
#two times that agree.
breakties.ctimes <- function(ctimevec){
	cdiffs <- diff(ctimevec)
	mindiff <- min(cdiffs[cdiffs > 10^(-10)])
	if(ctimevec[1] < 10^(-10)){
		ctimevec[1] <- mindiff/5
		cdiffs <- diff(ctimevec)
	}
	if(sum(cdiffs < 10^(-10)) == 0){return(ctimevec)}
	tie.inds <- which(cdiffs < 10^(-10)) + 1
	ctimevec[tie.inds] <- ctimevec[tie.inds] + mindiff/5
	if(sum(diff(ctimevec) < 10^(-10)) > 0){return(breakties.ctimes(ctimevec))}
	ctimevec
}


#Take in a list of three vectors of coalescence times,
#one for the whole tree (element [[1]]), 
#one for the "ref" subtree (element [[2]]), 
#and one for the "alt" subtree (element [[3]]).
#Also take in a list of times (in the same units as the vector of 
#coalescent times).
#returns estimates of alt allele frequency and estimated variance of 
#frequency estimate. 
#This version assumes that the alt allele is derived and assigns alt
#frequency to 0 before place proportion on the branch on which the mutation
#must have occurred. 
p_ests_wait_AICpartition <- function(ctime.list, time.eval, place = 0.5, weight.K.ref = 1, weight.K.alt = 1){
	lvec.ref <- find.lvec.AIC.lump(ctime.list[[2]], weight.K.ref)
	lvec.alt <- find.lvec.AIC.lump(ctime.list[[3]], weight.K.alt)
	Ns_ref <- estN_waittimes_partition(breakties.ctimes(sort(ctime.list[[2]])), lvec.ref)
	Ns_alt <- estN_waittimes_partition(breakties.ctimes(sort(ctime.list[[3]])), lvec.alt)
	#tree join time is the time in full tree that doesn't appear in either subtree.
	tj <- tree_join_time(ctime.list[[1]], ctime.list[[2]][-length(ctime.list[[2]])], ctime.list[[3]][-length(ctime.list[[3]])], NULL)
	lca <- Ns_alt[nrow(Ns_alt), 1]	
	if(tj >= lca){
		Ns_alt <- rbind(Ns_alt, c(lca + (tj - lca)*place, 0, 0))
	}
	if(tj < lca){
		Ns_alt <- rbind(Ns_alt, c(lca + 0.001, 0, 0))
	}
	p.ests <- numeric(length(time.eval))
	var.ests <- p.ests
	for(i in 1:length(time.eval)){
		Ns.r.t <- getN_estNmat(Ns_ref, time.eval[i])
		Ns.a.t <- getN_estNmat(Ns_alt, time.eval[i])
		p.ests[i] <- Ns.a.t[1]/(Ns.a.t[1] + Ns.r.t[1])
		var.ests[i] <- ts_var_quotient(Ns.a.t[1], Ns.a.t[2], Ns.a.t[1] + Ns.r.t[1], Ns.a.t[2] + Ns.r.t[2], Ns.a.t[2])
	}
	cbind(p.ests, var.ests)
}



#Take in a list of three vectors of coalescence times,
#one for the whole tree (element [[1]]), 
#one for the "ref" subtree (element [[2]]), 
#and one for the "alt" subtree (element [[3]]).
#Also take in a list of times (in the same units as the vector of 
#coalescent times).
#returns estimates of alt allele frequency and estimated variance of 
#frequency estimate. 
#This version assumes that the alt allele is derived and assigns alt
#frequency to 0 before place proportion on the branch on which the mutation
#must have occurred. 
p_ests_wait_AIC_l <- function(ctime.list, time.eval, place = 0.5, weight.K.ref = 1, weight.K.alt = 1){
	l.ref <- find.l.AIC(ctime.list[[2]], weight.K.ref)
	l.alt <- find.l.AIC(ctime.list[[3]], weight.K.alt)
	lvec.ref <- get.fixedl.part(length(ctime.list[[2]]), l.ref)
	lvec.alt <- get.fixedl.part(length(ctime.list[[3]]), l.alt)
	Ns_ref <- estN_waittimes_partition(breakties.ctimes(sort(ctime.list[[2]])), lvec.ref)
	Ns_alt <- estN_waittimes_partition(breakties.ctimes(sort(ctime.list[[3]])), lvec.alt)
	#tree join time is the time in full tree that doesn't appear in either subtree.
	tj <- tree_join_time(ctime.list[[1]], ctime.list[[2]][-length(ctime.list[[2]])], ctime.list[[3]][-length(ctime.list[[3]])], NULL)
	lca <- Ns_alt[nrow(Ns_alt), 1]	
	if(tj >= lca){
		Ns_alt <- rbind(Ns_alt, c(lca + (tj - lca)*place, 0, 0))
	}
	if(tj < lca){
		Ns_alt <- rbind(Ns_alt, c(lca + 0.001, 0, 0))
	}
	p.ests <- numeric(length(time.eval))
	var.ests <- p.ests
	for(i in 1:length(time.eval)){
		Ns.r.t <- getN_estNmat(Ns_ref, time.eval[i])
		Ns.a.t <- getN_estNmat(Ns_alt, time.eval[i])
		p.ests[i] <- Ns.a.t[1]/(Ns.a.t[1] + Ns.r.t[1])
		var.ests[i] <- ts_var_quotient(Ns.a.t[1], Ns.a.t[2], Ns.a.t[1] + Ns.r.t[1], Ns.a.t[2] + Ns.r.t[2], Ns.a.t[2])
	}
	cbind(p.ests, var.ests)
}







#Function that returns a list of three vectors of coalescence times,
#one for the whole tree (element [[1]] of the output), 
#one for the "ref" subtree (element [[2]] of the output), 
#and one for the "alt" subtree (element [[3]] of the output).
#non-null values of mut.time override the tree_join_time() function. Can 
#be used with place = 1 to dictate the specific time at which mutation occurred.
#units_in and units_out can be used to adjust times. For example, if trees are from
#ms, which uses units of 4N, and units of 2N are desired for output times, then
#set units_in = 4 and units_out = 2.
trees_to_times <- function(tree.all, tree.ref, tree.alt, times, sure.alt.is.derived = FALSE, place = 0.5, mut.time = NULL, units_in = 2, units_out = 2){
	allt <- coal.times(tree.all)
	reft <- coal.times(tree.ref)
	altt <- coal.times(tree.alt)
	if(units_in != units_out){
		allt <- allt*units_in/units_out
		reft <- reft*units_in/units_out
		altt <- altt*units_in/units_out
	}
	tj <- tree_join_time(allt, reft, altt, mut.time)
	add_mutation(allt, reft, altt, tj, times, sure.alt.is.derived, place)
}


#Takes a list of three vectors of coalescence times,
#one for the whole tree (element [[1]] of the output), 
#one for the "ref" subtree (element [[2]] of the output), 
#and one for the "alt" subtree (element [[3]] of the output).
#and a vector of times at which to compute estimator.
#at time 0 (or less), nothing has coalesced, so add one to
#the length of each vector to get the number of tips per tree.
#Returns a matrix with the numbers of lineages of ref allele 
#(column 1) and alt allele (column 2) at each time in times.
#(each row of output is for the corresponding entry in times.)
times_to_lins <- function(tree.times, times = seq(0.005, 8.005, by = 0.01)){
	reft <- tree.times[[2]]
	altt <- tree.times[[3]]
	count.lins <- function(x, vec){
		if(x <= 0){
			return(length(vec) + 1)		
		}
		sum(vec > x)
	}
	lins_ref <- sapply(times, FUN = count.lins, vec = reft )
	lins_alt <- sapply(times, FUN = count.lins, vec = altt )
	cbind(lins_ref, lins_alt)
}


#Takes a matrix of lineage numbers per time and returns
#the proportion at each timepoint in the second (alt) column.
est_af_traj_neut <- function(lins){
	lins[,2]/rowSums(lins)
}


#Takes a matrix of lineage numbers per time and returns
#the binomial sampling variance of the neutral MLE of the allele frequency.
est_af_var_neut_bin <- function(lins){
	j <- lins[,2]
	r <- rowSums(lins)	
	j*(r-j)/(r^3)
}

#Takes a matrix of lineage numbers per time and returns
#the posterior variance of the neutral estimator when viewed in a Bayesian way.
est_af_var_neut_post <- function(lins){
	j <- lins[,2]
	r <- rowSums(lins)	
	j*(r-j)/((r^2)*(r+1))
}




#Take in a tree (as an ape/phylo object) and stretch/shrink 
#its branch lengths randomly, in a way simulating random observation
#of distances due to Poisson mutation.
#theta.region is a theta parameter for a region, giving mutation rate
#along branch in coal units. Larger values mean the tree gets warped less.
#force.ultra forces the tree to be ultrametric.
#tol is a number any distances less than which are replaced with a larger (but still small)
#random number.
#Returns a phylo object with randomized distances.
warp_a_tree <- function(tree, theta.region, force.ultra = TRUE, tol = 1e-10){
	edges.unique <- unique(tree$edge.length)
	observed.edges <- rnorm(length(edges.unique), edges.unique, sqrt(edges.unique / theta.region))
	#observed.edges[observed.edges < tol] <- runif(sum(observed.edges < tol), tol, min(observed.edges[observed.edges > tol]))
	observed.edges[observed.edges < tol] <- edges.unique[observed.edges < tol]
	tree.warped <- tree	
	tree.warped$edge.length <- observed.edges[match(tree$edge.length, edges.unique)]
	if(force.ultra){ #tree is forced to be ultrametric
		dismat <- cophenetic(tree.warped)
		tree.warped <- nnls.tree(dismat,tree.warped,rooted=TRUE,trace=0)
	}
	tree.warped
}






#####Functions for hypothesis testing


#Given two vectors of (possibly estimated) allele frequencies and associated effect sizes
#(assumed to remain constant), estimate selection gradient. 
#No timing is input here, but if known, divide this estimate by the number
#of generations to estimate the per-generation gradient.
#change is calculated forward in time, assuming that af.1 is more recent.
#Here's the theory:
#assume that delta.p is distributed as N(delta.t * 2N * s * p *(1-p), delta.t * p * (1-p))
#where p is allele frequency (at either timepoint, should be equivalent by reversibility?)
#delta.t is in coalescent (1 = 2N gens) units, and 2N is # of chromosomes in pop (i.e. diploids).
#This is a short-time approx to diffusion.
#then we estimate 
#beta.est = sum(eff_sizes * delta.p) / sum(eff_sizes^2 * p *(1-p))
#Conditional on p, we get expectation by plugging in expectation of delta.p, or
#E(beta.est) = sum(eff_sizes * (delta.t *2N * beta * eff_sizes * p *(1-p))) / sum(eff_sizes^2 * p *(1-p))
#            = beta * delta.t * 2N
#That is, beta * the number of generations that have passed.
#To make an estimator of beta itself, we can divide by the number of generations that have passed.
#To get the variance (conditional on p), we use
#Var(beta.est) = (1 / sum(eff_sizes^2 * p * (1 - p)))^2 * sum(eff_sizes^2 * Var(delta.p)) 
#              = delta.t / sum(eff_sizes^2 * p * (1 - p))
#if we divide beta.est by (delta.t * 2N), then the variance is
#Var(beta.est / (delta.t * 2N)) = 
est.sel.grad <- function(eff_sizes, af.1, af.2){
	delta.p <- af.1 - af.2
	sum(eff_sizes * delta.p)/sum(eff_sizes^2 * af.1 * (1 - af.1))
}


#delta.t is in coalescent units.
var.sel.grad <- function(eff_sizes, af.1, delta.t){
	delta.t / sum(eff_sizes^2 * af.1 * (1-af.1))
}



#Returns a permutation distribution of the estimated selection gradient.
#effect sizes are permuted, which entails an assumption
#that there's never been selection on any of these loci, but probably
#most sensitive to selection between the two timepoints.
perm.sel.grad <- function(eff_sizes, af.1, af.2, n.perms = 10000){
	reps <- matrix(rep(eff_sizes, n.perms), ncol = n.perms)
	perms <- apply(reps, 2, sample)
	apply(perms, 2, est.sel.grad, af.1 = af.1, af.2 = af.2)
}



#Returns a permutation distribution of the estimated selection gradient.
#effect sizes are sign-flipped.
perm.sel.grad.signflip <- function(eff_sizes, af.1, af.2, n.perms = 10000){
	reps <- matrix(rep(eff_sizes, n.perms), ncol = n.perms)
	sign.flip <- function(vec){
		signs <- rbinom(length(vec), 1, 0.5)*2 - 1
		signs * vec
	}	
	perms <- apply(reps, 2, sign.flip)
	apply(perms, 2, est.sel.grad, af.1 = af.1, af.2 = af.2)
}




#est.sel.grad(eff_sizes, trajs_neut[50,], trajs_neut[100,])
#psg <- perm.sel.grad(eff_sizes, trajs_neut[50,], trajs_neut[100,])

#fits a weighted ls with no intercept, relating allele-frequency
#change to effect size
weighted.lm.sel.grad.ni <- function(eff_sizes, af.1, af.2){
	af.chgs <- af.1 - af.2
	lm(af.chgs ~ eff_sizes - 1, weights = 1/(af.1*(1-af.1)))
}

#summary(weighted.lm.sel.grad(eff_sizes, trajs_neut[50,], trajs_neut[100,]))
#summary(weighted.lm.sel.grad(eff_sizes, trajs_est_wt_l1[50,], trajs_est_wt_l1[100,]))


#fits a weighted ls with intercept, relating allele-frequency
#change to effect size
weighted.lm.sel.grad <- function(eff_sizes, af.1, af.2){
	af.chgs <- af.1 - af.2
	lm(af.chgs ~ eff_sizes, weights = 1/(af.1*(1-af.1)))
}


#Runs a spearman's rho test relating allele frequency change to effect size.
spearman.seltest <- function(eff_sizes, af.1, af.2){
	af.chgs <- af.1 - af.2
	cor.test(af.chgs, eff_sizes, method = "spearman")
}

#spearman.seltest(eff_sizes, trajs_neut[50,], trajs_neut[100,])

#computes the F matrix for a modified version of Q_x. The present
#population replaces the "ancestral" population (by reversibility
#argument). The times vector contains the times (in coalescent units)
#at which the allele frequencies to be tested are assessed.
#each entry f_i_ and f__i contains the ith time.
#note that the first entry of times should not be 0.
fmat.time <- function(times){
	fmat <- matrix(times[length(times)], nrow = length(times), ncol = length(times))
	for(i in (length(times)-1):1){
		fmat[i,] <- times[i]
		fmat[,i] <- times[i]
	}
	fmat
}


#Estimates coalescent time between two timepoints given two vectors of allele frequencies
#corresponding to the two timepoints. The allele frequency vector af.1 is treated as fixed.
#The estimate is based on the short-time approximation that given af.1, af.2 will be normal with
#expectation af.1 and variance t * af.1 * (1 - af.1), where t is coalescent time.
#Thus, time can be estimated as the sample variance of (af.1 - af.2) / sqrt(af.1 * (1 - af.1)).
#allele frequencies that are 0 or 1 in the first timepoint (af.1) are excluded from the calculation.
est.time <- function(af.1, af.2){
	chgs <- af.1 - af.2
	chgs <- chgs[af.1 > 0 & af.1 < 1]
	scale <- sqrt(af.1*(1 - af.1))
	scale <- scale[af.1 > 0 & af.1 < 1]
	var(chgs/scale)
}

#Takes a matrix with allele frequencies for distinct loci in different columns
#and allele frequencies for distinct time points in different rows.
#Estimates the coalescent time passed between time points (rows) using
#est.time(). Each time estimate is independent of all the others.
est.times.seq <- function(trajs){
	n.times <- nrow(trajs) - 1
	times <- numeric(n.times)	
	for(i in 1:n.times){
		times[i] <- est.time(trajs[i,], trajs[i+1,] )	
	}
	times
}

#Takes a matrix with allele frequencies for distinct loci in different columns
#and allele frequencies for distinct time points in different rows.
#Estimates the coalescent time passed between time points (rows) using
#est.time(). Rather than estimating times between adjacent measurements,
#all times are estimated with respect to the first timepoint.
est.times.fromstart <- function(trajs){
	n.times <- nrow(trajs) - 1
	times <- numeric(n.times)	
	for(i in 1:n.times){
		times[i] <- est.time(trajs[1,], trajs[i+1,] )	
	}
	c(times[1], diff(times)) #express as times between points, like est.times.seq()
}


#Compute additive genetic variance at the present.
add.var <- function(eff_sizes, freqs, ploidy = 2){
	ploidy*sum((eff_sizes^2) * freqs * (1 - freqs))
}

#add.var(eff_sizes, curr.freqs)


#Computes modified Qx for one trajectory of estimated phenotypes.
#traj is a trajectory of estimated avg trait values going back into 
#the past (not including present), *with the value at the present subtracted off*.
#times are the times (in coal units) at which the mean phenotype is estimated.
#You can pass in an fmat that's already been inverted to save a little time.
#if the fmat is inverted already, set fmat.inv = TRUE; otherwise use the default.
compute_Qx_traj <- function(traj, fmat, Va, f.inverted = FALSE){
	fma <- fmat	
	if(f.inverted == FALSE){
		fma <- solve(fmat)
	}	
	(t(traj) %*% fma %*% traj) / (2*Va)
}

#compute_Qx_traj(traj.phen.wt_l1[-1] - traj.phen.neut[1], fm, 1)
#fm.inv <- solve(fm)
#compute_Qx_traj(traj.phen.wt_l1[-1] - true.per.time[1], fm.inv, 1, f.inverted = TRUE)


#function to compute a Qx test given a matrix of allele frequency trajectories
#and a vector of effect sizes. If timevec = NULL, the times between observations
#are estimated according to allele-frequency variation. Otherwise they must be
#supplied. If perms == 0, then p values are computed according to chisq distribution;
#otherwise effect sizes are permuted to get a null distribution. 
#It's assumed that each row of trajmat
#has allele freqs (in the same order as the eff_sizes vector) at some point in
#time, with the most recent time in the first row and subsequent rows
#in time order into the past. The most recent time is treated as the "ancestral"
#measurement and subtracted from everything else. (Va is also computed
#on the basis of the most recent observation). trajmat must have
#at least two rows.
#If time is estimated, it can be estimated either from the start (timeest = "seq") 
#or in sequence ("seq", the default),
#meaning that each interval is estimated independently of the others.
Qx_test <- function(trajmat, eff_sizes, timevec = NULL, perms = 1000, timeest = "seq"){
	if(is.null(timevec)){
		if(timeest == "seq"){		
			timevec <- cumsum(est.times.seq(trajmat))
		}
		if(timeest == "fromstart"){
			timevec <- cumsum(est.times.fromstart(trajmat))
		}			
	}
	fm <- fmat.time(timevec)
	fma <- solve(fm)
	Va <- add.var(eff_sizes, trajmat[1,])
	traj.phen <- as.numeric(2 * eff_sizes %*% t(trajmat))
	traj.phen.adj <- traj.phen[-1] - traj.phen[1]	
	Q_stat <- as.numeric(compute_Qx_traj(traj.phen.adj, fma, Va, f.inverted = TRUE))
	if(perms == 0){
		return(c(Q_stat, length(timevec), pchisq(Q_stat, length(timevec), lower.tail = FALSE)))	
	}
	permdist <- numeric(perms)
	for(i in 1:perms){
		perm.effs <- sample(eff_sizes)
		Va.p <- add.var(perm.effs, trajmat[1,])
		traj.phen <- as.numeric(2 * perm.effs %*% t(trajmat))
		traj.phen.adj <- traj.phen[-1] - traj.phen[1]
		permdist[i] <- as.numeric(compute_Qx_traj(traj.phen.adj, fma, Va.p, f.inverted = TRUE))
	}
	c(Q_stat, length(timevec), mean(permdist > Q_stat))
}

#Qx_test(trajs_neut[time %in% ((0:10)/100),], eff_sizes, perms = 10000)




#Using Fu and Li, Genetics (1993), Eq. 14
#computes the variance of the sum of terminal branch lengths
#for a tree of n leaves from a population of constant size 2N.
#variance is for the sum terminal branch length in
#coalescent units of 2N generations.
#ntot is an optional argument for the total number of 
var.term.mean <- function(n){
	if(n == 1){
		return(NA)	
	}	
	if(n == 2){c_n <- 1}
	if(n > 2){
		a_n <- sum(1/(1:(n-1)))
		c_n <- 2*(n * a_n - 2*(n-1))/((n-1)*(n-2))  	
	}
	varn <- 4 * c_n / n^2
}


#Computes an analogue of tSDS from trees. In particular, computes difference in
#mean singleton branch length in the ref trees and alt trees, scales them
#according to an approximate variance. These are signed so that they're positive
#if the alt allele has shorter terminal branch lengths (i.e. the alt allele has been selected up)
#Uses Fu and Li variance, but multiplies by p^2 (est by allele freq in sample) to 
#keep in in units of 2N gens instead of 2Np gens.
tSDS_analogue <- function(ref_trees, alt_trees){
	extract.n <- function(x){length(x$tip.label)}
	extract.tbl.sum <- function(x){sing.bl(branch.length.by.n.descendants(x))}
	nref <- sapply(ref_trees, FUN = extract.n)
	nalt <- sapply(alt_trees, FUN = extract.n)
	termsum_ref <- sapply(ref_trees, FUN = extract.tbl.sum)
	termsum_alt <- sapply(alt_trees, FUN = extract.tbl.sum)
	tm_ref <- termsum_ref / nref
	tm_alt <- termsum_alt / nalt
	var.tmref <- sapply(nref, var.term.mean)*(nref/(nref+nalt))^2
	var.tmalt <- sapply(nalt, var.term.mean)*(nalt/(nref+nalt))^2
	n <- (nref + nalt)[1]	
	var.1coal <- 4 * (2*(sum(1/(1:n)) - 1)/(n*(n-1)) - 1/(n^2)) #from Fu & Li 1993, var coal time for one branch.
	var.tmref[is.na(var.tmref)] <- var.1coal
	var.tmalt[is.na(var.tmalt)] <- var.1coal
	#var.tmalt[nref < 5 | nalt < 5] <- NA
	#tm_ref[is.na(var.tmref)] <- tm_ref[is.na(var.tmref)] / 2
	(tm_ref - tm_alt)/sqrt(var.tmref + var.tmalt)
}


#ts <- tSDS_analogue(anc_trees_ms, der_trees_ms)

#sum(ts * sign(eff_sizes))

#ts.perm <- numeric(10000)
#for(rep in 1:10000){
#	ts.perm[rep] <- sum(ts * sign(sample(eff_sizes)))
#}


#like tSDS, but normalize by means and sds from neutral sims, which are provided by the user.
tSDS_analogue_meansd <- function(ref_trees, alt_trees, means.daf, sds.daf){
	extract.n <- function(x){length(x$tip.label)}
	extract.tbl.sum <- function(x){sing.bl(branch.length.by.n.descendants(x))}
	nref <- sapply(ref_trees, FUN = extract.n)
	nalt <- sapply(alt_trees, FUN = extract.n)
	termsum_ref <- sapply(ref_trees, FUN = extract.tbl.sum)
	termsum_alt <- sapply(alt_trees, FUN = extract.tbl.sum)
	tm_ref <- termsum_ref / nref
	tm_alt <- termsum_alt / nalt
	(tm_ref - tm_alt - means.daf) / sds.daf
}







##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
#No-longer-used and unused functions.



#Compute probability that 
#we shift from m0 to mt lineages of type 1
#we shift from n0 to nt lineages of type 2
#AND
#the frequency of type 1 changes from p0 to pt, all in t generations.
#assumes fixed haploid population size N and constant selection coefficient
#s (forward in time; this computes probability of a trajectory backward, and so
#the negative of the forward-in-time s is used) on alleles of type 1. 
#form, frac.swit, and func.user are as in time.scaling().
logposterior.pt <- function(m0, mt, n0, nt, p0, pt, N, t, s, form = "linear", frac.swit = 0.5, func.user = NULL){
	log_g_m <- p.ancs(m0, mt, N*p0, N*pt, t, form, frac.swit, func.user, logp = TRUE)
	log_g_n <- p.ancs(n0, nt, N*(1-p0), N*(1-pt), t, form, frac.swit, func.user, logp = TRUE)
	mu <- p0 - s*p0*(1-p0)*t
	sig2 <- p0*(1-p0)*t/N
	log_p_shift <- dnorm(pt, mu, sqrt(sig2), log = TRUE)
	log_g_m + log_g_n + log_p_shift
}

#density of the Laplace distribution.
#code from the function of the same name in package rmutil.
#just using this directly so that we don't have to load in the package.
dlaplace <- function (y, m = 0, s = 1, log = FALSE) 
{
    if (any(s <= 0)) 
        stop("s must be positive")
    tmp <- -abs(y - m)/s - log(2 * s)
    if (!log) 
        tmp <- exp(tmp)
    tmp
}

#like logposterior.pt(), but puts a Laplace distribution on s, centered at 0
#and with a user-supplied dispersion.
logposterior.pt.s <- function(m0, mt, n0, nt, p0, pt, N, t, s, s.disp = 0.01/sqrt(2), form = "linear", frac.swit = 0.5, func.user = NULL){
	log_p_s <- dlaplace(s, 0, s.disp, log = TRUE)	
	log_g_m <- p.ancs(m0, mt, N*p0, N*pt, t, form, frac.swit, func.user, logp = TRUE)
	log_g_n <- p.ancs(n0, nt, N*(1-p0), N*(1-pt), t, form, frac.swit, func.user, logp = TRUE)
	mu <- p0 - s*p0*(1-p0)*t
	sig2 <- p0*(1-p0)*t/N
	log_p_shift <- dnorm(pt, mu, sqrt(sig2), log = TRUE)
	log_p_s + log_g_m + log_g_n + log_p_shift
}


#Find the value of pt that maximizes the "posterior" computed in logposterior.pt().
max.posterior.pt <- function(m0, mt, n0, nt, p0, N, t, s, form = "linear", frac.swit = 0.5, func.user = NULL){
	post <- function(pt, m0, mt, n0, nt, p0, N, t, s, form = "linear", frac.swit = 0.5, func.user = NULL){
		-logposterior.pt(m0, mt, n0, nt, p0, pt, N, t, s, form, frac.swit, func.user)
	}
	optim(p0, post, method = "Brent", lower = 0, upper = 1, m0 = m0, mt = mt, n0 = n0,
		nt = nt, p0 = p0, N = N, t = t, s = s, form = form, frac.swit = frac.swit,
		func.user = func.user)$par
}




#draw a sample from the posterior distribution of pt using rejection sampling.
#candidates drawn from a uniform--probably could be made more efficient.
#possible to provide the pt that maximizes the posterior and the value
#of the max at that point to save time.
post.rej.samp <- function(m0, mt, n0, nt, p0, N, t, s, max.post = NULL, argmax.post = NULL, form = "linear", frac.swit = 0.5, func.user = NULL){
	if(is.null(argmax.post)){
		argmax.post <- max.posterior.pt(m0, mt, n0, nt, p0, N, t, s, form, frac.swit, func.user)
	}
	if(is.null(max.post)){
		max.post <- exp(logposterior.pt(m0, mt, n0, nt, p0, argmax.post, N, t, s, form, frac.swit, func.user))
	}
	samp <- NULL
	while(is.null(samp)){
		cand <- runif(1, 0, 1)
		cand.post <- exp(logposterior.pt(m0, mt, n0, nt, p0, cand, N, t, s, form, frac.swit, func.user))
		comp <- runif(1,0,1)
		if(cand.post/max.post >= comp){
			samp <- cand
		}
	}
	samp
}




#Estimate allele frequency trajectory by maximizing posterior.
#Entry i in ss is the selection coefficient (advantage of heterozygote)
#between times[i] and times[i+1]. N is the fixed haploid pop size.
est_af_traj_max.posterior <- function(lins, N, times = seq(0.005, 8.005, by = 0.01), p0=NULL, ss = rep(0,length(times-1)), form = "linear", frac.swit = 0.5, func.user = NULL){
	traj <- rep(-1, length(lins[,1]))
	traj[lins[,1] == 0] <- 1
	traj[lins[,2] == 0] <- 0
	traj[1] <- lins[1,2]/(lins[1,2] + lins[1,1])
	if(!is.null(p0)){traj[1] <- p0}
	for(i in which(traj < 0)){
		traj[i] <- max.posterior.pt(lins[i-1,2], lins[i,2], lins[i-1,1], lins[i,1], traj[i-1], N, N*(times[i] - times[i-1]), ss[i-1], form,frac.swit, func.user)
	}
	traj
}


##Computes a trajectory by sampling from the posterior at each timepoint, going back
#into the past.
est_af_traj_samp.posterior <- function(lins, N, times = seq(0.005, 8.005, by = 0.01), p0=NULL, ss = rep(0,length(times-1)), form = "linear", frac.swit = 0.5, func.user = NULL){
	traj <- rep(-1, length(lins[,1]))
	traj[lins[,1] == 0] <- 1
	traj[lins[,2] == 0] <- 0
	traj[1] <- lins[1,2]/(lins[1,2] + lins[1,1])
	if(!is.null(p0)){traj[1] <- p0}
	for(i in which(traj < 0)){
		traj[i] <- post.rej.samp(lins[i-1,2], lins[i,2], lins[i-1,1], lins[i,1], traj[i-1], N, N*(times[i] - times[i-1]), ss[i-1], form=form, frac.swit=frac.swit, func.user=func.user)
	}
	traj
}





#Produce 1-generation transition matrix for Wright-Fisher model.
#N is haploid population size (so use 2*N for diploids)
#s is selection coefficient
wf_transition_1gen <- function(N, s = 0){
	ps <- rep((0:N)/N, N+1)
	qs <- (ps*(1+s))/(1 + ps*s)
	ns <- rep(0:N, each = N+1)
	matrix(dbinom(ns, N, qs), nrow = N+1, ncol = N+1)
}

#Given an eigendecomposition of a 1-generation transition matrix
#(produced by calling eigen()), produces the t-generation
#transition matrix
wf_transition_tgen <- function(eig, t, P.inv = NULL){
	P <- eig$vectors
	D <- diag(eig$values^t)
	if(is.null(P.inv)){P.inv <- solve(P)}
	trans <- P %*% D %*% P.inv
	Re(trans)
}


#Compute the probability that a new mutation is lost before ever reaching/exceeding
#a target frequency.
p_quickloss <- function(N, s = 0, targ = 1/1000){
	if(1/N >= targ){return(0)}
	targ.n <- ceiling(N*targ)
	ps <- rep((0:targ.n)/N, targ.n+1)
	qs <- (ps*(1+s))/(1 + ps*s)
	ns <- rep(0:targ.n, each = targ.n + 1)
	t.1gen <- matrix(dbinom(ns, N, qs), nrow = targ.n+1, ncol = targ.n+1)
	t.1gen[nrow(t.1gen),] <- c(rep(0, targ.n), 1)
	t.1gen[,ncol(t.1gen)] <- 1 - rowSums(t.1gen[,1:targ.n])	
	eig.ob <- eigen(t.1gen)
	wf_transition_tgen(eig.ob, targ.n*10)[2,1]
}



p_samp_g0 <- function(N, t, n, s){
	wf.1g <- wf_transition_1gen(N, s)
	eig <- eigen(wf.1g)
	wf.tg <- wf_transition_tgen(eig, t)
	p_pop <- wf.tg[2,]
	i <- 0:N
	1 - sum(p_pop * (1 - i/N)^n)
}










