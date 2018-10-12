#Doc Edge, 4/13/17
#Goal: Run rent+ on 1700 height-associated SNPs from 1KG.

library(psych)
library(data.table)
library(ape)

id_file <- "GBR_IDS.txt"
vcf_file <- "test.vcf"
out_pref <- "hap"
rent_file <- "rent_in.txt"
rentplus_fn <- "/home/medge/Data/coal_selection_sims/RentPlus.jar"

height_fn <- "/home/medge/Data/coal_selection_sims/Biobank_height_013018/biobank_height_maxBF_PickrellWindows.txt"

helper_fn <- "/home/medge/Data/coal_selection_sims/helper_functions_coal_sel.R"

source(helper_fn) #read in helper functions and load ape package


mark.dat <- fread(height_fn, head = TRUE)
nmarks <- nrow(mark.dat)

wid <- 100000 #half-window size around focal locus to extract

n_ders <- numeric(0)
thetas <- numeric(0)

bls.by.ndesc.all <- list()
bls.by.ndesc.anc <- list()
bls.by.ndesc.der <- list()

trees.all <- list()
trees.anc <- list()
trees.der <- list()


total.bls.all <- numeric(0)
total.bls.anc <- numeric(0)
total.bls.der <- numeric(0)

splits <- strsplit(mark.dat$variant, ":")

chr.dat <- numeric()
pos.dat <- numeric()
anc.allele.dat <- character()
der.allele.dat <- character()
for(i in 1:length(splits)){
	chr.dat[i] <- as.numeric(splits[[i]][1])
	pos.dat[i] <- as.numeric(splits[[i]][2])
	anc.allele.dat[i] <- splits[[i]][3]
	der.allele.dat[i] <- splits[[i]][4]
}


seg.sites <- numeric()
meas.length <- numeric()
n.chrs <- numeric()

monomorphic_loci <- numeric()




for(i in 1:dim(mark.dat)[1]){
   snp <- mark.dat$rsid[i]
   pos.snp <- pos.dat[i]
   chr <- chr.dat[i]
   lb <- pos.snp - wid
   ub <- pos.snp + wid

   #Write a string that will call tabix to extract a vcf from 1000 genomes.
   tabix.str <- paste("tabix -h /home/medge/Data/1000gen_rel_20130502/ALL.chr", as.character(chr), ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ", as.character(chr), ":", as.character(lb), "-", as.character(ub), " | vcf-subset -c ", id_file, " > ", vcf_file, sep = "")
   system(tabix.str)

   #Write to IMPUTE format, which is easy to change to RENT+ format.
   vcftools.str <- paste("vcftools --vcf ", vcf_file, " --out ", out_pref, " --IMPUTE", sep = "")
   system(vcftools.str)

   #Read in the IMPUTE format haplotypes, keep polymorphic loci, and 
   #rearrange into RENT+ format. We shuffle so that the derived alleles
   #are at the end.

   hap.fn <- paste(out_pref, ".impute.hap",sep = "")
   haps <- try(as.matrix(read.table(hap.fn)))
   if (inherits(haps, 'try-error')){
       print(paste("SNP", as.character(i), "has a problem")) 
       next
   }
   info.fn <- paste(out_pref, ".impute.legend",sep = "")
   info <- read.table(info.fn, head = TRUE)
   rsIDs.char <- as.character(info[,1])
   
   loci.poly <- rowMeans(haps) > 0
   haps.poly <- haps[loci.poly,]
   info.poly <- info[loci.poly,]
   col.keysnp <- which(info.poly$pos == as.numeric(pos.snp))
   if(length(col.keysnp) == 0){
      print(paste("SNP", as.character(i), "missing from 1KG"))    
      next 
   }
   
   if(i == 1){
       info.keysites <- info[info$pos == pos.snp,]
       info.keysites[,1] <- rsIDs.char[info$pos == pos.snp]
   }
   info.keysites[i,] <- info[info$pos == pos.snp,]
   info.keysites[i,1] <- rsIDs.char[info$pos == pos.snp]
   #print(i)
   #print(snp)
   #print(pos.snp)
   #print(splits[[i]])
   #print(info[info$pos == pos.snp,])



   n_der <- sum(haps.poly[col.keysnp,])
   haps.poly.der.at.end <- haps.poly[, order(haps.poly[col.keysnp,])]
   haps.charvec <- apply(t(haps.poly.der.at.end), 1, paste, sep = "", collapse = "")
   haps.char <- paste(haps.charvec, sep = "", collapse = "\n")
   pos.char <- paste(as.character(info.poly$pos), sep = "", collapse = " ")
   rentfile.str <- paste(pos.char, haps.char, sep = "\n", collapse = "")
   write(rentfile.str, rent_file)

   seg.sites[i] <- sum(loci.poly)
   n.chrs[i] <- ncol(haps)
   n_ders[i] <- n_der
   
   meas.length[i] <- max(info.poly$pos) - min(info.poly$pos)
   print(paste("segsites =", as.character(seg.sites[i])))
   print(paste("measured length =", as.character(meas.length[i])))

   if(n_der == n.chrs[i] | n_der == 0){
      monomorphic_loci <- c(monomorphic_loci, i)
      print(paste("SNP", as.character(i), "is monomorphic in the sample")) 
      next
   }

   #run rent+
   #rent.string <- paste("java -jar ", rentplus_fn, " -t ",  rent_file, sep = "")
   #system(rent.string)

   #read in theta estimate
   rent_theta_fn <- paste(rent_file, ".theta", sep = "")
   thetas[i] <- suppressWarnings(as.numeric(read.table(rent_theta_fn, nrows = 1)$V1))

   #Read in trees at selected site in rent and rent+
   rent_trees_fn <- paste(rent_file, ".trees", sep = "")
   rent_trees <- readLines(rent_trees_fn)

   #which tree is for the focal site?
   rent_tr_ind <- grep(paste(as.character(pos.snp), "\t", sep = ""), rent_trees)

   #get tree at focal site in newick format.
   rent_tree_newick <- sub(paste(as.character(pos.snp), "\t", sep = ""), "", rent_trees[rent_tr_ind])
   rent_tree <- read.tree(text = paste(rent_tree_newick, ";", sep = ""))

   n_chroms <- dim(haps.poly)[2]
   anc_tips <- which(rent_tree$tip.label %in% as.character(1:(n_chroms - n_der))) 
   der_tips <- which(as.numeric(rent_tree$tip.label) %in% (n_chroms - n_der + 1):n_chroms)
   anc_tree <- drop.tip(rent_tree, der_tips)
   der_tree <- drop.tip(rent_tree, anc_tips)

   
   trees.all[[i]] <- rent_tree
   trees.anc[[i]] <- anc_tree
   trees.der[[i]] <- der_tree
   bls.by.ndesc.all[[i]] <- branch.length.by.n.descendants(rent_tree)
   bls.by.ndesc.anc[[i]] <- branch.length.by.n.descendants(anc_tree)
   bls.by.ndesc.der[[i]] <- branch.length.by.n.descendants(der_tree)

   total.bls.all[i] <- total.bl(rent_tree)
   total.bls.anc[i] <- total.bl(anc_tree)
   total.bls.der[i] <- total.bl(der_tree)
   print(paste("round", as.character(i), "complete"))
}

#save.image(file = "rent_1700_trees_1KG_BiobankHeight_013018.RData")


n_ancs <- n_chroms- n_ders

worked <- !(is.na(n_ders) | n_ders == n_chroms | n_ders == 0)

mark.dat.all <- mark.dat
mark.dat <- mark.dat.all[worked,]

beta <- mark.dat$beta
Z <- mark.dat$tstat
rsID <- mark.dat$rsid

thetas.w <- thetas[worked]



#Goal: adjust time information using locus-wise thetas.
#fit neutral, waiting-time, and lineages-remaining estimators
#test for selection.
#Rent+ times are computed as t = 2*d/(L*theta)
#where d is a Hamming distance between haplotypes, L is the length of the haplotype
#in base pairs, and theta is estimated value of 4Nu (est. by Watterson).
#We'll multiply by theta estimated per site to get distance-matrix units,
#then rescale so that the mean tmrca is 2.

trees.der.rescale <- trees.der[worked]
trees.anc.rescale <- trees.anc[worked]
trees.all.rescale <- trees.all[worked]


for(i in 1:length(trees.all.rescale)){
	trees.der.rescale[[i]]$edge.length <- trees.der.rescale[[i]]$edge.length*thetas.w[i]
	trees.anc.rescale[[i]]$edge.length <- trees.anc.rescale[[i]]$edge.length*thetas.w[i]
	trees.all.rescale[[i]]$edge.length <- trees.all.rescale[[i]]$edge.length*thetas.w[i]
}



tmrcas <- numeric(0)
for(i in 1:length(trees.all.rescale)){
	tmrcas[i] <- max(branching.times(trees.all.rescale[[i]]))
}

for(i in 1:length(trees.all.rescale)){
	trees.der.rescale[[i]]$edge.length <- trees.der.rescale[[i]]$edge.length*2/mean(tmrcas)
	trees.anc.rescale[[i]]$edge.length <- trees.anc.rescale[[i]]$edge.length*2/mean(tmrcas)
	trees.all.rescale[[i]]$edge.length <- trees.all.rescale[[i]]$edge.length*2/mean(tmrcas)
}

time <- seq(0, 4, by = 0.001)

times.c <- list()
for(i in 1:length(trees.all.rescale)){
	times.c[[i]] <- trees_to_times(trees.all.rescale[[i]], trees.anc.rescale[[i]], trees.der.rescale[[i]], time, sure.alt.is.derived = FALSE)
}

lins.list <- list()
for(i in 1:length(trees.all.rescale)){
	lins.list[[i]] <- times_to_lins(times.c[[i]], time)
}



#neutral
trajs_neut <- matrix(nrow = length(time), ncol = length(trees.all.rescale))
vars_neut_bin <- matrix(nrow = length(time), ncol = length(trees.all.rescale))
vars_neut_post <- matrix(nrow = length(time), ncol = length(trees.all.rescale))
for(i in 1:length(trees.all.rescale)){	
	trajs_neut[,i] <- est_af_traj_neut(lins.list[[i]])
	vars_neut_bin[,i] <- est_af_var_neut_bin(lins.list[[i]])
	vars_neut_post[,i] <- est_af_var_neut_post(lins.list[[i]])
}
traj.phen.neut <- 2 * trajs_neut %*%  beta
var.phen.neut.bin <- 4 * vars_neut_bin %*% beta^2
var.phen.neut.post <- 4 * vars_neut_post %*% beta^2




#waiting time-based estimates and variance
trajs_est_wt_l1 <- matrix(nrow = length(time), ncol = length(trees.all.rescale))
trajs_var_wt_l1 <- matrix(nrow = length(time), ncol = length(trees.all.rescale))
for(i in 1:length(trees.all.rescale)){
	wt.estvar <- p_ests_wait(times.c[[i]], time, ell.ref = 1, ell.alt = 1)
	trajs_est_wt_l1[,i] <- wt.estvar[,1]	
	trajs_var_wt_l1[,i] <- wt.estvar[,2]	
}
traj.phen.wt_l1 <- 2 * trajs_est_wt_l1 %*%  beta 
var.phen.wt_l1 <- 4 * trajs_var_wt_l1 %*%  beta^2 
traj.phen.wt_l1[time == 0] <- traj.phen.neut[time == 0]
var.phen.wt_l1[time == 0] <- var.phen.neut.post[time == 0]








###########################################################
##############################################################################
##############################################################################



#time <- seq(0, 100000, by = 10)


#lineage number based estimates and variance.
trajs_mom_smoothtime <- matrix(nrow = length(time), ncol = length(trees.all.rescale))
for(i in 1:length(trees.all.rescale)){		
	trajs_mom_smoothtime[,i] <- est_af_traj_mom.smoothtime(lins.list[[i]], time)
}
traj.phen.mom_smoothtime <- 2 * trajs_mom_smoothtime %*%  beta 
traj.phen.mom_smoothtime[time == 0] <- traj.phen.neut[time == 0] 


#Smooth MOM variance
trajs_var_mom_smoothtime <- matrix(nrow = length(time), ncol = length(trees.all.rescale))
for(i in 1:length(trees.all.rescale)){		
	trajs_var_mom_smoothtime[,i] <- est_af_var_mom.smoothtime(lins.list[[i]], time*2*10000)
}
var.phen.mom_smoothtime <- 4 * trajs_var_mom_smoothtime %*%  beta^2 
traj.phen.mom_smoothtime[time == 0] <- traj.phen.neut[time == 0]
var.phen.mom_smoothtime[time == 0] <- var.phen.neut.bin[time == 0]















permdist.waittraj <- matrix(nrow = 1000, ncol = length(traj.phen.wt_l1))
for(i in 1:1000){
	permdist.waittraj[i,] <- 2 * trajs_est_wt_l1 %*%  sample(beta)
	permdist.waittraj[i,] <- permdist.waittraj[i,] + traj.phen.wt_l1[1] - permdist.waittraj[i,1]
}


permdist.neuttraj <- matrix(nrow = 1000, ncol = length(traj.phen.neut))
for(i in 1:1000){
	permdist.neuttraj[i,] <- 2 * trajs_neut %*%  sample(beta)
	permdist.neuttraj[i,] <- permdist.neuttraj[i,] + traj.phen.neut[1] - permdist.neuttraj[i,1]
}


permdist.waittraj.sf <- matrix(nrow = 1000, ncol = length(traj.phen.wt_l1))
for(i in 1:1000){
	signflip <- 2*rbinom(length(beta), 1, 0.5) - 1
	permdist.waittraj.sf[i,] <- 2 * trajs_est_wt_l1 %*%  (beta*signflip)
	permdist.waittraj.sf[i,] <- permdist.waittraj.sf[i,] + traj.phen.wt_l1[1] - permdist.waittraj.sf[i,1]
}

allele.freqs <- n_ders/182
plot(allele.freqs, beta)
cor(allele.freqs, abs(beta))

head(lins.list[[3]],30)
plot(time, trajs_est_wt_l1[,3], ylim = c(0,1), xlim = c(0, 500), type = "l", col = "gray")
lines(time, trajs_neut[,3])
mean(trajs_neut[ncol(trajs_neut),])
median(trajs_neut[1,])


for(i in 2:(ncol(trajs_est_wt_l1))){
	lines(time, trajs_est_wt_l1[,i], col = "grey")
}







tiff("heightlines.tif")
plot(time, traj.phen.neut, type = "l", xlab = "Time", ylab = "Estimated mean polygenic height", ylim = c(-1, 2), xlim = c(0, 100000))
#plot(time, traj.phen.neut, type = "l", xlab = "Time", ylab = "Estimated mean polygenic height", ylim = c(-2, 2.5), xlim = c(0, 300))
#plot(time, traj.phen.neut, type = "l", xlab = "Time", ylab = "Estimated mean polygenic height", ylim = c(-2, 2.5), xlim = c(0, .5))
polygon(c(time, rev(time)), c(traj.phen.neut + 1.96*sqrt(var.phen.neut.post), rev(traj.phen.neut - 1.96*sqrt(var.phen.neut.post))) , col = "grey", border = NA)
polygon(c(time, rev(time)), c(traj.phen.wt_l1 + 1.96*sqrt(var.phen.wt_l1), rev(traj.phen.wt_l1  - 1.96*sqrt(var.phen.wt_l1 ))) , col = "light blue", border = NA)
polygon(c(time, rev(time)), c(traj.phen.mom_smoothtime + 1.96*sqrt(var.phen.mom_smoothtime), rev(traj.phen.mom_smoothtime  - 1.96*sqrt(var.phen.mom_smoothtime ))) , col = "pink", border = NA)
lines(time, traj.phen.neut)
lines(time, traj.phen.wt_l1, col = "blue")
lines(time, traj.phen.mom_smoothtime, col = "red")

dev.off()



pal = c(rgb(0,1,1), rgb(1,0,0), rgb(.5,0,.5))
pal.int = c(rgb(0,1,1,.33), rgb(1,0,0,.33), rgb(.5,0,.5,.33))

tiff("height_lr_neutOnly.tif", units = "in", width = 6, height = 4, res = 300, compression = "lzw")
par(cex = 0.85)
plot(-time, traj.phen.neut, type = "l", xlab = "Time (estimated years)", ylab = "Estimated mean polygenic height (UK 1000G)", ylim = c(-2, 2.5), xlim = c(-0.5, 0), xaxt = "n")
axis(1, at = c(0, -1/6, -1/3, -1/2), labels = c("0", "100k", "200k", "300k"))
polygon(c(-time, rev(-time)), c(traj.phen.neut + 1.96*sqrt(var.phen.neut.post), rev(traj.phen.neut - 1.96*sqrt(var.phen.neut.post))) , col = pal.int[1], border = NA)
#polygon(c(-time, rev(-time)), c(traj.phen.wt_l1 + 1.96*sqrt(var.phen.wt_l1), rev(traj.phen.wt_l1  - 1.96*sqrt(var.phen.wt_l1 ))) , col = pal.int[2], border = NA)
#polygon(c(-time, rev(-time)), c(traj.phen.mom_smoothtime + 1.96*sqrt(var.phen.mom_smoothtime), rev(traj.phen.mom_smoothtime  - 1.96*sqrt(var.phen.mom_smoothtime ))) , col = "pink", border = NA)
lines(-time, traj.phen.neut, col = pal[1], lwd = 2)
#lines(-time, traj.phen.wt_l1, col = pal[2], lwd = 2)
#lines(time, traj.phen.mom_smoothtime, col = "red")
legend("topleft", lwd = 1, col = pal[1], legend = c("Neutral est."))
dev.off()


tiff("height_lr_2ests.tif", units = "in", width = 6, height = 4, res = 300, compression = "lzw")
par(cex = 0.85)
plot(-time, traj.phen.neut, type = "l", xlab = "Time (estimated years)", ylab = "Estimated mean polygenic height (UK 1000G)", ylim = c(-2, 2.5), xlim = c(-0.1, 0), xaxt = "n")
axis(1, at = c(0, -1/6, -1/3, -1/2), labels = c("0", "100k", "200k", "300k"))
polygon(c(-time, rev(-time)), c(traj.phen.neut + 1.96*sqrt(var.phen.neut.post), rev(traj.phen.neut - 1.96*sqrt(var.phen.neut.post))) , col = pal.int[1], border = NA)
polygon(c(-time, rev(-time)), c(traj.phen.wt_l1 + 1.96*sqrt(var.phen.wt_l1), rev(traj.phen.wt_l1  - 1.96*sqrt(var.phen.wt_l1 ))) , col = pal.int[2], border = NA)
#polygon(c(-time, rev(-time)), c(traj.phen.mom_smoothtime + 1.96*sqrt(var.phen.mom_smoothtime), rev(traj.phen.mom_smoothtime  - 1.96*sqrt(var.phen.mom_smoothtime ))) , col = "pink", border = NA)
lines(-time, traj.phen.neut, col = pal[1], lwd = 2)
lines(-time, traj.phen.wt_l1, col = pal[2], lwd = 2)
#lines(time, traj.phen.mom_smoothtime, col = "red")
legend("topleft", lwd = 1, col = pal[1:2], legend = c("Neutral est.", "Waiting-time est."))
dev.off()


#Qx test

timevec <- seq(5,20, by = 5)/100
Va <- add.var(beta, trajs_neut[1,])
time.ests.neut <- cumsum(est.times.seq(trajs_neut[time %in% c(max(time[time < min(timevec)]), timevec),]))
fm.neut <- fmat.time(time.ests.neut)

Qx_neut <- compute_Qx_traj(traj.phen.neut[time %in% timevec] - traj.phen.neut[time == max(time[time < min(timevec)])]  , fm.neut, Va, f.inverted = FALSE)
1- pchisq(Qx_neut, length(timevec) - 1)







for(i in 1:1000){
	lines(time, permdist.waittraj[i,], col = "grey")
}

for(i in 1:1000){
	lines(time, permdist.neuttraj[i,], col = "grey")
}


for(i in 1:1000){
	lines(time, permdist.waittraj.sf[i,], col = "grey")
}




time.orig <- seq(0.01, 20, by = .01)
times.c.orig <- list()
for(i in 1:length(trees.all)){
	times.c.orig[[i]] <- trees_to_times(trees.all[[i]], trees.anc[[1]], trees.der[[i]], time.orig, sure.alt.is.derived = FALSE)
}

lins.list.orig <- list()
for(i in 1:length(trees.all)){
	lins.list.orig[[i]] <- times_to_lins(times.c.orig[[i]], time.orig)
}


#neutral
trajs_neut.orig <- matrix(nrow = length(time.orig), ncol = length(trees.all.rescale))
vars_neut_bin.orig <- matrix(nrow = length(time.orig), ncol = length(trees.all.rescale))
vars_neut_post.orig <- matrix(nrow = length(time.orig), ncol = length(trees.all.rescale))
for(i in 1:length(trees.all.rescale)){	
	trajs_neut.orig[,i] <- est_af_traj_neut(lins.list.orig[[i]])
	vars_neut_bin.orig[,i] <- est_af_var_neut_bin(lins.list.orig[[i]])
	vars_neut_post.orig[,i] <- est_af_var_neut_post(lins.list.orig[[i]])
}
traj.phen.neut.orig <- 2 * trajs_neut.orig %*%  beta
var.phen.neut.bin.orig <- 4 * vars_neut_bin.orig %*% beta^2
var.phen.neut.post.orig <- 4 * vars_neut_post.orig %*% beta^2

plot(time.orig, traj.phen.neut.orig, type = "l", xlim = c(0,4),xlab = "Time", ylab = "Estimated mean polygenic height", xaxt = "n")




#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
plot(coal.times(trees.anc[[2]]))
trees.anc[[1]]$edge.length

dert <- coal.times(trees.der[[1]])
anct <- coal.times(trees.anc[[1]])
allt <- coal.times(trees.all[[1]])


time <- seq(0.005, 100.005, by = 0.01)
trees.all.w <- trees.all[worked & n_ders > 1]
trees.anc.w <- trees.anc[worked & n_ders > 1]
trees.der.w <- trees.der[worked & n_ders > 1]

trajs <- matrix(nrow = length(time), ncol = length(trees.all.w))

for(i in 1:length(trees.all.w)){
	trajs[,i] <- est_alt_af_traj(trees.all.w[[i]], trees.anc.w[[i]], trees.der.w[[i]], time)
}

traj.ht <- trajs %*%  beta[worked & n_ders > 1]

trajs.var <- matrix(nrow = length(time), ncol = length(trees.all.w))

for(i in 1:length(trees.all.w)){
	trajs.var[,i] <- est_alt_af_var(trees.all.w[[i]], trees.anc.w[[i]], trees.der.w[[i]], time)
}

traj.ht.var.rent <- trajs.var %*%  beta[worked & n_ders > 1]^2

tmrcas <- numeric(0)
for(i in 1:length(trees.all.w)){
	tmrcas[i] <- max(branching.times(trees.all.w[[i]]))
}


#mean tmrcas = 604.8. Quick and dirty adjustment is to divide times by 604.8/4 = 151.2
#to get units of N generations, then multiply by N=10000 * 29-year generation time to get
#approximate years 

time.adj <- (time / 151.2) * 10000 * 29

pdf("Height_polyscore_UK_rent_withCredInt.pdf")
plot(rev(time.adj) - max(time.adj), traj.ht, type = "l", xlim = c(-30000, 0 ), xlab = "Time (thousands of years before present)", ylab = "Estimated mean polygenic height", ylim = c(0.15, 0.7), xaxt = "n")

axis(1, at = c(0, -10000, -20000, -30000, -40000, -50000), labels = c("0", "10", "20", "30", "40", "50"))

#lines(time, traj.phen.rent + 1.96*sqrt(traj.phen.var.rent), col = pal[3])
#lines(time, traj.phen.rent - 1.96*sqrt(traj.phen.var.rent), col = pal[3])
polygon(c(rev(time.adj), time.adj) - max(time.adj), c(traj.ht + 1.96*sqrt(traj.ht.var.rent), rev(traj.ht - 1.96*sqrt(traj.ht.var.rent)) ), border = NA, col = "grey")
lines(rev(time.adj) - max(time.adj), traj.ht)
dev.off()




n.perms <- 1000
perms <- matrix(ncol = n.perms, nrow = length(traj.ht))
for(i in 1:n.perms){
	perms[,i] <- trajs %*%  sample(beta[worked & n_ders > 1])
}


perms.sign <- matrix(ncol = n.perms, nrow = length(traj.ht))
for(i in 1:n.perms){
	shuff.beta <- beta[worked & n_ders > 1] * (2*rbinom(length(beta[worked & n_ders > 1]), 1, 0.5) - 1)
	perms.sign[,i] <- trajs %*%  shuff.beta
}


tiff("Height_polyscore_UK_rent_withPerms.tif")
plot(time, traj.ht, type = "l", lwd = 2, ylim = c(0,1), ylab = "estimated genetic height")
for(i in 1:n.perms){
	lines(time, perms[,i] + traj.ht[1] - perms[1,i], col = "gray", lwd = 0.5)
}
lines(time, traj.ht, lwd = 2)
dev.off()


tiff("Height_polyscore_UK_rent_withSignPerms.tif")
plot(time, traj.ht, type = "l", lwd = 2, ylim = c(0,1), ylab = "estimated genetic height")
for(i in 1:n.perms){
	lines(time, perms.sign[,i] + traj.ht[1] - perms.sign[1,i], col = "gray", lwd = 0.5)
}
lines(time, traj.ht, lwd = 2)
dev.off()



#Plot that marks changes that are larger than typical changes between timepoints in 
#permutation distribution.
perms.diffs <- diff(perms.sign)
traj.ht.diffs <- diff(traj.ht)
lq.diffs <- apply(perms.diffs, 1, quantile, probs = 0.025)
hq.diffs <- apply(perms.diffs, 1, quantile, probs = 0.975)

plot(time, traj.ht, type = "l", lwd = 2, ylim = c(0.15,0.7), ylab = "estimated genetic height")
rug(time[which(traj.ht.diffs > hq.diffs | traj.ht.diffs < lq.diffs)])












time.ests.neut <- cumsum(est.times.seq(trajs_neut[time %in% c(max(time[time < min(timevec)]),timevec),]))
	fm.neut <- fmat.time(time.ests.neut)
	Qx_neut <- compute_Qx_traj(traj.phen.neut[time %in% timevec] - traj.phen.neut[time == max(time[time < min(timevec)])]  , fm.neut, Va, f.inverted = FALSE)

