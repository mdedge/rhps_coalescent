#Goal: Analyze the simulation pulses to show approximate tSDS power


N <- 10000
herit <- 1
out_dir <- "out_061318/"

#nsamps <- 20
sel.intenses <- .01
n.locis <- c(100)
n_chromss <- c(20, 50, 100, 200, 500, 1000)
t.offs <- c(0, .01, .02, .03, .04, .05, .06, .07, .08, .09)
ts <- t.offs + .005
phen_nums <- c(1:100) #1 number for each rep we want to do at each combination of parameters

#traj.fn <- "temp_out1/temp.txt"
#msout.fn <- "temp_out1/ms_out.txt"
#rent_in_fn <- "temp_out1/rent_in.txt"


helper_fn <- "../../helper_functions_coal_sel.R"
source(helper_fn) #read in helper functions and load ape package


#ms_dir <- "../../msseldir/"
#rentplus_fn <- "../../RentPlus.jar"
#len_hap <- 200000 #the length (in base pairs) of the haplotype -- short because 
#we are not worrying about recombination--just need the sel site.
#sel_site <- 100000 #the position of the selected site in the haplotype
#u <- 2e-8 #the neutral mutation rate per base pair/generation
#r <- 2.5e-8 #the recombination rate per base pair/generation
options(scipen = 999) #disable scientific notation so that parameters passed to ms are given as numbers
#sd.trait <- 1

#time <- c(seq(0, 4, by = 0.001))
pars <- expand.grid(n_chromss, phen_nums, t.offs, sel.intenses, n.locis)
tSDS_mat <- matrix(nrow = nrow(pars), ncol = 6)



for(k in 1:dim(pars)[1]){
	n_chroms <- pars[k,1]
	sel.intens <- sel.intenses
	n.loci <- n.locis	
	t <- pars[k,3] + .005
	t.off <- pars[k,3] 
	phen_num <- pars[k,2]
	#if((k - 1) %% length(n_chromss) == 0){
	#	print(k)
	#	print(paste("simulating polygenic score trajectory", as.character(t.off)))
	#	source("pheno_sim_aftrajs.R")
	#	fn_str <- paste("trajs_loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")
	#save(list = c("trajs","time", "eff_sizes"),file = paste(out_dir, fn_str, ".RData", sep = ""))
	#}	
	#source("pheno_sim_trees.R")
	
	load(paste("out_061318/trees_loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), ".RData", sep = ""))

	#split into ancestral and derived trees and retrieve coalescence times.
	anc_trees_ms <- list()
	der_trees_ms <- list()

	for(i in 1:n.loci){
		anc_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
		der_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
		anc_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], der_tips_ms)
		der_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], anc_tips_ms)
	}

	ts <- tSDS_analogue(anc_trees_ms, der_trees_ms)
	tSDS.stat <- sum(ts * sign(eff_sizes), na.rm = TRUE)
	tSDS.stat.beta <- sum(ts * eff_sizes, na.rm = TRUE)
	tSDS.stat.sign <- sum(sign(ts) * sign(eff_sizes), na.rm = TRUE)
	ts.perm <- numeric(10000)	
	ts.perm.beta <- numeric(10000)
	ts.perm.sign <- numeric(10000)
	for(rep in 1:10000){
		perm.eff <- sample(eff_sizes)		
		ts.perm[rep] <- sum(ts * sign(perm.eff), na.rm = TRUE)
		ts.perm.beta[rep] <- sum(ts * perm.eff, na.rm = TRUE)
		ts.perm.sign[rep] <- sum(sign(ts) * sign(perm.eff), na.rm = TRUE)
	}
	sds.p <- mean(abs(ts.perm) >= abs(tSDS.stat))
	sds.p.beta <- mean(abs(ts.perm.beta) >= abs(tSDS.stat.beta))
	sds.p.sign <- mean(abs(ts.perm.sign) >= abs(tSDS.stat.sign))
	tSDS_mat[k,] <- c(tSDS.stat, sds.p, tSDS.stat.beta, sds.p.beta, tSDS.stat.sign, sds.p.sign)

	print(pars[k,])
	print(tSDS_mat[k,])
	print(paste("trial", as.character(k), "complete."))
}

save("tSDS_mat", file = "tSDS_mat.RData")

aggregate(tSDS_mat[,2], FUN = function(x){mean(x<.05)}, by = list(pars[,3], pars[,1]))
aggregate(tSDS_mat[,4], FUN = function(x){mean(x<.05)}, by = list(pars[,3], pars[,1]))
aggregate(tSDS_mat[,6], FUN = function(x){mean(x<.05)}, by = list(pars[,3], pars[,1]))






pows <- aggregate(tSDS_mat[,4], FUN = function(x){mean(x<.05)}, by = list(pars[,3], pars[,1]))

bmp("sds_pow.bmp")
pal <- c("purple", "pink", "blue", "orange", "brown", "grey")
plot(pows[pows[,2] == 1000,1], pows[pows[,2] == 1000,3], ylim = c(0,1), type = "o", col = pal[1],
xlim = c(0,.1), ylab = "Power", xlab = "End of selection interval (coalescent units before present)", bty = "n", las = 1, pch = 20)
points(pows[pows[,2] == 500,1], pows[pows[,2] == 500,3], col = pal[2], type = "o", pch = 20)
points(pows[pows[,2] == 200,1], pows[pows[,2] == 200,3], col = pal[3], type = "o", pch = 20)
points(pows[pows[,2] == 100,1], pows[pows[,2] == 100,3], col = pal[4], type = "o", pch = 20)
points(pows[pows[,2] == 50,1], pows[pows[,2] == 50,3], col = pal[5], type = "o", pch = 20)
points(pows[pows[,2] == 20,1], pows[pows[,2] == 20,3], col = pal[6], type = "o", pch = 20)
legend("topright", bty = "n", lty = 1, col = pal, legend = c(" n = 1000","n = 500","n = 200","n = 100","n = 50","n = 20"))
dev.off()






#redo with empirical means and SDs from simulations.
tSDS_mat_emp <- matrix(nrow = nrow(pars), ncol = 6)


e <- new.env()
load("../neutral_sims_for_tSDS/neutral_SDS_sims_summary.RData", envir = e)
dafs.sim <- unique(e$pars[,2])

for(k in 1:dim(pars)[1]){
	n_chroms <- pars[k,1]
	sel.intens <- sel.intenses
	n.loci <- n.locis	
	t <- pars[k,3] + .005
	t.off <- pars[k,3] 
	phen_num <- pars[k,2]

	if((k - 1) %% length(n_chromss) == 0){
		#print(paste("simulating polygenic score trajectory", as.character(t.off)))
		#source("pheno_sim_aftrajs.R")
		fn_str <- paste("trajs_loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")
		traj.env <- new.env()		
		load(file = paste(out_dir, fn_str, ".RData", sep = ""), envir = traj.env)
		curr.daf <- sapply(traj.env$trajs, FUN = function(x){x[1,2]})
		nearest.dafs <- sapply(curr.daf, FUN = function(x, y){y[which.min(abs(y-x))[1]]}, y = dafs.sim)
	}		
	
	
	load(paste("out_061318/trees_loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), ".RData", sep = ""))

	#split into ancestral and derived trees and retrieve coalescence times.
	anc_trees_ms <- list()
	der_trees_ms <- list()

	for(i in 1:n.loci){
		anc_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character(1:(n_chroms - n_ders[i])))
		der_tips_ms <- which(ms_trees_list[[i]]$tip.label %in% as.character((n_chroms - n_ders[i] + 1):n_chroms))
		anc_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], der_tips_ms)
		der_trees_ms[[i]] <- drop.tip(ms_trees_list[[i]], anc_tips_ms)
	}
	
	neutpars.ind <- apply(cbind(n_chroms, nearest.dafs), 1, function(x, y){which(y[,1] == x[1] & y[,2] == x[2])}, y = e$pars)

	means.daf <- e$mean.sds[neutpars.ind]
	sds.daf <- e$sd.sds[neutpars.ind]
	ts <- tSDS_analogue_meansd(anc_trees_ms, der_trees_ms, means.daf, sds.daf)
	tSDS.stat <- sum(ts * sign(eff_sizes), na.rm = TRUE)
	tSDS.stat.beta <- sum(ts * eff_sizes, na.rm = TRUE)
	tSDS.stat.sign <- sum(sign(ts) * sign(eff_sizes), na.rm = TRUE)
	ts.perm <- numeric(10000)	
	ts.perm.beta <- numeric(10000)
	ts.perm.sign <- numeric(10000)
	for(rep in 1:10000){
		perm.eff <- sample(eff_sizes)		
		ts.perm[rep] <- sum(ts * sign(perm.eff), na.rm = TRUE)
		ts.perm.beta[rep] <- sum(ts * perm.eff, na.rm = TRUE)
		ts.perm.sign[rep] <- sum(sign(ts) * sign(perm.eff), na.rm = TRUE)
	}
	sds.p <- mean(abs(ts.perm) >= abs(tSDS.stat))
	sds.p.beta <- mean(abs(ts.perm.beta) >= abs(tSDS.stat.beta))
	sds.p.sign <- mean(abs(ts.perm.sign) >= abs(tSDS.stat.sign))
	tSDS_mat_emp[k,] <- c(tSDS.stat, sds.p, tSDS.stat.beta, sds.p.beta, tSDS.stat.sign, sds.p.sign)

	print(pars[k,])
	print(tSDS_mat[k,])
	print(paste("trial", as.character(k), "complete."))
}


save("tSDS_mat_emp", file = "tSDS_mat_emp.RData")




aggregate(tSDS_mat_emp[,2], FUN = function(x){mean(x<.05, na.rm = TRUE)}, by = list(pars[,3], pars[,1]))
aggregate(tSDS_mat_emp[,4], FUN = function(x){mean(x<.05, na.rm = TRUE)}, by = list(pars[,3], pars[,1]))
aggregate(tSDS_mat_emp[,6], FUN = function(x){mean(x<.05, na.rm = TRUE)}, by = list(pars[,3], pars[,1]))






pows <- aggregate(tSDS_mat_emp[,4], FUN = function(x){mean(x<.05)}, by = list(pars[,3], pars[,1]))

pdf("sds_pow_emp.pdf", width = 6, height = 4.5)
pal <- c("purple", "pink", "blue", "orange", "brown", "grey")
plot(pows[pows[,2] == 1000,1], pows[pows[,2] == 1000,3], ylim = c(0,1), type = "o", col = pal[1],
xlim = c(0,.1), ylab = "Power", xlab = "End of selection interval (coalescent units before present)", bty = "n", las = 1, pch = 20)
points(pows[pows[,2] == 500,1], pows[pows[,2] == 500,3], col = pal[2], type = "o", pch = 20)
points(pows[pows[,2] == 200,1], pows[pows[,2] == 200,3], col = pal[3], type = "o", pch = 20)
points(pows[pows[,2] == 100,1], pows[pows[,2] == 100,3], col = pal[4], type = "o", pch = 20)
points(pows[pows[,2] == 50,1], pows[pows[,2] == 50,3], col = pal[5], type = "o", pch = 20)
points(pows[pows[,2] == 20,1], pows[pows[,2] == 20,3], col = pal[6], type = "o", pch = 20)
legend("topright", bty = "n", lty = 1, col = pal, legend = c(" n = 1000","n = 500","n = 200","n = 100","n = 50","n = 20"))
dev.off()





#add T_X power
N <- 10000
herit <- 1
n_chromss <- c(20, 50, 100, 200, 500, 1000)
t.offs <- c(0, .01, .02, .03, .04, .05, .06, .07, .08, .09)
ts <- t.offs + .005
phen_nums <- 1:100


e <- new.env()
fn_str <- paste("qxmat_ton", as.character(ts[1]), "_toff", as.character(t.offs[1]), "phen_nums_", as.character(min(phen_nums)), "_to_", max(phen_nums), sep = "")
load(paste(fn_str, ".RData", sep = ""), envir = e)
all_pars <- e$pars
all_qxtest_mat <- e$qxtest_mat



for(i in 2:length(t.offs)){
	e <- new.env()
	fn_str <- paste("qxmat_ton", as.character(ts[i]), "_toff", as.character(t.offs[i]), "phen_nums_", as.character(min(phen_nums)), "_to_", max(phen_nums), sep = "")
	load(paste(fn_str, ".RData", sep = ""), envir = e)
	all_pars <- rbind(all_pars, e$pars)
	all_qxtest_mat <- rbind(all_qxtest_mat, e$qxtest_mat)
}

e$pows <- aggregate(all_qxtest_mat[,4], all_pars[,c(5,1)], function(x){mean(x < .05)})




pdf("sds_pow_emp_wTx.pdf", width = 6, height = 4.5)
pal <- c("purple", "pink", "blue", "orange", "brown", "grey")
plot(pows[pows[,2] == 1000,1], pows[pows[,2] == 1000,3], ylim = c(0,1), type = "o", col = pal[1],
xlim = c(0,.1), ylab = "Power", xlab = "End of selection interval (coalescent units before present)", bty = "n", las = 1, pch = 20)
points(pows[pows[,2] == 500,1], pows[pows[,2] == 500,3], col = pal[2], type = "o", pch = 20)
points(pows[pows[,2] == 200,1], pows[pows[,2] == 200,3], col = pal[3], type = "o", pch = 20)
points(pows[pows[,2] == 100,1], pows[pows[,2] == 100,3], col = pal[4], type = "o", pch = 20)
points(pows[pows[,2] == 50,1], pows[pows[,2] == 50,3], col = pal[5], type = "o", pch = 20)
points(pows[pows[,2] == 20,1], pows[pows[,2] == 20,3], col = pal[6], type = "o", pch = 20)
points(e$pows[e$pows[,2] == 1000,1], e$pows[e$pows[,2] == 1000,3], col = "black", type = "o", pch = 20, lty = 2)
legend("topright", bty = "n", lty = c(1,1,1,1,1,1,2), col = c(pal,"black"), legend = c("n = 1000","n = 500","n = 200","n = 100","n = 50","n = 20",expression(paste(italic(T[X])," w. n = 1000")) ))
dev.off()

bmp("sds_pow_emp_wTx.bmp")
pal <- c("purple", "pink", "blue", "orange", "brown", "grey")
plot(pows[pows[,2] == 1000,1], pows[pows[,2] == 1000,3], ylim = c(0,1), type = "o", col = pal[1],
xlim = c(0,.1), ylab = "Power", xlab = "End of selection interval (coalescent units before present)", bty = "n", las = 1, pch = 20)
points(pows[pows[,2] == 500,1], pows[pows[,2] == 500,3], col = pal[2], type = "o", pch = 20)
points(pows[pows[,2] == 200,1], pows[pows[,2] == 200,3], col = pal[3], type = "o", pch = 20)
points(pows[pows[,2] == 100,1], pows[pows[,2] == 100,3], col = pal[4], type = "o", pch = 20)
points(pows[pows[,2] == 50,1], pows[pows[,2] == 50,3], col = pal[5], type = "o", pch = 20)
points(pows[pows[,2] == 20,1], pows[pows[,2] == 20,3], col = pal[6], type = "o", pch = 20)
points(e$pows[e$pows[,2] == 1000,1], e$pows[e$pows[,2] == 1000,3], col = "black", type = "o", pch = 20, lty = 2)
legend("topright", bty = "n", lty = c(1,1,1,1,1,1,2), col = c(pal,"black"), legend = c("n = 1000","n = 500","n = 200","n = 100","n = 50","n = 20",expression(paste(italic(T[X])," w. n = 1000")) ))
dev.off()


