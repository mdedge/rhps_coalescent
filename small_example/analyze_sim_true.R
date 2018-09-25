
#Doc Edge, 6/26/18

#once a trait dataset is loaded, analyze all the true trees.



#neutral
trajs_neut <- matrix(nrow = length(time), ncol = length(ms_trees_list))
vars_neut_bin <- matrix(nrow = length(time), ncol = length(ms_trees_list))
#vars_neut_post <- matrix(nrow = length(time), ncol = length(ms_trees_list))
for(i in 1:length(ms_trees_list)){	
	trajs_neut[,i] <- est_af_traj_neut(lins.list[[i]])
	vars_neut_bin[,i] <- est_af_var_neut_bin(lins.list[[i]])
#	vars_neut_post[,i] <- est_af_var_neut_post(lins.list[[i]])
}
trajs_neut[time == 0,] <- (n_ders/n_chroms) #in the present, just use sample allele frequency.
#This is the same as the output of est_af_traj_neut() if no coalescent times get rounded to 0.
vars_neut_bin[time == 0,] <- (n_ders/n_chroms)*(1 - (n_ders/n_chroms))/n_chroms
traj.phen.neut <- 2 * trajs_neut %*%  eff_sizes 
var.phen.neut.bin <- 4 * vars_neut_bin %*% eff_sizes^2
#var.phen.neut.post <- 4 * vars_neut_post %*% eff_sizes^2



#Method of moments from smoothed coalescent time estimates.
trajs_mom_smoothtime <- matrix(nrow = length(time), ncol = length(ms_trees_list))
trajs_var_mom_smoothtime <- matrix(nrow = length(time), ncol = length(ms_trees_list))
for(i in 1:length(ms_trees_list)){		
	trajs_mom_smoothtime[,i] <- est_af_traj_mom.smoothtime(lins.list[[i]], time)
	trajs_var_mom_smoothtime[,i] <- est_af_var_mom.smoothtime(lins.list[[i]], time*2*N)
}
traj.phen.mom_smoothtime <- 2 * trajs_mom_smoothtime %*%  eff_sizes 
var.phen.mom_smoothtime <- 4 * trajs_var_mom_smoothtime %*%  eff_sizes^2 
traj.phen.mom_smoothtime[time == 0] <- traj.phen.neut[time == 0]
var.phen.mom_smoothtime[time == 0] <- var.phen.neut.bin[time == 0]


#waiting time-based estimates and variance
trajs_est_wt_l1 <- matrix(nrow = length(time), ncol = length(ms_trees_list))
trajs_var_wt_l1 <- matrix(nrow = length(time), ncol = length(ms_trees_list))
for(i in 1:length(ms_trees_list)){
	wt.estvar <- p_ests_wait(times.c[[i]], time, ell.ref = 1, ell.alt = 1)
	trajs_est_wt_l1[,i] <- wt.estvar[,1]	
	trajs_var_wt_l1[,i] <- wt.estvar[,2]	
}
traj.phen.wt_l1 <- 2 * trajs_est_wt_l1 %*%  eff_sizes 
var.phen.wt_l1 <- 4 * trajs_var_wt_l1 %*%  eff_sizes^2 
traj.phen.wt_l1[time == 0] <- traj.phen.neut[time == 0]
var.phen.wt_l1[time == 0] <- var.phen.neut.bin[time == 0]


#flat line backward estimator
traj.flatline <- rep(traj.phen.neut[time == 0], length(traj.phen.neut))

#straight line to ancestral value at 2 coal units.
traj.straightline <- traj.phen.neut[time == 0] - (traj.phen.neut[time == 0]/2) * time
traj.straightline[time > 2] <- 0


true.per.time <- numeric(0)
for(j in 1:length(time)){
	true.per.time[j] <- phen.traj[which.min(abs(pt.time - time[j]))]
}

mat.trajs <- matrixify.list.of.trajs(trajs)

true.afs.per.time <- matrix(nrow = length(time), ncol = n.loci)
for(j in 1:length(time)){
	true.afs.per.time[j,] <- mat.trajs[which.min(abs(pt.time - time[j])),]
}


#save errors, unscaled and scaled by estimated se, for each method at all times.
err.neut <- traj.phen.neut - true.per.time
err.smoothmom <- traj.phen.mom_smoothtime - true.per.time
err.wt_l1 <- traj.phen.wt_l1 - true.per.time
err.flat <- traj.flatline - true.per.time
err.straight <- traj.straightline - true.per.time

err.neut.std.bin <- err.neut / sqrt(var.phen.neut.bin)
err.smoothmom.std <- err.smoothmom / sqrt(var.phen.mom_smoothtime)
err.wt_l1.std <- err.wt_l1 / sqrt(var.phen.wt_l1)

err.mat <- cbind(err.neut, err.smoothmom, err.wt_l1, err.flat, err.straight)
err.mat.std <- cbind(err.neut.std.bin, err.smoothmom.std, err.wt_l1.std)

err.array[,,iter] <- err.mat
err.std.array[,,iter] <- err.mat.std
mat.true.phentrajs[,iter] <- true.per.time

#Tests
qxtest_mat[iter,1:3] <- Qx_test(trajs_neut[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat[iter,4] <- Qx_test(trajs_neut[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat[iter,5:7] <- Qx_test(trajs_mom_smoothtime[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat[iter,8] <- Qx_test(trajs_mom_smoothtime[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat[iter,9:11] <- Qx_test(trajs_est_wt_l1[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat[iter,12] <- Qx_test(trajs_est_wt_l1[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]

print("true tree T_X statistic, number of timepoints, and permutation p")
print(Qx_test(trajs_neut[time %in% ((0:10)/100),], eff_sizes, perms = 10000))




