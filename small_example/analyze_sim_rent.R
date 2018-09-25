
#Doc Edge, 6/26/18

#once a trait dataset is loaded, analyze all the rent trees.


#neutral - rent
trajs_neut_rent <- matrix(nrow = length(time), ncol = length(rent_trees_list))
vars_neut_bin_rent <- matrix(nrow = length(time), ncol = length(rent_trees_list))
#vars_neut_post_rent <- matrix(nrow = length(time), ncol = length(rent_trees_list))
for(i in 1:length(rent_trees_list)){	
	trajs_neut_rent[,i] <- est_af_traj_neut(lins.list.rent[[i]])
	vars_neut_bin_rent[,i] <- est_af_var_neut_bin(lins.list.rent[[i]])
#	vars_neut_post_rent[,i] <- est_af_var_neut_post(lins.list.rent[[i]])
}
trajs_neut_rent[time == 0,] <- (n_ders/n_chroms) #in the present, just use sample allele frequency.
#This is the same as the output of est_af_traj_neut() if no coalescent times get rounded to 0.
vars_neut_bin_rent[time == 0,] <- (n_ders/n_chroms)*(1 - (n_ders/n_chroms))/n_chroms
traj.phen.neut.rent <- 2 * trajs_neut_rent %*%  eff_sizes 
var.phen.neut.bin.rent <- 4 * vars_neut_bin_rent %*% eff_sizes^2
#var.phen.neut.post.rent <- 4 * vars_neut_post_rent %*% eff_sizes^2




#Method of moments from smoothed coalescent time estimates---RENT.
trajs_mom_smoothtime_rent <- matrix(nrow = length(time), ncol = length(rent_trees_list))
trajs_var_mom_smoothtime_rent <- matrix(nrow = length(time), ncol = length(rent_trees_list))
for(i in 1:length(rent_trees_list)){		
	trajs_mom_smoothtime_rent[,i] <- est_af_traj_mom.smoothtime(lins.list.rent[[i]], time)
	trajs_var_mom_smoothtime_rent[,i] <- est_af_var_mom.smoothtime(lins.list.rent[[i]], time*2*N)
}
traj.phen.mom_smoothtime_rent <- 2 * trajs_mom_smoothtime_rent %*%  eff_sizes 
var.phen.mom_smoothtime_rent <- 4 * trajs_var_mom_smoothtime_rent %*%  eff_sizes^2 
traj.phen.mom_smoothtime_rent[time == 0] <- traj.phen.neut.rent[time == 0]
var.phen.mom_smoothtime_rent[time == 0] <- var.phen.neut.bin.rent[time == 0]





#waiting time-based estimates and variance---RENT
trajs_est_wt_l1_rent <- matrix(nrow = length(time), ncol = length(rent_trees_list))
trajs_var_wt_l1_rent <- matrix(nrow = length(time), ncol = length(rent_trees_list))
for(i in 1:length(rent_trees_list)){
	wt.estvar.rent <- p_ests_wait(times.c.rent[[i]], time, ell.ref = 1, ell.alt = 1)
	trajs_est_wt_l1_rent[,i] <- wt.estvar.rent[,1]	
	trajs_var_wt_l1_rent[,i] <- wt.estvar.rent[,2]	
}
traj.phen.wt_l1_rent <- 2 * trajs_est_wt_l1_rent %*%  eff_sizes 
var.phen.wt_l1_rent <- 4 * trajs_var_wt_l1_rent %*%  eff_sizes^2 
traj.phen.wt_l1_rent[time == 0] <- traj.phen.neut.rent[time == 0]
var.phen.wt_l1_rent[time == 0] <- var.phen.neut.bin.rent[time == 0]



true.per.time <- numeric(0)
for(j in 1:length(time)){
	true.per.time[j] <- phen.traj[which.min(abs(pt.time - time[j]))]
}

mat.trajs <- matrixify.list.of.trajs(trajs)

true.afs.per.time <- matrix(nrow = length(time), ncol = n.loci)
for(j in 1:length(time)){
	true.afs.per.time[j,] <- mat.trajs[which.min(abs(pt.time - time[j])),]
}






#save errors, unscaled and scaled by estimated se, for each method at all times---RENT.

err.neut.rent <- traj.phen.neut.rent - true.per.time
err.smoothmom.rent <- traj.phen.mom_smoothtime_rent - true.per.time
err.wt_l1.rent <- traj.phen.wt_l1_rent - true.per.time


err.neut.std.bin.rent <- err.neut.rent / sqrt(var.phen.neut.bin.rent)
err.smoothmom.std.rent <- err.smoothmom.rent / sqrt(var.phen.mom_smoothtime_rent)
err.wt_l1.std.rent <- err.wt_l1.rent / sqrt(var.phen.wt_l1_rent)

err.mat.rent <- cbind(err.neut.rent, err.smoothmom.rent, err.wt_l1.rent)
err.mat.std.rent <- cbind(err.neut.std.bin.rent, err.smoothmom.std.rent, err.wt_l1.std.rent)

err.array.rent[,,iter] <- err.mat.rent
err.std.array.rent[,,iter] <- err.mat.std.rent
mat.true.phentrajs[,iter] <- true.per.time


#Tests
qxtest_mat_rent[iter,1:3] <- Qx_test(trajs_neut_rent[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_rent[iter,4] <- Qx_test(trajs_neut_rent[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat_rent[iter,5:7] <- Qx_test(trajs_mom_smoothtime_rent[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_rent[iter,8] <- Qx_test(trajs_mom_smoothtime_rent[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat_rent[iter,9:11] <- Qx_test(trajs_est_wt_l1_rent[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_rent[iter,12] <- Qx_test(trajs_est_wt_l1_rent[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]

print("rent tree T_X statistic, number of timepoints, and permutation p")
print(Qx_test(trajs_neut_rent[time %in% ((0:10)/100),], eff_sizes, perms = 10000))







