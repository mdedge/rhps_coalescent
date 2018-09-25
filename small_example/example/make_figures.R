#Doc Edge, 7/5/18
#Extract data from simulations and make figures.
#to be run in maintext_sims_rent

load("analyzed_trees.RData")

weak0204.truetraj <- mat.true.phentrajs
weak0204.err.array <- err.array
weak0204.err.std.array <- err.std.array
weak0204.qxtest_mat <- qxtest_mat

weak0204.rent.err.array <- err.array.rent
weak0204.rent.err.std.array <- err.std.array.rent
weak0204.rent.qxtest_mat <- qxtest_mat_rent


weak0204.avg.err.mat <- apply(weak0204.err.array, 1:2, mean, na.rm = TRUE)
weak0204.avg.sq.err.mat <- apply(weak0204.err.array^2, 1:2, mean, na.rm = TRUE)

weak0204.rent.avg.err.mat <- apply(weak0204.rent.err.array, 1:2, mean, na.rm = TRUE)
weak0204.rent.avg.sq.err.mat <- apply(weak0204.rent.err.array^2, 1:2, mean, na.rm = TRUE)



#Confidence intervals
prop.covered <- function(vec, q = 1.96, ...){
	mean(abs(vec) <= q, ...)
}

weak0204.covered.mat <- apply(weak0204.err.std.array, 1:2, prop.covered, na.rm = TRUE)

plot(time, weak0204.covered.mat[,1], type = "l", xlim = c(0,1/2), ylim = c(.8,1))
lines(time, weak0204.covered.mat[,2], col = "orange")
lines(time, weak0204.covered.mat[,3], col = "blue")
lines(c(0,10), c(.95,.95), col = "grey", lty = 2)

weak0204.rent.covered.mat <- apply(weak0204.rent.err.std.array, 1:2, prop.covered, na.rm = TRUE)

plot(time, weak0204.rent.covered.mat[,1], type = "l", xlim = c(0,1/2), ylim = c(.8,1))
lines(time, weak0204.rent.covered.mat[,2], col = "orange")
lines(time, weak0204.rent.covered.mat[,3], col = "blue")
lines(c(0,10), c(.95,.95), col = "grey", lty = 2)



#tests---power / Type I error
apply(weak0204.qxtest_mat[,c(3,4,7,8,11,12)], 2, function(x){mean(x<.05)})
apply(weak0204.rent.qxtest_mat[,c(3,4,7,8,11,12)], 2, function(x){mean(x<.05)})


pal <- c('#1b9e77','#d95f02','#7570b3','#e7298a')
#############################################################
# --- plot one trial 

i <- 1 #trial to show

pdf("Figure_1run.pdf", width = 5, height = 4)
par(las = 1)
plot(rev(time) - max(time), weak0204.truetraj[,i], type = "l", xlim = c(-.1, 0),
 bty = "n", xlab = "Time (coalescent units, past to left)", 
ylab = "Population-average polygenic score")
neut.est <- weak0204.truetraj[,i] + weak0204.err.array[,1,i]
lines(rev(time) - max(time), neut.est, col = pal[1])

mom.est <- weak0204.truetraj[,i] + weak0204.err.array[,2,i]
lines(rev(time) - max(time), mom.est, col = pal[2])

lines(rev(time) - max(time), weak0204.truetraj[,i], col = "black")
lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
legend("topleft", bty = "n", lty = 1, col = c("black", pal[1:2]), legend = c("true", "proportion-of-lineages", "lineages remaining"))
dev.off()

#############################################################
#--- multi-panel plot of bias and MSE --- both true and rent+

xl <- c(-.25, 0)
yl.b <- c(-1,2)
yl.mse <- c(0,4)
wid <- 3.5
hei <- 2.5

pdf("Figure_biasMSE.pdf", width = wid*2, height = hei)
par(mfrow = c(1,2), mar = c(4,4,2,1), las = 1, mgp = c(2,.7,0), cex = .85)


plot(rev(time) - max(time), weak0204.avg.err.mat[,1], type = "l", xlim = xl, ylim = yl.b,
 bty = "n", xlab = "Time", 
ylab = "Bias", col = pal[1])
lines(rev(time) - max(time), weak0204.avg.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), weak0204.avg.err.mat[,3] , col = pal[3])
lines(rev(time) - max(time), weak0204.avg.err.mat[,5] , col = pal[4])

lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,1] , col = pal[1], lty = 2)
lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,2] , col = pal[2], lty = 2)
lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,3] , col = pal[3], lty = 2)

lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)

mtext("A", side = 3, line = .5, at = min(xl))


plot(rev(time) - max(time), weak0204.avg.sq.err.mat[,1], type = "l", xlim = xl, ylim = yl.mse,
 bty = "n", xlab = "Time", 
ylab = "Mean squared error", col = pal[1])
lines(rev(time) - max(time), weak0204.avg.sq.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), weak0204.avg.sq.err.mat[,3] , col = pal[3])
lines(rev(time) - max(time), weak0204.avg.sq.err.mat[,5] , col = pal[4])

lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,1] , col = pal[1], lty = 2)
lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,2] , col = pal[2], lty = 2)
lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,3] , col = pal[3], lty = 2)

lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
mtext("F", side = 3, line = .5, at = min(xl))
legend("topleft", bty = "n", lty = 1, col = c(pal[c(1,3,2,4)]), legend = c("proportion-of-lineages", "waiting time", "lineages remaining", "straight line to ancestral state"))
mtext("B", side = 3, line = .5, at = min(xl))

dev.off()



#############################################################
# multi-panel plot of  confidence-interval coverage.

xl <- c(-.25, 0)
yl.conf <- c(.6,1)

wid <- 3.5
hei <- 2.5

pdf("Figure_CI.pdf", width = wid*2, height = hei)
par(mfrow = c(1,2), mar = c(4,4,2,1), las = 1, mgp = c(2.5,.7,0), cex = .8)


plot(rev(time) - max(time), weak0204.covered.mat[,1], type = "l", xlim = xl, ylim = yl.conf,
 bty = "n", xlab = "Time", 
ylab = "Coverage (true trees)", col = pal[1])
lines(rev(time) - max(time), weak0204.covered.mat[,2] , col = pal[2])
lines(rev(time) - max(time), weak0204.covered.mat[,3] , col = pal[3])

lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)

mtext("E", side = 3, line = .5, at = min(xl))

plot(rev(time) - max(time), weak0204.rent.covered.mat[,1], type = "l", xlim = xl, ylim = yl.conf,
 bty = "n", xlab = "Time", 
ylab = "Coverage (RENT+)", col = pal[1], lty = 2)
lines(rev(time) - max(time), weak0204.rent.covered.mat[,2] , col = pal[2], lty = 2)
lines(rev(time) - max(time), weak0204.rent.covered.mat[,3] , col = pal[3], lty = 2)

lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)

mtext("F", side = 3, line = .5, at = min(xl))

dev.off()






