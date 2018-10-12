#Doc Edge, 7/5/18
#Extract data from simulations and make figures.
#to be run in maintext_sims_rent

load("analyzed_true_trees_100118.RData")

saw.truetraj <- mat.true.phentrajs
saw.err.array <- err.array
saw.err.std.array <- err.std.array
saw.qxtest_mat <- qxtest_mat



saw.avg.err.mat <- apply(saw.err.array, 1:2, mean, na.rm = TRUE)
saw.avg.sq.err.mat <- apply(saw.err.array^2, 1:2, mean, na.rm = TRUE)

pal <- c('#1b9e77','#d95f02','#7570b3','#e7298a')

plot(time, saw.avg.sq.err.mat[,1], type = "l", ylim = c(0,2), xlim = c(0,1/2))
lines(time, saw.avg.sq.err.mat[,2], col = "orange")
lines(time, saw.avg.sq.err.mat[,3], col = "blue")
#lines(time, neut.avg.sq.err.mat[,4], col = "red")
lines(time, saw.avg.sq.err.mat[,5], col = "green")




#Confidence intervals
prop.covered <- function(vec, q = 1.96, ...){
	mean(abs(vec) <= q, ...)
}

saw.covered.mat <- apply(saw.err.std.array, 1:2, prop.covered, na.rm = TRUE)

plot(time, saw.covered.mat[,1], type = "l", xlim = c(0,1/2), ylim = c(.8,1))
lines(time, saw.covered.mat[,2], col = "orange")
lines(time, saw.covered.mat[,3], col = "blue")
lines(c(0,10), c(.95,.95), col = "grey", lty = 2)


#tests
apply(saw.qxtest_mat[,c(3,4,7,8,11,12)], 2, function(x){mean(x<.05)})


#There was a mistake in the Qx_test code when I first ran this; it ought
#to be a 10-df chisq and not 9-df. The following gives the p values when 10-df chisq
#is used.
apply(saw.qxtest_mat[,c(1,5,9)], 2, function(x, df){mean(pchisq(x, df, lower.tail = FALSE) < .05)}, df = 10)



pal <- c('#1b9e77','#d95f02','#7570b3','#e7298a')
#############################################################
#Figure 3 --- pull out one trial with weak selection
#from .04 to .02 and show result

i <- 1 #trial to show

pdf("Figure1trial.pdf", width = 5, height = 4)
par(las = 1)
plot(rev(time) - max(time), saw.truetraj[,i], type = "l", xlim = c(-.1, 0), ylim = c(-1,3),
 bty = "n", xlab = "Time (coalescent units, past to left)", 
ylab = "Population-average polygenic score")
neut.est <- saw.truetraj[,i] + saw.err.array[,1,i]

lines(rev(time) - max(time), neut.est, col = pal[1])


mom.est <- saw.truetraj[,i] + saw.err.array[,2,i]
lines(rev(time) - max(time), mom.est, col = pal[2])

#lines(rev(time) - max(time), weak0204.truetraj[,i] + weak0204.err.array[,3,i], col = pal[3])
lines(rev(time) - max(time), saw.truetraj[,i], col = "black")
#lines(rev(time) - max(time), weak0204.truetraj[,i] + weak0204.err.array[,5,i], col = pal[4])
lines(-c(0, 0), c(-100, 100), col = "gray", lty = 2)
lines(-c(.01, .01), c(-100, 100), col = "gray", lty = 2)
lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
legend("topleft", bty = "n", lty = 1, col = c("black", pal[1:2]), legend = c("true", "proportion-of-lineages", "lineages remaining"))
dev.off()

#############################################################
#Figure 4 --- multi-panel plot of bias and MSE --- both true and rent+

xl <- c(-.25, 0)
yl.b <- c(-1,2)
yl.mse <- c(0,4)
wid <- 3.5
hei <- 2.5

pdf("Figurebias.pdf", width = wid*2, height = hei*1)
par(mfrow = c(1,2), mar = c(4,4,2,1), las = 1, mgp = c(2,.7,0), cex = .85)

plot(rev(time) - max(time), saw.avg.err.mat[,1], type = "l", xlim = xl, ylim = yl.b,
 bty = "n", xlab = "Time (coalescent units, past to left)", 
ylab = "Bias", col = pal[1])
lines(rev(time) - max(time), saw.avg.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), saw.avg.err.mat[,3] , col = pal[3])
lines(rev(time) - max(time), saw.avg.err.mat[,5] , col = pal[4])

legend("topleft", bty = "n", lty = 1, col = c(pal[c(1,3,2,4)]), legend = c("proportion-of-lineages", "waiting time", "lineages remaining", "straight line to ancestral state"))
mtext("A", side = 3, line = .5, at = min(xl))


plot(rev(time) - max(time), saw.avg.sq.err.mat[,1], type = "l", xlim = xl, ylim = yl.mse,
 bty = "n", xlab = "Time", 
ylab = "Mean squared error", col = pal[1])
lines(rev(time) - max(time), saw.avg.sq.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), saw.avg.sq.err.mat[,3] , col = pal[3])
lines(rev(time) - max(time), saw.avg.sq.err.mat[,5] , col = pal[4])

mtext("B", side = 3, line = .5, at = min(xl))

dev.off()






