#Doc Edge, 7/5/18
#Extract data from simulations and make figures.
#to be run in maintext_sims_rent

load("neutral/analyzed_true_trees_072718.RData")

neut.truetraj <- mat.true.phentrajs
neut.err.array <- err.array
neut.err.std.array <- err.std.array
neut.qxtest_mat <- qxtest_mat


load("neutral/analyzed_rent_trees_072718.RData")

neut.rent.err.array <- err.array.rent
neut.rent.err.std.array <- err.std.array.rent
neut.rent.qxtest_mat <- qxtest_mat_rent

load("weaksel_02/analyzed_true_trees_072718.RData")

weak02.truetraj <- mat.true.phentrajs
weak02.err.array <- err.array
weak02.err.std.array <- err.std.array
weak02.qxtest_mat <- qxtest_mat

load("weaksel_02/analyzed_rent_trees_072718.RData")

weak02.rent.err.array <- err.array.rent
weak02.rent.err.std.array <- err.std.array.rent
weak02.rent.qxtest_mat <- qxtest_mat_rent

load("weaksel_0204/analyzed_true_trees_072718.RData")

weak0204.truetraj <- mat.true.phentrajs
weak0204.err.array <- err.array
weak0204.err.std.array <- err.std.array
weak0204.qxtest_mat <- qxtest_mat

load("weaksel_0204/analyzed_rent_trees_072718.RData")

weak0204.rent.err.array <- err.array.rent
weak0204.rent.err.std.array <- err.std.array.rent
weak0204.rent.qxtest_mat <- qxtest_mat_rent



neut.avg.err.mat <- apply(neut.err.array, 1:2, mean, na.rm = TRUE)
neut.avg.sq.err.mat <- apply(neut.err.array^2, 1:2, mean, na.rm = TRUE)

neut.rent.avg.err.mat <- apply(neut.rent.err.array, 1:2, mean, na.rm = TRUE)
neut.rent.avg.sq.err.mat <- apply(neut.rent.err.array^2, 1:2, mean, na.rm = TRUE)

pal <- c('#1b9e77','#d95f02','#7570b3','#e7298a')

plot(time, neut.avg.sq.err.mat[,1], type = "l", ylim = c(0,2), xlim = c(0,1/2))
lines(time, neut.avg.sq.err.mat[,2], col = "orange")
lines(time, neut.avg.sq.err.mat[,3], col = "blue")
#lines(time, neut.avg.sq.err.mat[,4], col = "red")
lines(time, neut.avg.sq.err.mat[,5], col = "green")

plot(time, neut.avg.err.mat[,1], type = "l", ylim = c(-1/2,1/2), xlim = c(0,1/2))
lines(time, neut.avg.err.mat[,2], col = "orange")
lines(time, neut.avg.err.mat[,3], col = "blue")
#lines(time, neut.avg.err.mat[,4], col = "red")
lines(time, neut.avg.err.mat[,5], col = "green")



plot(time, neut.rent.avg.sq.err.mat[,1], type = "l", ylim = c(0,2), xlim = c(0,1/2))
lines(time, neut.rent.avg.sq.err.mat[,2], col = "orange")
lines(time, neut.rent.avg.sq.err.mat[,3], col = "blue")
#lines(time, neut.rent.avg.sq.err.mat[,4], col = "red")
lines(time, neut.avg.sq.err.mat[,5], col = "green")

plot(time, neut.rent.avg.err.mat[,1], type = "l", ylim = c(-1/2,1/2), xlim = c(0,1/2))
lines(time, neut.rent.avg.err.mat[,2], col = "orange")
lines(time, neut.rent.avg.err.mat[,3], col = "blue")
#lines(time, neut.rent.avg.err.mat[,4], col = "red")
lines(time, neut.avg.err.mat[,5], col = "green")


weak02.avg.err.mat <- apply(weak02.err.array, 1:2, mean, na.rm = TRUE)
weak02.avg.sq.err.mat <- apply(weak02.err.array^2, 1:2, mean, na.rm = TRUE)

weak02.rent.avg.err.mat <- apply(weak02.rent.err.array, 1:2, mean, na.rm = TRUE)
weak02.rent.avg.sq.err.mat <- apply(weak02.rent.err.array^2, 1:2, mean, na.rm = TRUE)


plot(time, weak02.avg.sq.err.mat[,1], type = "l", xlim = c(0,1/2))
lines(time, weak02.avg.sq.err.mat[,2], col = "orange")
lines(time, weak02.avg.sq.err.mat[,3], col = "blue")
#lines(time, weak02.avg.sq.err.mat[,4], col = "red")
lines(time, weak02.avg.sq.err.mat[,5], col = "green")

plot(time, weak02.avg.err.mat[,1], type = "l", xlim = c(0,1/2))
lines(time, weak02.avg.err.mat[,2], col = "orange")
lines(time, weak02.avg.err.mat[,3], col = "blue")
#lines(time, weak02.avg.err.mat[,4], col = "red")
lines(time, weak02.avg.err.mat[,5], col = "green")



plot(time, weak02.rent.avg.sq.err.mat[,1], type = "l", xlim = c(0,1/2))
lines(time, weak02.rent.avg.sq.err.mat[,2], col = "orange")
lines(time, weak02.rent.avg.sq.err.mat[,3], col = "blue")
#lines(time, weak02.avg.sq.err.mat[,4], col = "red")
lines(time, weak02.avg.sq.err.mat[,5], col = "green")

plot(time, weak02.rent.avg.err.mat[,1], type = "l", xlim = c(0,1/2))
lines(time, weak02.rent.avg.err.mat[,2], col = "orange")
lines(time, weak02.rent.avg.err.mat[,3], col = "blue")
#lines(time, weak02.avg.err.mat[,4], col = "red")
lines(time, weak02.avg.err.mat[,5], col = "green")



weak0204.avg.err.mat <- apply(weak0204.err.array, 1:2, mean, na.rm = TRUE)
weak0204.avg.sq.err.mat <- apply(weak0204.err.array^2, 1:2, mean, na.rm = TRUE)

weak0204.rent.avg.err.mat <- apply(weak0204.rent.err.array, 1:2, mean, na.rm = TRUE)
weak0204.rent.avg.sq.err.mat <- apply(weak0204.rent.err.array^2, 1:2, mean, na.rm = TRUE)


plot(time, weak0204.avg.sq.err.mat[,1], type = "l", xlim = c(0,1/2))
lines(time, weak0204.avg.sq.err.mat[,2], col = "orange")
lines(time, weak0204.avg.sq.err.mat[,3], col = "blue")
#lines(time, weak0204.avg.sq.err.mat[,4], col = "red")
lines(time, weak0204.avg.sq.err.mat[,5], col = "green")

plot(time, weak0204.avg.err.mat[,1], type = "l", xlim = c(0,1/2))
lines(time, weak0204.avg.err.mat[,2], col = "orange")
lines(time, weak0204.avg.err.mat[,3], col = "blue")
#lines(time, weak0204.avg.err.mat[,4], col = "red")
lines(time, weak0204.avg.err.mat[,5], col = "green")



plot(time, weak0204.rent.avg.sq.err.mat[,1], type = "l", xlim = c(0,1/2))
lines(time, weak0204.rent.avg.sq.err.mat[,2], col = "orange")
lines(time, weak0204.rent.avg.sq.err.mat[,3], col = "blue")
#lines(time, weak0204.avg.sq.err.mat[,4], col = "red")
lines(time, weak0204.avg.sq.err.mat[,5], col = "green")

plot(time, weak0204.rent.avg.err.mat[,1], type = "l", xlim = c(0,1/2))
lines(time, weak0204.rent.avg.err.mat[,2], col = "orange")
lines(time, weak0204.rent.avg.err.mat[,3], col = "blue")
#lines(time, weak0204.avg.err.mat[,4], col = "red")
lines(time, weak0204.avg.err.mat[,5], col = "green")











#Confidence intervals
prop.covered <- function(vec, q = 1.96, ...){
	mean(abs(vec) <= q, ...)
}

neut.covered.mat <- apply(neut.err.std.array, 1:2, prop.covered, na.rm = TRUE)

plot(time, neut.covered.mat[,1], type = "l", xlim = c(0,1/2), ylim = c(.8,1))
lines(time, neut.covered.mat[,2], col = "orange")
lines(time, neut.covered.mat[,3], col = "blue")
lines(c(0,10), c(.95,.95), col = "grey", lty = 2)

neut.rent.covered.mat <- apply(neut.rent.err.std.array, 1:2, prop.covered, na.rm = TRUE)

plot(time, neut.rent.covered.mat[,1], type = "l", xlim = c(0,1/2), ylim = c(.8,1))
lines(time, neut.rent.covered.mat[,2], col = "orange")
lines(time, neut.rent.covered.mat[,3], col = "blue")
lines(c(0,10), c(.95,.95), col = "grey", lty = 2)





weak02.covered.mat <- apply(weak02.err.std.array, 1:2, prop.covered, na.rm = TRUE)

plot(time, weak02.covered.mat[,1], type = "l", xlim = c(0,1/2), ylim = c(.8,1))
lines(time, weak02.covered.mat[,2], col = "orange")
lines(time, weak02.covered.mat[,3], col = "blue")
lines(c(0,10), c(.95,.95), col = "grey", lty = 2)

weak02.rent.covered.mat <- apply(weak02.rent.err.std.array, 1:2, prop.covered, na.rm = TRUE)

plot(time, weak02.rent.covered.mat[,1], type = "l", xlim = c(0,1/2), ylim = c(.8,1))
lines(time, weak02.rent.covered.mat[,2], col = "orange")
lines(time, weak02.rent.covered.mat[,3], col = "blue")
lines(c(0,10), c(.95,.95), col = "grey", lty = 2)



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



#tests
apply(neut.qxtest_mat[,c(3,4,7,8,11,12)], 2, function(x){mean(x<.05)})
apply(neut.rent.qxtest_mat[,c(3,4,7,8,11,12)], 2, function(x){mean(x<.05)})

apply(weak02.qxtest_mat[,c(3,4,7,8,11,12)], 2, function(x){mean(x<.05)})
apply(weak02.rent.qxtest_mat[,c(3,4,7,8,11,12)], 2, function(x){mean(x<.05)})

apply(weak0204.qxtest_mat[,c(3,4,7,8,11,12)], 2, function(x){mean(x<.05)})
apply(weak0204.rent.qxtest_mat[,c(3,4,7,8,11,12)], 2, function(x){mean(x<.05)})


#There was a mistake in the Qx_test code when I first ran this; it ought
#to be a 10-df chisq and not 9-df. The following gives the p values when 10-df chisq
#is used.
apply(neut.qxtest_mat[,c(1,5,9)], 2, function(x, df){mean(pchisq(x, df, lower.tail = FALSE) < .05)}, df = 10)
apply(neut.rent.qxtest_mat[,c(1,5,9)], 2, function(x, df){mean(pchisq(x, df, lower.tail = FALSE) < .05)}, df = 10)

apply(weak02.qxtest_mat[,c(1,5,9)], 2, function(x, df){mean(pchisq(x, df, lower.tail = FALSE) < .05)}, df = 10)
#apply(weak02.rent.qxtest_mat[,c(1,5,9)], 2, function(x, df){mean(pchisq(x, df, lower.tail = FALSE) < .05)}, df = 10)

apply(weak0204.qxtest_mat[,c(1,5,9)], 2, function(x, df){mean(pchisq(x, df, lower.tail = FALSE) < .05)}, df = 10)
#apply(weak0204.rent.qxtest_mat[,c(1,5,9)], 2, function(x, df){mean(pchisq(x, df, lower.tail = FALSE) < .05)}, df = 10)



pal <- c('#1b9e77','#d95f02','#7570b3','#e7298a')
#############################################################
#Figure 3 --- pull out one trial with weak selection
#from .04 to .02 and show result

i <- 3 #trial to show

pdf("Figure3.pdf", width = 5, height = 4)
par(las = 1)
plot(rev(time) - max(time), weak0204.truetraj[,i], type = "l", xlim = c(-.1, 0), ylim = c(-2.5,1),
 bty = "n", xlab = "Time (coalescent units, past to left)", 
ylab = "Population-average polygenic score")
neut.est <- weak0204.truetraj[,i] + weak0204.err.array[,1,i]
#neut.se <- weak0204.err.array[,1,i] / weak0204.err.std.array[,1,i]
#up.neut <- neut.est + 1.96*neut.se
#dn.neut <- neut.est - 1.96*neut.se
lines(rev(time) - max(time), neut.est, col = pal[1])
#lines(rev(time) - max(time), up.neut, col = pal[1])
#lines(rev(time) - max(time), dn.neut, col = pal[1])

#wt.est <- weak0204.truetraj[,i] + weak0204.err.array[,3,i]
#wt.se <- weak0204.err.array[,3,i] / weak0204.err.std.array[,3,i]
#up.wt <- wt.est + 1.96*wt.se
#dn.wt <- wt.est - 1.96*wt.se
#lines(rev(time) - max(time), mom.est, col = pal[3])
#lines(rev(time) - max(time), up.mom, col = pal[3])
#lines(rev(time) - max(time), dn.mom, col = pal[3])

mom.est <- weak0204.truetraj[,i] + weak0204.err.array[,2,i]
#mom.se <- weak0204.err.array[,2,i] / weak0204.err.std.array[,2,i]
#up.mom <- mom.est + 1.96*mom.se
#dn.mom <- mom.est - 1.96*mom.se
lines(rev(time) - max(time), mom.est, col = pal[2])
#lines(rev(time) - max(time), up.mom, col = pal[2])
#lines(rev(time) - max(time), dn.mom, col = pal[2])

#lines(rev(time) - max(time), weak0204.truetraj[,i] + weak0204.err.array[,3,i], col = pal[3])
lines(rev(time) - max(time), weak0204.truetraj[,i], col = "black")
#lines(rev(time) - max(time), weak0204.truetraj[,i] + weak0204.err.array[,5,i], col = pal[4])
lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
legend("topleft", bty = "n", lty = 1, col = c("black", pal[1:2]), legend = c("true", "proportion-of-lineages", "lineages remaining"))
dev.off()

#############################################################
#Figure 4 --- multi-panel plot of bias and MSE --- both true and rent+

xl <- c(-.25, 0)
yl.b <- c(-1,2)
yl.mse <- c(0,4)
wid <- 3.5
hei <- 2.5

pdf("Figure4.pdf", width = wid*2, height = hei*3)
par(mfrow = c(3,2), mar = c(4,4,2,1), las = 1, mgp = c(2,.7,0), cex = .85)

plot(rev(time) - max(time), neut.avg.err.mat[,1], type = "l", xlim = xl, ylim = yl.b,
 bty = "n", xlab = "Time (coalescent units, past to left)", 
ylab = "Bias", col = pal[1])
lines(rev(time) - max(time), neut.avg.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), neut.avg.err.mat[,3] , col = pal[3])
lines(rev(time) - max(time), neut.avg.err.mat[,5] , col = pal[4])

lines(rev(time) - max(time), neut.rent.avg.err.mat[,1] , col = pal[1], lty = 2)
lines(rev(time) - max(time), neut.rent.avg.err.mat[,2] , col = pal[2], lty = 2)
lines(rev(time) - max(time), neut.rent.avg.err.mat[,3] , col = pal[3], lty = 2)

legend("topleft", bty = "n", lty = 1, col = c(pal[c(1,3,2,4)]), legend = c("proportion-of-lineages", "waiting time", "lineages remaining", "straight line to ancestral state"))
mtext("A", side = 3, line = .5, at = min(xl))


plot(rev(time) - max(time), neut.avg.sq.err.mat[,1], type = "l", xlim = xl, ylim = yl.mse,
 bty = "n", xlab = "Time", 
ylab = "Mean squared error", col = pal[1])
lines(rev(time) - max(time), neut.avg.sq.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), neut.avg.sq.err.mat[,3] , col = pal[3])
lines(rev(time) - max(time), neut.avg.sq.err.mat[,5] , col = pal[4])

lines(rev(time) - max(time), neut.rent.avg.sq.err.mat[,1] , col = pal[1], lty = 2)
lines(rev(time) - max(time), neut.rent.avg.sq.err.mat[,2] , col = pal[2], lty = 2)
lines(rev(time) - max(time), neut.rent.avg.sq.err.mat[,3] , col = pal[3], lty = 2)

mtext("B", side = 3, line = .5, at = min(xl))


plot(rev(time) - max(time), weak02.avg.err.mat[,1], type = "l", xlim = xl, ylim = yl.b,
 bty = "n", xlab = "Time", 
ylab = "Bias", col = pal[1])
lines(rev(time) - max(time), weak02.avg.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), weak02.avg.err.mat[,3] , col = pal[3])
lines(rev(time) - max(time), weak02.avg.err.mat[,5] , col = pal[4])

lines(rev(time) - max(time), weak02.rent.avg.err.mat[,1] , col = pal[1], lty = 2)
lines(rev(time) - max(time), weak02.rent.avg.err.mat[,2] , col = pal[2], lty = 2)
lines(rev(time) - max(time), weak02.rent.avg.err.mat[,3] , col = pal[3], lty = 2)

lines(-c(.0, .0), c(-100, 100), col = "gray", lty = 2)
lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)

mtext("C", side = 3, line = .5, at = min(xl))


plot(rev(time) - max(time), weak02.avg.sq.err.mat[,1], type = "l", xlim = xl, ylim = yl.mse,
 bty = "n", xlab = "Time", 
ylab = "Mean squared error", col = pal[1])
lines(rev(time) - max(time), weak02.avg.sq.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), weak02.avg.sq.err.mat[,3] , col = pal[3])
lines(rev(time) - max(time), weak02.avg.sq.err.mat[,5] , col = pal[4])

lines(rev(time) - max(time), weak02.rent.avg.sq.err.mat[,1] , col = pal[1], lty = 2)
lines(rev(time) - max(time), weak02.rent.avg.sq.err.mat[,2] , col = pal[2], lty = 2)
lines(rev(time) - max(time), weak02.rent.avg.sq.err.mat[,3] , col = pal[3], lty = 2)

lines(-c(.0, .0), c(-100, 100), col = "gray", lty = 2)
lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
mtext("D", side = 3, line = .5, at = min(xl))



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

mtext("E", side = 3, line = .5, at = min(xl))


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

dev.off()



#############################################################
#Figure 5 --- multi-panel plot of  confidence-interval coverage.

xl <- c(-.25, 0)
yl.conf <- c(.6,1)

wid <- 3.5
hei <- 2.5

pdf("Figure5.pdf", width = wid*2, height = hei*3)
par(mfrow = c(3,2), mar = c(4,4,2,1), las = 1, mgp = c(2.5,.7,0), cex = .8)

plot(rev(time) - max(time), neut.covered.mat[,1], type = "l", xlim = xl, ylim = yl.conf,
 bty = "n", xlab = "Time (coalescent units, past to left)", 
ylab = "Coverage (true trees)", col = pal[1])
lines(rev(time) - max(time), neut.covered.mat[,2] , col = pal[2])
lines(rev(time) - max(time), neut.covered.mat[,3] , col = pal[3])

legend("bottomleft", bty = "n", lty = 1, col = c(pal[c(1,3,2)]), legend = c("proportion-of-lineages", "waiting time", "lineages remaining"))
mtext("A", side = 3, line = .5, at = min(xl))

plot(rev(time) - max(time), neut.rent.covered.mat[,1], type = "l", xlim = xl, ylim = yl.conf,
 bty = "n", xlab = "Time", 
ylab = "Coverage (RENT+)", col = pal[1], lty = 2)
lines(rev(time) - max(time), neut.rent.covered.mat[,2] , col = pal[2], lty = 2)
lines(rev(time) - max(time), neut.rent.covered.mat[,3] , col = pal[3], lty = 2)

mtext("B", side = 3, line = .5, at = min(xl))

plot(rev(time) - max(time), weak02.covered.mat[,1], type = "l", xlim = xl, ylim = yl.conf,
 bty = "n", xlab = "Time", 
ylab = "Coverage (true trees)", col = pal[1])
lines(rev(time) - max(time), weak02.covered.mat[,2] , col = pal[2])
lines(rev(time) - max(time), weak02.covered.mat[,3] , col = pal[3])

lines(-c(.0, .0), c(-100, 100), col = "gray", lty = 2)
lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)

mtext("C", side = 3, line = .5, at = min(xl))

plot(rev(time) - max(time), weak02.rent.covered.mat[,1], type = "l", xlim = xl, ylim = yl.conf,
 bty = "n", xlab = "Time", 
ylab = "Coverage (RENT+)", col = pal[1], lty = 2)
lines(rev(time) - max(time), weak02.rent.covered.mat[,2] , col = pal[2], lty = 2)
lines(rev(time) - max(time), weak02.rent.covered.mat[,3] , col = pal[3], lty = 2)

lines(-c(.0, .0), c(-100, 100), col = "gray", lty = 2)
lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)

mtext("D", side = 3, line = .5, at = min(xl))

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






