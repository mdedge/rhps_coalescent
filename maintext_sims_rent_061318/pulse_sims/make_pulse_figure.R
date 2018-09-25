
#complete phenotype simulations by looping through parameter values and calling 
#pheno_sim_1iter.R for each set of parameters.

N <- 10000
herit <- 1
out_dir <- "out_061318/"

n_chromss <- c(20, 50, 100, 200, 500, 1000)
t.offs <- c(0, .01, .02, .03, .04, .05, .06, .07, .08, .09)
ts <- t.offs + .005
phen_nums <- 1:100

helper_fn <- "../../helper_functions_coal_sel.R"
source(helper_fn) #read in helper functions and load ape package

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

pows <- aggregate(all_qxtest_mat[,4], all_pars[,c(5,1)], function(x){mean(x < .05)})

pdf("Figure_pulse.pdf", width = 6, height = 4.5)
pal <- c("purple", "pink", "blue", "orange", "brown", "grey")
plot(pows[pows[,2] == 1000,1], pows[pows[,2] == 1000,3], ylim = c(0,1), type = "o", col = pal[1],
xlim = c(0,.1), ylab = "Power", xlab = "End of selection interval (coalescent units before present)", bty = "n", las = 1, pch = 20)
points(pows[pows[,2] == 500,1], pows[pows[,2] == 500,3], col = pal[2], type = "o", pch = 20)
points(pows[pows[,2] == 200,1], pows[pows[,2] == 200,3], col = pal[3], type = "o", pch = 20)
points(pows[pows[,2] == 100,1], pows[pows[,2] == 100,3], col = pal[4], type = "o", pch = 20)
points(pows[pows[,2] == 50,1], pows[pows[,2] == 50,3], col = pal[5], type = "o", pch = 20)
points(pows[pows[,2] == 20,1], pows[pows[,2] == 20,3], col = pal[6], type = "o", pch = 20)
legend("topright", bty = "n", lty = 1, col = pal, legend = c("n = 1000","n = 500","n = 200","n = 100","n = 50","n = 20"))
dev.off()


