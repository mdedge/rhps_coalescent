#1/30/18
#This is super slow, with lots of for() loops. Worth optimizing if we use more traits.
#Goal: Take biobank effect sizes, find the highest Bayes factor in each of Pickrell's ~1700
#regions. End up with a dataframe of basic statistics for these loci.

#Downloaded with
#wget https://www.dropbox.com/s/sbfgb6qd5i4cxku/50.assoc.tsv.gz?dl=0 -O 50.assoc.tsv.gz
#Commands for downloading each phenotype found here https://docs.google.com/spreadsheets/d/1b3oGI2lUt57BcuHttWaZotQcI0-mBRPyZihz87Ms_No/edit#gid=1209628142
#and referenced here: https://sites.google.com/broadinstitute.org/ukbbgwasresults/home?authuser=0

#Bayes factors per site computed as in Pickrell (2014) AJHG Eq. 4
#BF = sqrt(1-r)/exp((z^2) * r / 2)
#where r = w/(v+w) and w is a prior variance (set to 0.1) and v is var estimate (se^2)
#and z is a standard z stat (estimate over standard error)

#Joe Pickrell's breakpoints: https://bitbucket.org/nygcresearch/ldetect-data

library(data.table)
dat <- fread("standing_height_UKB.assoc.tsv")
breaks <- fread("breakpoints_Pickrell_EUR.bed", head = TRUE)


r <- 0.1/(0.1 + dat$se^2)
bf <- sqrt(1 - r)/exp(-(dat$tstat^2) * r / 2)

splits <- strsplit(dat$variant, ":")

chr.dat <- numeric()
pos.dat <- numeric()
for(i in 1:length(splits)){
	chr.dat[i] <- as.numeric(splits[[i]][1])
	pos.dat[i] <- as.numeric(splits[[i]][2])
}

break.chr <- as.numeric(gsub("chr", "", breaks$chr))

segment <- numeric(length(chr.dat))
for(i in 1:length(break.chr)){
	reg.chr <- break.chr[i]
	segment[chr.dat == reg.chr & pos.dat >= breaks$start[i] & pos.dat < breaks$stop[i]] <- i
}

write.table(max.bf.loci, file = "biobank_height_maxBF_PickrellWindows.txt", quote = FALSE, row.names = FALSE)

max.bf.loci <- dat[segment == 1,][which.max(bf[segment == 1]),]
for(i in 2:nrow(breaks)){
	max.bf.loci <- rbind(max.bf.loci, dat[segment == i,][which.max(bf[segment == i]),])	
}



