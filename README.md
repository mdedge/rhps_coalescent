# rhps_coalescent
This repository contains code associated with the paper "Reconstructing the history of polygenic scores using coalescent trees" by M.D. Edge and G. Coop.

The functions that encode the operations described in the main text are defined in helper_functions_coal_sel.R. This is the best place to look if you want to see how we implemented the main procedures, including allele-frequency trajectory simulation, handing of coalescent trees, estimation, and hypothesis testing. Most of the functions in the file are for handling phylo objects as used in package ape, and many of the functions depend on package ape.

msseldir contains files to run the version of mssel we used for simulations.

RentPlus.jar is a java jar file used to run RENT+. This version of RENT+ is slightly altered from the one at https://github.com/SajadMirzaei/RentPlus.

The height_analyses directory contains scripts used to run the height analyses appearing in figure 7, including data files with the SNP names and effect sizes. I ran the analyses using a copy of the 1000 Genomes data downloaded to my local machine. To reproduce the analyses, you would need to alter the script to access either a copy of the data on your machine or stored remotely.

The maintext_sims_rent_061318 directory contains all the R scripts used to simulate allele frequency trajectories, call mssel to generate trees and data, call RENT+ to infer trees, and use the procedures in the paper to estimate polygenic score trajectories and test them for selection. These scripts consist primarily of calls to mssel, rent+, and the functions defined in helper_functions_coal_sel.R. The code in this directory is rather sparsely documented and is uploaded mostly for completeness.
