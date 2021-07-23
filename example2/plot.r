#library(PIMAsummary, lib="~/Rlibs/");
source("/picb/humpopg-sim/mathew_prj/MALDmef_pack/PIMA-pkg/R/PIMA_time.r")
PIMA_time(
        ad_inp="ADM.ad",
        rawld_inp="ADM.rawld",
        denoise_cut=0.99, num_chr=10)



