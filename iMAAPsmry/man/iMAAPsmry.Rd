\name{iMAAPssmry}
\alias{iMAAPssmry}
\title{Time summary of iMAAPs}
\description{
	Output the summary spectrum and estimate the admxiture time interval from the results of iMAAPs.
}

\usage{
	iMAAPssmry(ad_inp, rawld_inp, num_chr, denoise_cut, method="median")
}

\arguments{
	\item{ad_inp}{is a matrix of fitting coeffients}
	\item{rawld_inp}{is a matrix of Weighted LD decays}
	\item{denoise_cut}{is a number between 0 and 1 that defines the noise cut threshold}
	\item{num_chr}{is a number of the chromosome used for Jacknife process}
	\item{method}{is a function used to estimate the average of the admixture signal, it has two choices: median and mean}
}

\examples{
library(iMAAPssmry)
data(ad_inp)#or#ad_inp<-read.csv(file="../example/ADM.ad",header=T,sep="\t");
data(rawld_inp)#or#rwald_inp<-read.csv(file="../example/ADM.rawld",header=T,sep="\t");
x<-iMAAPssmry(
        ad_inp=ad_inp,
        rawld_inp=rawld_inp,
        denoise_cut=0.99, num_chr=10);
#summary spectrum
x$sum_spectrum
#estimated admixture time
x$time
}
