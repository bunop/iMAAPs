# iMAAPs
To infer multiple-wave admixture by fitting ALD using a p-spectrum

iMAAPs is a powerful tool to estimate multiple-wave population admixed time, which is currently designed to infer the two-way, multiple-wave admixture based on admixture induced LD. This software can deal with genotype data, haplotype data and the data re-coded according to admixture ancestries.

Related Publication: Heredity. 118(5):503-510. Link:https://www.nature.com/articles/hdy20175


##########################################################################################
iMAAPs V1.0.0
##########################################################################################
Infer Multiple-wave Admixture by fitting Admixture LD using Polynomial spectrum

The iMAAPs is a powerful tool to estimate multiple-wave population admixed time, which is a software currently designed to infer the two-way, multiple-wave admixture based on admixture induced LD.

For any technical questions, please contact 
{yingzhou&yuankai}@picb.ac.cn
##########################################################################################
1. INSTALL
iMAAPs is currently developed for Linux environments. The source code is availble at website:
http://www.picb.ac.cn/PGG/resource.php
To compile iMAAPs, GSL, FFTW and OpenMP should be installed properly. Please see details at

GSL
http://www.gnu.org/software/gsl/
FFTW
http://www.fftw.org/
OpenMP
http://www.openmp.org/

Then go to the main directory of iMAAPs, "make" and then you get the executable file iMAAPs.
2, Quick start
The basic command follows as:
./iMAAPs -p [parameter file]

A simple example:
./iMAAPs adm.par > adm.log

The parameter file contains the inputs and parameters, and using “./iMAAPs -P” or “./iMAAPs -help para” to get a brief description of parameter file.

3. Inputs
a> Input lists
The file(*.list) contains all the data input of genotype files and SNP ananotation. These files are seperated into different chromosomes. 
	1st column: Chromosome ID
	2nd column: Genotype files
	3rd contains: SNP files

b> Genotype file
The genotype is in EIGENSTRAT format. Each columns represents a sample individual and each row represents position for that SNP. We use 0, 1 and 2 to code the binary data, where 9 is missing. Markers with more than two states should be deleted from analyis. Any characters different from 0,1 and 2 are dealled as missing in our programe without any warning.


c> SNP file
Each file contains the SNP information for one chromosome.
	1st column: Chromosome.  Use X for X chromosome.  	
	2nd column: SNP name
  	3rd column: Genetic position (in Morgans)
  	4th column: Physical position (in bp)
  	5th and 6th columns: (Optional) reference alleles and variant alleles

e> Time file
Time file contains the time points (generation) used to test admixture signals. The defalt time file is suggested to regular use. 

d> Individuals file
Individual files describe the basic information of the sample we used
	1st column: Individual ID
	2nd column: Gender
	3rd column: Population Label

f> Parameter file
Basic setting:
     Inputfilelist   List of input files.                           Default: -
     Indfile         Individual notation.                           Default: -
     Admpop          Label for admixed population.                  Default: -
     Refpops         Labels for reference populations.   	    Default: -
     Timefile        Time file for testing the admixture signals.   Default: -
     Outdir	     Output directory				    Default: ./

Advanced setting:
     Jackknife       Whether applied jack knife test.               Default: off
                     "1" "on"  "ON"  stands for open.
                     "0" "off" "OFF" stands for close.
     Mindis          Minimum distance of two WLD bins(in Morgan).   Defalut: 0.0002
     Maxdis          Maximum distance of two WLD bins(in Morgan).   Default: 0.3
     Binsize         Bin size of WLD bin(in Morgan)                 Default: 0.0002
     Iter            Number of iteration to the final fitting results.              
							            Default: 100000
     Num_threads     Number of threads used globally.    	    Default: 1
     Num_threads_wld Number of threads used for WLD calculation.    Default: equal to num_threads


4. Output files
iMAAPs has two main output: *.rawld and *.ad. *.rawld records the Weighted LD we calculated. *.ad records the fitting results for the weighted LD for both overall fitting and for fitting with each Jackknife. For each time point, the exponential coefficients at the time point are listed at that row. For a better summary, a R script is provided to deal with these data.

################
#Brief Tutorial#
################

Suppose you have installed iMAAPs and related package successfully, then go the example directory:
Step 1:
>../iMAAPs/iMAAPs -p ADM.par > ADM.log

You get the output file ADM.rawld and ADM.ad

Step 2: Summarize the output with R package iMAAPsmry

a) Install the package
install.package("iMAAPsmry_0.1.tar.gz");
library(iMAAPsmry);

b) Summary spectrum and admixture time interval
data(ad_inp)#or#ad_inp<-read.csv(file="../example/ADM.ad",header=T,sep="\t");
data(rawld_inp)#or#rwald_inp<-read.csv(file="../example/ADM.rawld",header=T,sep="\t");

x<-iMAAPsmry(ad_inp=ad_inp,rawld_inp=rawld_inp,denoise_cut=0.99, num_chr=10);
#summary spectrum
> x$sum_spectrum
      gen          mag
 [1,]   0 0.0017808938
 [2,]  47 0.0078082779
 [3,]  48 0.0166346230
 [4,]  49 0.0213023051
 [5,]  50 0.0241563812
 [6,]  51 0.0281666768
 [7,]  52 0.0283945524
 [8,]  53 0.0246115749
 [9,]  54 0.0195652641
[10,]  55 0.0140609798
[11,]  56 0.0093997124
[12,]  57 0.0003601793

#estimated admixture time
> x$time
       median       SD       pvalue         mag       mag_sd      prop
[1,]  0.00000 0.000000 1.152493e-08 0.001780894 0.0003210443 0.1919001
[2,] 51.48246 2.394118 2.834765e-08 0.028394552 0.0063921496 0.8080999

############################
#Tips for parameter setting#
############################

The parameters file(*.par) composes parameters to run iMAAPs:

Timefile lists the time points for testing the admixture signals. We strongly suggest you to use default value without any modification unless you knew well about the basic algorithm.
Mindis controls the starting distance of weigthed LD used for exponential fitting. This parameter should be properly set according to the marker density. Small Mindis requires dense markers.
Maxdis controls the upboundary of the distance between two markers used for exponential fitting. Maxdis should be smaller than the shortest chromosome used for inference.
Binsize controls the stepwise distance for calculating the average of weighted LD. Small Binsize requires dense markers.
Iter controls the number of Iteration to get the accurate fitting results. We suggest to use the default value.
Jackknife is used for calculating the significance of the admixture signals. For a quick scan of the data, you can put off Jackknife.
Num_threads controls the number of threads used globally, this value is no bigger than the number of the chromosomes used for inference.
Num_threads_wld controls the  number of threads used for calculating weighted LD, which is specified for the memory control when dealling NGS data. If you have large memory(>128G for each fitting process) or you should set it to a small value, say 3.


