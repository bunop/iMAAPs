/*
 * par.h
 *
 *  Created on: 2015.3.31
 *      Author: Yuan Kai
 *
 *  This file contains declaration of parameter classes to calculate the WLD and to fit the curves.
 *
 *  WldMode         Mode of WLD calculation.
 *                  5 mode are available now.
 *
 *  WldPar          Class contains parameters for WLD calculation and the following fitting process.
 *                  Could be initialization with a parameter file.
 *
 */

#ifndef PAR_H_
#define PAR_H_

# include <string>
# include <vector>

typedef enum WldMode
/*
 * Enumeration of WLD Mode
 * 5 mode available mow.
 */
{
	W,		//weight
	W2,		//weight square
	REF1,	//1 reference LD
	REF2,	//2 reference LD
	COM,	//combined LD
	Null	//Null
}WMODE;

class WldPar
/*
 * Class contains parameters for WLD calculation and the following fitting process.
 */
{
public:
	WldPar() : plist(""), ptime(""), pind(""),
	pout(""),
	pwld(""), pad(""),
	chr(),
	pops(),
	time(),
	refname(), adname(),
	popind(),
	pgeno(), psnp(),
	min_dis(-1), max_dis(-1), bin_size(-1),
	c_thd(-1),
	adp(0),
	burnin(-1),
	n_thread(-1),
	n_thread_wld(-1),
	jackknife(false), isadp(false),
	islong(false),
	mode(Null){}
	//Default construction function.

	void read(const std::string &path);
	//Read parameter from parameter file.
	void check() const;
	//Check all of these parameters.
	void print() const;
	//Print parameter on the screen.

	std::string plist, ptime, pind;					//Paths of different input files
	std::string pout;								//Path of output directory
	std::string pwld, pad;							//Paths of different output files
	std::vector < int > chr;						//Chromosome ID
	std::vector < int > pops;						//Population ID
	std::vector < int > time;						//Fitting time
	std::vector < std::string > refname, adname;	//population name of reference and admixed populations
	std::vector < int > popind;
	std::vector < std::string > pgeno, psnp;		//Paths of different genotype and SNP files
	double min_dis, max_dis, bin_size;				//minimum and maximum distance of each bin, bin size (Morgan)
	double c_thd;									//C threshold
	double adp;										//admixed proportion of pop1
	int burnin;										//burn in number of interation
	int n_thread;									//number of thread
	int n_thread_wld;								//number of thread for WLD
	bool jackknife, isadp;
	bool islong;									//whether there is chromosoeme longer than maxdis
	WMODE mode;
	/*
	 * plist       Path of input files list.
	 *             Contains path of all the genotype files and snp files.
	 * ptime       Path of time file.
	 *             Contains all of the time points used to fit the curve.
	 * pind        Path of individual file.
	 *             Contains ID, population and gender of each individuals.
	 * pout        Path of output directory
	 * pwld        Path of output WLD file.
	 * pad         Path of output admixture signal file.
	 * chr         Vector stores the chromosome ID.
	 * pops        Vector stores the population ID.
	 *             This ID stands for the position of these population stores during the program running.
	 *             Reference population is in front of admixed population
	 * time        Vector to store the time points used to fit the curves.
	 * refname     Vector of string to store population labels of reference populations.
	 * adname      Vector of string to store population labels of admixed population.
	 * popind      Vector to store the number of individuals of each populations.
	 * pgeno       Vector of string to store paths of all genotype files.
	 * psnp        Vector of string to store paths of all snp files.
	 * min_dis     Minimum distance of WLD.
	 *             In Morgan.
	 * max_dis     Maximum distance of WLD
	 *             In Morgan
	 * bin_size    Bin size during WLD calculation
	 *             In Morgan
	 * c_thd       Fitting parameter.
	 * adp         Admixed proportion of reference populations.
	 *             Combined LD model used only.
	 * burnin      Number of iteration of during fitting.
	 * n_thread    Number of thread used during WLD calculation and fitting.
	 * jackknife   Whether apply jackknife test.
	 * isadp       Whether given a prior admixed proportion.
	 */
private:
	void get_par(const std::string & line);
	//Read specific parameter into this class
	void split_str(std::vector < std::string > & names, const std::string & line);
	//Split a string by semicolon.
	void list_in();
	//Read list file
	void pop_in();
	//Read Individual file
	void time_in();
	//Read time file
};

#endif /* PAR_H_ */

