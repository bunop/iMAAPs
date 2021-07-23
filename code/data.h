/*
 * data.h
 *
 *  Created on: 2015.4.2
 *      Author: Yuan Kai
 *
 *  This head file contains declaration of data classes to calculate the Weighted LD (WLD).
 *
 *  WldData         Class used to calculate WLD.
 *                  Generic function.
 *
 *  WldChrDataBase  Base class to calculate WLD for each chromosomes.
 *                  Abstract class.
 *                  Generic function.
 *
 *  WldChrDataW     Derived class based from WldChrDataBase to calculate WLD under Weighted Model.
 *                  Generic function.
 *
 *  WldChrDataW2    Derived class based from WldChrDataBase to calculate WLD under Weighted Squared Model.
 *                  Generic function.
 *
 *  WldChrDataREF1  Derived class based from WldChrDataBase to calculate WLD under 1 Reference Model.
 *                  Generic function.
 *
 *  WldChrDataREF2  Derived class based from WldChrDataBase to calculate WLD under 2 Reference Model.
 *                  Generic function.
 *
 *  WldChrDataCOM   Derived class based from WldChrDataBase to calculate WLD under Combined LD Model.
 *                  Generic function.
 *
 */

#ifndef DATA_H_
#define DATA_H_

# include <cstdlib>
# include "par.h"

class WldData
/*
 * Class used to calculate WLD.
 *
 * A generic function could be used to calculate it automatically.
 *
 * 3 ways to use this generic function:
 *
 *     1:
 *         WldData data;
 *         WldPar p;
 *
 *         ... p initialization
 *
 *         data(p);
 *     2:
 *         WldPar p;
 *
 *         ... p initialization
 *
 *         WldData data(p);
 *         data();
 *     3:
 *         WldPar p;
 *
 *         ... p initialization
 *
 *         WLdData()(p);
 */
{
public:
	WldData() : par(NULL), wld(NULL), wld_fit(NULL), count(NULL), ad(NULL), ds(NULL), exwld(NULL), excount(NULL) ,adp(NULL){}
	//Default construction function
	WldData( WldPar & p ) throw();
	//Construction function constructed by parameters
	~WldData();
	//Deconstruction function

	//Pipeline for WLD calculation
	void operator()();
	//Initialization with WldPar
	void operator()(WldPar &p);
	//Could be the object without initialization

private:

	void read(WldPar & p) throw();
	//Read parameters
	void calculation();
	//Calculate the WLD
	void fit();
	//Fit the curve
	void jackknife();
	//Jackknife test
	void output();
	//Output results

	WldPar *par;
	double **wld, **wld_fit, **count, **ad, *ds, **exwld, **excount, *adp;
	/*
	 * par       Pointer point to an object of WldPar.
	 * wld       Two-dimensional array stores results of WLD.
	 *           [Chromosome ID][Number of bins]
	 * wld_fit   Two-dimensional array stores results of fitted WLD.
	 *           [Chromosome ID][Number of bins]
	 * count     Two-dimensional array stores results of number of count in each paired bins of different distance.
	 *           [Chromosome ID][Number of bins]
	 * ad        Two-dimensional array stores results of fitted admixed signals.
	 *           [Chromosome ID][Number of bins]
	 * ds        Array stores distance of paired bins.
	 *           cent-Morgan
	 * exwld     Two-dimensional array stores results of extra WLD.
	 *           Combined WLD Model only.
	 *           2 Reference Model WLD [Chromosome ID][Number of bins]
	 *           1 Reference Model WLD ref2->ref1 [Chromosome ID][Total number of bins + Number of bins]
	 *           1 Reference Model WLD ref2->ref1 [Chromosome ID][Total number of bins * 2 + Number of bins]
	 *           Weighted square Model WLD [Chromosome ID][Total number of bins * 3 + Number of bins]
	 * excount   Two-dimensional array stores results of extra number of count in each paired bins of different distance.
	 *           Combined WLD Model only.
	 *           2 Reference Model count [Chromosome ID][Number of bins]
	 *           1 Reference Model count ref2->ref1 [Chromosome ID][Total number of bins + Number of bins]
	 *           1 Reference Model count ref2->ref1 [Chromosome ID][Total number of bins * 2 + Number of bins]
	 *           Weighted square Model count [Chromosome ID][Total number of bins * 3 + Number of bins]
	 * adp       Array stores admixed proportion of numerator and denominator to calculate the admixed proportion
	 *           Combined WLD Model only.
	 *           numerator   [Chromosome ID * 2]
	 *           denominator [Chromosome ID * 2 + 1]
	 */
};

class WldChrDataBase
/*
 * Class used to calculate WLD for each chromosomes.
 * No dynamic space is managed by this class. Deconstructor function is not needed for this class and its derived classes in our soft.
 *
 * A generic function could be used to calculate it automatically.
 *
 */
{
public:
	WldChrDataBase() : par(NULL), chr(0), wld(NULL), count(NULL), exwld(NULL), excount(NULL), adp(NULL), pos(){}
	void operator()(WldPar &p, const int & chr_, double *wld_, double *count_, double *exwld_=NULL, double *excount_=NULL, double *adp_=NULL);
	//Pipeline to calculate the WLD of each chromosome
protected:
	virtual void calculation()=0;
	//Function for WLD calculation

	void pos_in();
	//Input position information from position file

	void ref1(double *n, double *af, double *sx, double *sx2, int *bin, char *geno, bool *skip, double *owld, double *ocount, const int & num, const int & nsite, const int & nsample,const int & nad) throw();
	//Function for calculation of WLD with 1 reference population and 1 admixed population
	void ref2(double *n, double* miss2, double *w, double *s, double *ka, char *geno, bool *skip, double *owld, double *ocount, const int & num, const int & nsite, const int & nsample) throw();
	//Function for calculation of WLD with 2 reference population and 1 admixed population.
	void weight(double *n, double *w, double *owld, double *ocount, const int & num, const int & nsample) throw();
	//Function for calculation of weight or weight square.

	WldPar *par;
	int chr;
	double *wld, *count, *exwld, *excount, *adp;
	std::vector < double > pos;
	/*
	 * par       Pointer point to an object of WldPar.
	 * chr       Number of chromosome.
	 * wld       Array to store WLD.
	 * count     Array to store count.
	 * exwld     Array to store extra WLD.
	 *           Combined LD Model only.
	 * excount   Array to store extra count.
	 *           Combined LD Model only.
	 * adp       Array to store admixed proportion.
	 *           Combined LD Modle only.
	 *
	 */
};

class WldChrDataW : public WldChrDataBase
{
private:
	void calculation() throw();
};

class WldChrDataW2 : public WldChrDataBase
{
private:
	void calculation() throw();
};

class WldChrDataREF1 : public WldChrDataBase
{
private:
	void calculation() throw();
};

class WldChrDataREF2 : public WldChrDataBase
{
private:
	void calculation() throw();
};

class WldChrDataCOM : public WldChrDataBase
{
private:
	void calculation() throw();
};

#endif /* DATA_H_ */
