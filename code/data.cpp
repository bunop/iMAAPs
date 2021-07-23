/*
 * data.cpp
 *
 *  Created on: 2015.4.3
 *      Author: Yuan Kai
 *
 *  This file contains definition of data classes to calculate the WLD.
 *
 *  WldData          void operator()()
 *                   void operator()(WldPar &)
 *                   void read(WldPar &)
 *                   void calculation()
 *                   void fit()
 *                   void jackknife()
 *                   void output()
 *
 *  WldChrDataBase   void operator()(WldPar &, const int &, double *, double *, double *, double * = NULL, double * = NULL)
 *                   void calculation()=0
 *                   void pos_in()
 *                   void ref1(double *, double *, double *, double *, int *, char *, bool *, double *, double *, const int &, const int &, const int, &const int &)
 *                   void ref2(double *, double *, double *, double *, double *, char *, bool *, double *, double *, const int &, const int &, const int &)
 *                   void weight(double *, double *, double *, double *, const int &, const int &)
 *
 *  WldChrDataW      void calculation()
 *
 *  WldChrDataW2     void calculation()
 *
 *  WldChrDartREF1   void calculation()
 *
 *  WldChrDataREF2   void calculation()
 *
 *  WldChrDataCOM    void calculation()
 *
 */

# include "data.h"
# include "omp.h"
# include "fitting.h"
# include <cstring>
# include <cmath>
# include "fftw3.h"
# include "err.h"
# include "convolution.h"
# include <iomanip>
# include <sstream>

using namespace std;

WldData::WldData(WldPar &p) throw()
/*
 * Construction function constructed by parameters
 * Allocation memory for the following computation.
 *
 * p   A reference to an object stores WLD calculation parameters.
 *
 */
{
	par=&p;
	int nchr=par->chr.size();
	int ntime=par->time.size();
	int nsample=(par->max_dis-par->min_dis)/par->bin_size;

	//Allocation memory for computation

	try
	{
		wld=new double*[nchr+1];			//0 for whole genome, 1-22... for each chromosomes
		wld_fit=new double*[nchr+1];
		count=new double*[nchr+1];
		ad=new double*[nchr+1];
		for(int i=0;i<=nchr;++i)
		{
			wld[i]=new double[nsample];
			wld_fit[i]=new double[nsample];
			count[i]=new double[nsample];
			ad[i]=new double[ntime];
		}
		ds=new double[nsample];

		//Special part for Combined WLD Model

		if(par->mode==COM)
		{
			exwld=new double*[nchr+1];
			excount=new double*[nchr+1];
			for(int i=0;i<=nchr;++i)
			{
				exwld[i]=new double[nsample*4];
				excount[i]=new double[nsample*4];
			}
			adp=new double[(nchr+1)*2];
		}
	}
	catch (std::bad_alloc &e)
	{
		err_print("Memory allocation failure!");
	}

	//Initialization of those memory

	for(int i=0;i<=nchr;++i)
	{
		memset(wld[i],0,sizeof(double)*nsample);
		memset(wld_fit[i],0,sizeof(double)*nsample);
		memset(count[i],0,sizeof(double)*nsample);
		memset(ad[i],0,sizeof(double)*ntime);
	}

	//Special part for Combined WLD Model

	if(par->mode==COM)
	{
		for(int i=0;i<=nchr;++i)
		{
			memset(exwld[i],0,sizeof(double)*nsample*4);
			memset(excount[i],0,sizeof(double)*nsample*4);
		}
		memset(adp,0,sizeof(double)*(nchr+1)*2);
	}
	for(int i=0;i<nsample;++i)
		ds[i]=par->bin_size*i*100;	//Multiply 100 converts the unit from M to cM
}

WldData::~WldData()
/*
 * Deconstruction function for WldData
 * Free memory
 */
{
	if(par!=NULL)
	{
		int nchr=par->chr.size();
		for(int i=0;i<=nchr;++i)
		{
			delete[]wld[i];
			delete[]wld_fit[i];
			delete[]count[i];
			delete[]ad[i];
		}
		delete[]wld;
		delete[]wld_fit;
		delete[]count;
		delete[]ad;
		delete[]ds;
		if(par->mode==COM)
		{
			delete[]adp;
			for(int i=0;i<=nchr;++i)
			{
				delete[]exwld[i];
				delete[]excount[i];
			}
			delete[]exwld;
			delete[]excount;
		}
	}
}

void WldData::operator ()()
/*
 * WLD calculation pipeline.
 * par must be initialized.
 */
{
	if(par==NULL)	//without initialization
		err_print("WldData type should be initialized firstly!");
	calculation();
	fit();
	output();
}

void WldData::operator ()(WldPar &p)
/*
 * WLD calculation pipeline.
 * par could be without initialization
 */
{
	read(p);
	calculation();
	fit();
	output();
}

void WldData::read(WldPar & p) throw()
/*
 * Read Parameters
 * If this->par and p are same, do nothing.
 * If this object is empty, initialize this object with parameter p.
 * If this object is not empty, clean this object and initialize it with p.
 * It equals to combine ~WldData() and WldData(WldPar).
 *
 * p   A reference to an object stores WLD calculation parameters.
 *
 */
{
	if(par==&p)
	//Same parameter
		return;

	//Deconstruction function

	if(par!=NULL)
	//Not empty
	{
		int nchr=par->chr.size();
		for(int i=0;i<=nchr;++i)
		{
			delete[]wld[i];
			delete[]wld_fit[i];
			delete[]count[i];
			delete[]ad[i];
		}
		delete[]wld;
		delete[]wld_fit;
		delete[]count;
		delete[]ad;
		delete[]ds;
		if(par->mode==COM)
		{
			delete[]adp;
			for(int i=0;i<=nchr;++i)
			{
				delete[]exwld[i];
				delete[]excount[i];
			}
			delete[]exwld;
			delete[]excount;
		}
	}

	//Construction function

	par=&p;
	int nchr=par->chr.size();
	int ntime=par->time.size();
	int nsample=(par->max_dis-par->min_dis)/par->bin_size;

	//Allocation memory for computation

	try
	{
		wld=new double*[nchr+1];			//0 for whole genome, 1-22... for each chromosomes
		wld_fit=new double*[nchr+1];
		count=new double*[nchr+1];
		ad=new double*[nchr+1];
		for(int i=0;i<=nchr;++i)
		{
			wld[i]=new double[nsample];
			wld_fit[i]=new double[nsample];
			count[i]=new double[nsample];
			ad[i]=new double[ntime];
		}
		ds=new double[nsample];

		//Special part for Combined WLD Model

		if(par->mode==COM)
		{
			exwld=new double*[nchr+1];
			excount=new double*[nchr+1];
			for(int i=0;i<=nchr;++i)
			{
				exwld[i]=new double[nsample*4];
				excount[i]=new double[nsample*4];
			}
			adp=new double[(nchr+1)*2];
		}
	}
	catch (std::bad_alloc &e)
	{
		err_print("Memory allocation failure!");
	}

	//Initialization of those memory

	for(int i=0;i<=nchr;++i)
	{
		memset(wld[i],0,sizeof(double)*nsample);
		memset(wld_fit[i],0,sizeof(double)*nsample);
		memset(count[i],0,sizeof(double)*nsample);
		memset(ad[i],0,sizeof(double)*ntime);
	}

	//Special part for Combined WLD Model

	if(par->mode==COM)
	{
		for(int i=0;i<=nchr;++i)
		{
			memset(exwld[i],0,sizeof(double)*nsample*4);
			memset(excount[i],0,sizeof(double)*nsample*4);
		}
		memset(adp,0,sizeof(double)*(nchr+1)*2);
	}
	for(int i=0;i<nsample;++i)
		ds[i]=par->bin_size*i*100;	//Multiply 100 converts the unit from M to cM
}

void WldData::calculation()
/*
 * Calculate the WLD for each chromosomes and whole genome.
 * WLD of different chromosomes are calculated in parallel with OpenMP.
 * WLD of whole genome are calculated based on results of each chromosomes.
 *
 */
{
	//Calculate the WLD for each chromosomes in parallel

	int nchr=par->chr.size();
	omp_set_num_threads(par->n_thread_wld);	//Set Number of thread to calculate the WLD
	cout<<"\nStart to calculate the WLD.\n"<<endl;
	cout<<"Chromosome:\t"<<flush;
#pragma omp parallel for
	for(int i=1;i<=nchr;++i)
	{
		//Here we used temporary variables to calculate the WLD for each chromosome to make sure the temporary resources were released in time.
		if(par->mode==W)			//Weighted Model
			WldChrDataW()(*par, i, wld[i], count[i]);
		else if(par->mode==W2)		//Weighted Square Model
			WldChrDataW2()(*par, i ,wld[i], count[i]);
		else if(par->mode==REF1)	//1 Reference Model
			WldChrDataREF1()(*par, i, wld[i], count[i]);
		else if(par->mode==REF2)	//2 Reference Model
			WldChrDataREF2()(*par, i, wld[i], count[i]);
		else if(par->mode==COM)		//Combined LD Model
			WldChrDataCOM()(*par, i, wld[i], count[i], exwld[i], excount[i], adp+(i*2));
#pragma omp critical
		cout<<par->chr[i-1]<<'\t'<<flush;
	}
	cout<<"\nWLD calculation done!\n"<<endl;

	if(!par->islong)
		err_print("No chromoseome's length is longer than maximum distance!");

	//Calculate the WLD for whole genome

	int nsample=(par->max_dis-par->min_dis)/par->bin_size;
	if(par->mode==COM)	//Combined LD
	{
		for(int i=1;i<=nchr;++i)
		{
			for(int j=0;j<4*nsample;++j)
			{
				exwld[0][j]+=exwld[i][j];
				excount[0][j]+=excount[i][j];
			}
			adp[0]+=adp[i*2];
			adp[1]+=adp[i*2+1];
		}
		adp[0]/=adp[1];
		if(par->isadp)
			adp[0]=par->adp;
		cout<<"Admixed proportion"<<endl;
		cout<<par->refname[0]<<":"<<par->refname[1]<<"="<<adp[0]<<":"<<1-adp[0]<<endl;
		for(int i=0;i<nsample;++i)
		{
			if(excount[0][i]<0.01 || excount[0][i+nsample]<0.01 || excount[0][i+nsample*2]<0.01 || excount[0][i+nsample*3]<0.01 )	//Any counts of 4 models is 0, Combined WLD for this bins is set to be 0.
				wld[0][i]=0;
			else
				wld[0][i]=(exwld[0][i]/excount[0][i]-exwld[0][i+nsample]/excount[0][i+nsample]*adp[0]-exwld[0][i+nsample*2]/excount[0][i+nsample*2]*(1-adp[0]))/(exwld[0][i+nsample*3]/excount[0][i+nsample*3]);
				//Combined LD = ( LD ( ref1&ref2 -> ad ) - LD ( ref2 -> ref1 ) * adp ( ref1 ) - LD ( ref1 -> ref2 ) * adp ( ref2 ) ) / LD-W2 ( ref1&ref2 )
			count[0][i]=1;
		}
		//Calculate for jackknife test
		for(int j=1;j<=nchr;++j)
		{
			for(int i=0;i<nsample;++i)
			{
				double tc1,tc2,tc3,tc4;
				double tw1,tw2,tw3,tw4;
				tc1=excount[0][i]-excount[j][i];
				tc2=excount[0][i+nsample]-excount[j][i+nsample];
				tc3=excount[0][i+nsample*2]-excount[j][i+nsample*2];
				tc4=excount[0][i+nsample*3]-excount[j][i+nsample*3];
				tw1=exwld[0][i]-exwld[j][i];
				tw2=exwld[0][i+nsample]-exwld[j][i+nsample];
				tw3=exwld[0][i+nsample*2]-exwld[j][i+nsample*2];
				tw4=exwld[0][i+nsample*3]-exwld[j][i+nsample*3];
				if(tc1<0.01 || tc2<0.01 || tc3<0.01 || tc4<0.01)
					wld[j][i]=0;
				else
					wld[j][i]=(tw1/tc1-tw2/tc2*adp[0]-tw3/tc3*(1-adp[0]))/(tw4/tc4);
				count[j][i]=1;
			}
		}
		for(int i=0;i<4*nsample;++i)
			if(excount[0][i]<0.01)
				exwld[0][i]=0;
			else
				exwld[0][i]/=excount[0][i];
	}
	else	//Non-Combined LD
	{
		for(int i=1;i<=nchr;++i)
		{
			for(int j=0;j<nsample;++j)
			{
				wld[0][j]+=wld[i][j];
				count[0][j]+=count[i][j];
			}
		}
		for(int i=1;i<=nchr;++i)
		{
			for(int j=0;j<nsample;++j)
			{
				wld[i][j]=wld[0][j]-wld[i][j];
				count[i][j]=count[0][j]-count[i][j];
				if(count[i][j]<0.01)
					wld[i][j]=0;
				else
					wld[i][j]/=count[i][j];
			}
		}
		for(int i=0;i<nsample;++i)
			if(count[0][i]<0.01)
				wld[0][i]=0;
			else
				wld[0][i]/=count[0][i];
	}
}

void WldData::fit()
/*
 * Fit the WLD decay curve.
 *
 */
{
	int nbin=(par->max_dis-par->min_dis)/par->bin_size;
	int ntime=par->time.size();
	int cut=ceil(par->min_dis/par->bin_size);
	cout<<"Start to fit the curve.\n"<<endl;
	if(par->jackknife)	//jackknife test
		jackknife();
	else				//Only for the whole genome
	{
		fitting(wld[0]+cut,ds+cut,wld_fit[0]+cut,nbin-cut,par->time, ad[0],ntime, par->burnin, par->c_thd);
		cout<<"Whole genome fitting done!"<<endl;
	}
	cout<<"Fitting done!\n"<<endl;
}

void WldData::jackknife()
/*
 * Jackknife test.
 *
 */
{
	int nchr=par->chr.size();
	int nsample=(par->max_dis-par->min_dis)/par->bin_size;
	int ntime=par->time.size();
	int cut=ceil(par->min_dis/par->bin_size);
	omp_set_num_threads(par->n_thread);
#pragma omp parallel for
	for(int i=0;i<=nchr;++i)
	{
		fitting(wld[i]+cut,ds+cut,wld_fit[i]+cut,nsample-cut,par->time, ad[i],ntime, par->burnin, par->c_thd);
#pragma omp critical
		if(i==0)
			cout<<"Whole genome fitting done!"<<endl;
		else
			cout<<"Chromosome "<<par->chr[i-1]<<" done!"<<endl;
	}
}

void WldData::output()
/*
 * Output results
 * Both WLD and admixed signal.
 *
 */
{
	int nchr=par->chr.size();
	int nsample=(par->max_dis-par->min_dis)/par->bin_size;
	int ntime=par->time.size();

	//Output WLD file

	cout<<"WLD file: "<<par->pwld<<endl;
	ofstream fpwld(par->pwld.c_str());
	if(par->mode==COM)
	{
		fpwld<<"Distance\tCombined_LD\tFitted\tWLD_2ref\tCount\tWLD_1ref_"<<par->refname[1]<<"->"<<par->refname[0]<<
		"\tCount\tWLD_1ref_"<<par->refname[0]<<"->"<<par->refname[1]<<"\tCount\tWeight_Square\tCount";
		if(par->jackknife)
			for(int i=1;i<=nchr;++i)
				fpwld<<"\tJack"<<i;
		fpwld<<endl;
		for(int i=1;i<nsample;++i)
		{
			double dis=i*par->bin_size;
			if(dis<par->min_dis)
				fpwld<<'#';
			fpwld<<setprecision(6)<<fixed<<dis<<'\t'<<setprecision(16)<<
			wld[0][i]<<'\t'<<wld_fit[0][i]<<'\t'<<
			setprecision(16)<<exwld[0][i]<<'\t'<<static_cast<int>(excount[0][i]+0.5)<<'\t'<<
			exwld[0][i+nsample]<<'\t'<<static_cast<int>(excount[0][i+nsample]+0.5)<<'\t'<<
			exwld[0][i+nsample*2]<<'\t'<<static_cast<int>(excount[0][i+nsample*2]+0.5)<<'\t'<<
			exwld[0][i+nsample*3]<<'\t'<<static_cast<int>(excount[0][i+nsample*3]+0.5);
			if(par->jackknife)
			{
				for(int j=1;j<=nchr;++j)
					fpwld<<'\t'<<wld[j][i];
			}
			fpwld<<endl;
		}
	}
	else
	{
		fpwld<<"Distance\tWLD\tCount\tFitted";
		if(par->jackknife)
			for(int i=1;i<=nchr;++i)
				fpwld<<"\tJack"<<i;
		fpwld<<endl;
		for(int i=1;i<nsample;++i)
		{
			double dis=i*par->bin_size;
			if(dis<par->min_dis)
				fpwld<<'#';
			fpwld<<setprecision(6)<<fixed<<dis<<'\t'<<setprecision(16)<<
			wld[0][i]<<'\t'<<static_cast<int>(count[0][i])<<'\t'<<wld_fit[0][i];
			if(par->jackknife)
			{
				for(int j=1;j<=nchr;++j)
					fpwld<<'\t'<<wld[j][i];
			}
			fpwld<<endl;
		}
	}

	//Output Admixed signal file

	cout<<"Admixed file: "<<par->pad<<endl;
	ofstream fpad(par->pad.c_str());
	fpad<<setprecision(10)<<setiosflags(ios::fixed);
	fpad<<"Time(generation)\tAll";
	if(par->jackknife)
		for(int i=0;i<nchr;++i)
			fpad<<"\tJack"<<i+1;
	fpad<<endl;
	for(int i=0;i<ntime;++i)
	{
		fpad<<par->time[i]<<'\t'<<ad[0][i];
		if(par->jackknife)
			for(int j=1;j<=nchr;++j)
				fpad<<'\t'<<ad[j][i];
		fpad<<endl;
	}
}

void WldChrDataBase::operator()(WldPar &p, const int &chr_, double *wld_, double *count_, double *exwld_, double *excount_, double * adp_)
/*
 * Pipeline for WLD calculation of each chromosomes.
 *
 * p          Parameters for calculation
 * chr_       Chromosome number. Used to get corrected path of geno and position files
 * wld_       Array used to store WLD
 * count_     Array used to store count of paired bins
 * exwld_     Array to store extra WLD information
 *            Combined LD Model Only
 * excount_   Array to store extra count information
 *            Combined LD Model Only
 * adp_       Array to store admixed proportion information
 *            Combined LD Model Only
 */
{
	par=&p;
	chr=chr_;
	wld=wld_;
	count=count_;
	exwld=exwld_;
	excount=excount_;
	adp=adp_;
	pos_in();
	calculation();
}

void WldChrDataBase::pos_in()
/*
 * Input position information from position file
 */
{
	ifstream fpi(par->psnp[chr].c_str());
	pos.clear();
	double tpos;
	string tmp,line;
	int tchr;
	while(getline(fpi,line))
	{
		istringstream cl(line);
		cl>>tmp>>tchr>>tpos;
		if(tchr!=chr)
		{
			cerr<<par->psnp[chr]<<"\nLine:"<<pos.size()+1<<endl;
			cerr<<line<<endl;
			err_print("Chromosome ID not corrected!");
		}
		if(tpos<0)
		{
			cerr<<par->psnp[chr]<<"\nLine:"<<pos.size()+1<<endl;
			cerr<<line<<endl;
			err_print("Position cannot be negative!");
		}
		pos.push_back(tpos);
	}
	if(pos.back() > par->max_dis)
		par->islong=true;
}

void WldChrDataBase::ref1(double *n, double *af, double *sx, double *sx2, int *bin, char *geno, bool *skip, double *owld, double *ocount, const int & num, const int & nsite, const int & nsample,const int &nad) throw()
/*
 * Function for calculation of WLD with 1 reference population and 1 admixed population.
 *
 * n       Array to store number sites for each bins.
 *         Length: Magnitude of FFTW plan.
 * af      Array to store allele frequency of reference population.
 *         Length: Total number of segregating sites of each chromosome.
 * sx      Array to store summation of reference(alternative) alleles of one site.
 *         Length: ToTal number of segregating sites of each chromosome.
 * sx2     Array to store summation of square of reference(alternative) alleles of one site.
 *         Length: ToTal number of segregating sites of each chromosome.
 * bin     Array to store which bins belong of each sites.
 *         Length: ToTal number of segregating sites of each chromosome.
 * geno    Array to store genotype information of admixed population.
 *         Length: ToTal number of segregating sites of each chromosome multiply number of individuals in admixed population.
 * skip    Array to store whether this site will be ignore.
 *         Length: ToTal number of segregating sites of each chromosome.
 * owld    Array to store the output results of WLD.
 *         Length: Number of samples to fit the curve.
 * ocount  Array to store the count number of WLD.
 *         Length: Number of samples to fit the curve.
 * num     Magnitude of FFTW plan.
 * nsite   Number of segregating sites of each chromosome.
 * nsample Number of samples to fit the curve.
 * nad     Number of individuals in admixed population.
 *
 * The simplification of FFTW calculation refer to Alder.
 *
 */
{
	double *in1, *in2, *out;
	fftw_complex *rev, *t1, *t2;
	try
	{
		in1=new double[num];
		in2=new double[num];
		out=new double[num];
		rev=new fftw_complex[num];
		t1=new fftw_complex[num];
		t2=new fftw_complex[num];
	}
	catch (std::bad_alloc &e)
	{
		err_print("Memory allocation failure!");
	}
	memset(out,0,sizeof(double)*num);
	memset(rev,0,sizeof(fftw_complex)*num);
	memset(t1,0,sizeof(fftw_complex)*num);
	memset(t2,0,sizeof(fftw_complex)*num);

	fftw_plan plan[2];
	fftw_plan revd;
	fftw_plan splan;

#pragma omp critical
	{
		plan[0]=fftw_plan_dft_r2c_1d(num,in1,t1,FFTW_ESTIMATE);
		plan[1]=fftw_plan_dft_r2c_1d(num,in2,t2,FFTW_ESTIMATE);
		revd=fftw_plan_dft_c2r_1d(num,rev,out,FFTW_ESTIMATE);
		splan=fftw_plan_dft_r2c_1d(num,n,t1,FFTW_ESTIMATE);
	}
	int p(0);

	double so1(nad), so2(so1*(nad-1)), so3(so2*(nad-2)), so4(so3*(nad-3));

//-pAx * pAy * S10 * S01 / SO2
//1

	memset(in1,0,sizeof(double)*num);
	for(int i=0;i<nsite;++i)
		if(!skip[i])
			in1[bin[i]]+=af[i]*sx[i];
	self_convolution(plan,rev,t1,num,-0.5/so2);

//pAx *	S10 * ( S01 *S01 - S02 ) / SO3
//3,5,7,8

	memset(in2,0,sizeof(double)*num);
	for(int i=0;i<nsite;++i)
		if(!skip[i])
			in2[bin[i]]+=sx[i]*sx[i]-sx2[i];
	convolution(plan,rev,t1,t2,num,0.5/so3);

//(2 * S20 - S10 * S10 ) * S01 * S01 / SO4
//11,12,15

	memset(in1,0,sizeof(double)*num);
	memset(in2,0,sizeof(double)*num);
	for(int i=0;i<nsite;++i)
		if(!skip[i])
		{
			in1[bin[i]]+=2*sx2[i]-sx[i]*sx[i];
			in2[bin[i]]+=sx[i]*sx[i];
		}
	convolution(plan,rev,t1,t2,num,0.125/so4);

//-S02 * S20 / SO4
//18

	memset(in1,0,sizeof(double)*num);
	for(int i=0;i<nsite;++i)
		if(!skip[i])
			in1[bin[i]]+=sx2[i];
	self_convolution(plan,rev,t1,num,-0.125/so4);

	for(int i=0;i<nad;++i)
	{
//pAx * pAy * S11 * ( SO2 + SO1 ) / SO2 / SO1
//2

		p=i*nsite;
		memset(in1,0,sizeof(double)*num);
		for(int j=0;j<nsite;++j)
		{
			if(!skip[j])
				in1[bin[j]]+=af[j]*geno[p];
			++p;
		}
		self_convolution(plan,rev,t1,num,0.5*(so2+so1)/so1/so2);

//pAx * ( S12 - S11 * S01 ) * ( 2 * SO2 + SO3 ) / SO2 / SO3
//4,6,9,10

		p=i*nsite;
		memset(in2,0,sizeof(double)*num);
		for(int j=0;j<nsite;++j)
		{
			if(!skip[j])
				in2[bin[j]]+=geno[p]*(geno[p]-sx[j]);
			++p;
		}
		convolution(plan,rev,t1,t2,num,0.5*(2*so2+so3)/so3/so2);

//S01 * ( S10 * S11 - 2 * S21 ) * (4 * SO3 + SO4 ) / SO3 / SO4
//13,14,17

		p=i*nsite;
		memset(in1,0,sizeof(double)*num);
		memset(in2,0,sizeof(double)*num);
		for(int j=0;j<nsite;++j)
		{
			if(!skip[j])
			{
				in1[bin[j]]+=sx[j]*geno[p];
				in2[bin[j]]+=geno[p]*(sx[j]-2*geno[p]);
			}
			++p;
		}
		convolution(plan,rev,t1,t2,num,0.125*(4*so3+so4)/so3/so4);

//2 * S22 * ( 3 * SO3 + SO4 ) / SO3 / SO4
//19

		p=i*nsite;
		memset(in1,0,sizeof(double)*num);
		for(int j=0;j<nsite;++j)
		{
			if(!skip[j])
				in1[bin[j]]+=geno[p]*geno[p];
			++p;
		}
		self_convolution(plan,rev,t1,num,0.125 *(4 *so3+so4)/so3/so4);
	}

//-S11 * S11 * ( 2 * SO3 +SO4 ) / SO3 / SO4
//16

	for(int i=0;i<nad;++i)
		for(int j=i+1;j<nad;++j)
		{
			int p1(i*nsite), p2(j*nsite);
			memset(in1,0,sizeof(double)*num);
			for(int k=0;k<nsite;++k)
			{
				if(!skip[k])
					in1[bin[k]]+=geno[p1]*geno[p2];
				++p1;
				++p2;
			}
			self_convolution(plan,rev,t1,num,-0.25*(2*so3+so4)/so3/so4);
		}

	fftw_execute(revd);
	for(int i=1;i<nsample && i<=num;++i)
		owld[i]=(out[i]+out[num-i])/num;

	memset(rev,0,sizeof(double)*num);

	self_convolution( &splan, rev, t1, num, 1.0);
	fftw_execute(revd);
	for(int i=1;i<nsample && i<=num;++i)
		ocount[i]=out[i]/num;

#pragma omp critical
	{
		fftw_destroy_plan(plan[0]);
		fftw_destroy_plan(plan[1]);
		fftw_destroy_plan(splan);
		fftw_destroy_plan(revd);
	}

	delete[]in1;
	delete[]in2;
	delete[]out;
	delete[]rev;
	delete[]t1;
	delete[]t2;
}

void WldChrDataBase::ref2(double *n, double *miss2, double *w, double *s, double *ka, char *geno, bool *skip, double *owld, double *ocount, const int & num, const int & nsite, const int & nsample) throw()
/*
 * Function for calculation of WLD with 2 reference population and 1 admixed population.
 *
 * n       Array to store number sites for each bins.
 *         Length: Magnitude of FFTW plan.
 * miss2   Array to store number of missing sites for each bins.
 *         Length: Magnitude of FFTW plan.
 * w       Array to store allele frequency differences of two reference populations.
 *         Length: Total number of segregating sites of each chromosome.
 * s       Array to store summation of reference(alternative) alleles of one site.
 *         Length: ToTal number of segregating sites of each chromosome.
 * ka      Array to store number of missing alleles of one site.
 *         Length: ToTal number of segregating sites of each chromosome.
 * geno    Array to store genotype information of admixed population.
 *         Length: ToTal number of segregating sites of each chromosome multiply number of individuals in admixed population.
 * skip    Array to store whether this site will be ignore.
 *         Length: ToTal number of segregating sites of each chromosome.
 * owld    Array to store the output results of WLD.
 *         Length: Number of samples to fit the curve.
 * ocount  Array to store the count number of WLD.
 *         Length: Number of samples to fit the curve.
 * num     Magnitude of FFTW plan.
 * nsite   Number of segregating sites of each chromosome.
 * nsample Number of samples to fit the curve.
 *
 * The simplification of FFTW calculation refer to Alder.
 *
 */
{
	int nad=par->popind[3];

	double *in1, *in2, *out;
	fftw_complex *rev, *t1, *t2;
	try
	{
		in1=new double[num];
		in2=new double[num];
		out=new double[num];
		rev=new fftw_complex[num];
		t1=new fftw_complex[num];
		t2=new fftw_complex[num];
	}
	catch (std::bad_alloc &e)
	{
		err_print("Memory allocation failure!");
	}
	memset(out,0,sizeof(double)*num);
	memset(rev,0,sizeof(fftw_complex)*num);
	memset(t1,0,sizeof(fftw_complex)*num);
	memset(t2,0,sizeof(fftw_complex)*num);

	fftw_plan plan[2];
	fftw_plan revd;
	fftw_plan splan;
	fftw_plan splan2;

#pragma omp critical
	{
		plan[0]=fftw_plan_dft_r2c_1d(num,in1,t1,FFTW_ESTIMATE);
		plan[1]=fftw_plan_dft_r2c_1d(num,in2,t2,FFTW_ESTIMATE);
		revd=fftw_plan_dft_c2r_1d(num,rev,out,FFTW_ESTIMATE);
		splan=fftw_plan_dft_r2c_1d(num,n,t1,FFTW_ESTIMATE);
		splan2=fftw_plan_dft_r2c_1d(num,miss2,t1,FFTW_ESTIMATE);
	}

	int p(0);

	for(int i=0;i<nad;++i)
	{
		memset(in1,0,sizeof(double)*num);
		memset(in2,0,sizeof(double)*num);
		for(int j=0;j<nsite;++j)
		{
			int bin=floor(pos[j]/par->bin_size);
			if(skip[j])
			{
				++p;
				continue;
			}
			if(ka[j])
			{
				if(geno[p]==9)
					in2[bin]+=s[j]*w[j]/(nad-ka[j])/(nad-ka[j]-1);
				else
					in2[bin]+=geno[p]*w[j]/(nad-ka[j]-1);
			}
			else
			{
				in1[bin]+=geno[p]*w[j];
				in2[bin]+=0.5*geno[p]*w[j]/(nad-1);
			}
			++p;
		}
		convolution(plan,rev,t1,t2,num,1.0);
	}
	memset(in1,0,sizeof(double)*num);
	memset(in2,0,sizeof(double)*num);
	for(int i=0;i<nsite;++i)
	{
		int bin=floor(pos[i]/par->bin_size);
		if(skip[i])
			continue;
		if(ka[i])
			in2[bin]+=s[i]*w[i]/(nad-ka[i])/(nad-ka[i]-1);
		else
		{
			in1[bin]+=s[i]*w[i];
			in2[bin]+=0.5*s[i]*w[i]/nad/(nad-1);
		}
	}
	convolution(plan,rev,t1,t2,num,-1.0);

	fftw_execute(revd);
	for(int i=1;i<nsample && i<=num;++i)
		owld[i]=(out[i]+out[num-i])/num;

	memset(rev,0,sizeof(fftw_complex)*num);

	self_convolution( &splan, rev, t1, num, 1.0);
	self_convolution( &splan2, rev, t1, num, -1.0);
	fftw_execute(revd);
	for(int i=1;i<nsample && i<=num;++i)
		ocount[i]=out[i]/num;

#pragma omp critical
	{
		fftw_destroy_plan(plan[0]);
		fftw_destroy_plan(plan[1]);
		fftw_destroy_plan(splan);
		fftw_destroy_plan(splan2);
		fftw_destroy_plan(revd);
	}

	delete[]in1;
	delete[]in2;
	delete[]out;
	delete[]rev;
	delete[]t1;
	delete[]t2;
}

void WldChrDataBase::weight(double *n, double *w, double *owld, double *ocount, const int & num, const int & nsample) throw()
/*
 * Function for calculation of weight or weight square.
 *
 * n       Array to store number sites for each bins.
 *         Length: Magnitude of FFTW plan.
 * w       Array to store allele frequency differences of two reference populations.
 *         Length: Magnitude of FFTW plan.
 * owld    Array to store the output results of WLD.
 *         Length: Number of samples to fit the curve.
 * ocount  Array to store the count number of WLD.
 *         Length: Number of samples to fit the curve.
 * num     Magnitude of FFTW plan.
 * nsample Number of samples to fit the curve.
 *
 */
{
	double *out;
	fftw_complex *rev, *t;
	try
	{
		out=new double[num];
		rev=new fftw_complex[num];
		t=new fftw_complex[num];
	}
	catch (std::bad_alloc &e)
	{
		err_print("Memory allocation failure!");
	}
	memset(out,0,sizeof(double)*num);
	memset(rev,0,sizeof(fftw_complex)*num);
	memset(t,0,sizeof(fftw_complex)*num);

	fftw_plan wplan;
	fftw_plan splan;
	fftw_plan revd;

#pragma omp critical
	{
		wplan=fftw_plan_dft_r2c_1d(num,w,t,FFTW_ESTIMATE);
		splan=fftw_plan_dft_r2c_1d(num,n,t,FFTW_ESTIMATE);
		revd=fftw_plan_dft_c2r_1d(num,rev,out,FFTW_ESTIMATE);
	}

	self_convolution( &wplan, rev, t, num, 1.0);
	fftw_execute(revd);
	for(int i=1;i<nsample && i<=num;++i)
		owld[i]=(out[i]+out[num-i])/num;

	memset(rev,0,sizeof(fftw_complex)*num);

	self_convolution( &splan, rev, t, num, 1.0);
	fftw_execute(revd);
	for(int i=1;i<nsample && i<=num;++i)
		ocount[i]=out[i]/num;

#pragma omp critical
	{
		fftw_destroy_plan(wplan);
		fftw_destroy_plan(splan);
		fftw_destroy_plan(revd);
	}

	delete[]out;
	delete[]rev;
	delete[]t;
}

void WldChrDataW::calculation() throw()
/*
 * Calculation function for derived class WldChrDataW (Weight Model)
 *
 */
{
	int nbin;
	nbin=ceil((pos.back()-pos.front())/par->bin_size);	//Number of bins
	int shift=0;
	while(1<<shift < nbin)
		++shift;
	++shift;
	int num=1<<shift;	//Magnitude of FFTW plan

	double *w, *n;		//Weight and Number of sites in each bins
	try
	{
		w=new double[num];
		n=new double[num];
	}
	catch (std::bad_alloc &e)
	{
		err_print("Memory allocation failure!");
	}
	memset(w,0,sizeof(double)*num);
	memset(n,0,sizeof(double)*num);

	ifstream fpi(par->pgeno[chr].c_str());
	string snp;
	int nsite=pos.size();
	int nind=par->pops.size();
	vector<int> & tpops(par->pops);
	for(int i=0;i<nsite;++i)
	{
		getline(fpi,snp);
		int bin=floor((pos[i]-pos.front())/par->bin_size);
		double count1(0),nallele1(0),count2(0),nallele2(0);
		if(snp.size()!=nind)
		{
			cerr<<par->pgeno[chr]<<"\nLine:"<<i+1<<endl;
			cerr<<"There are "<<snp.size()<<" SNPs, while your dataset contains "<<nind<<" indiviudlas!"<<endl;
			err_print("Pls check your genotype files!");
		}
		for(int j=0;j<nind;++j)
		{
			int tsnp=snp[j]-'0';
			if(tpops[j]==1 && tsnp!=9)
			{
				++count1;
				nallele1+=tsnp;
			}
			else if(tpops[j]==2 && tsnp!=9)
			{
				++count2;
				nallele2+=tsnp;
			}
		}
		if(count1 && count2)
		{
			w[bin]+=(nallele1/count1/2-nallele2/count2/2);
			n[bin]++;
		}
	}
	getline(fpi,snp);
	if(fpi && snp!="")
	{
		cerr<<par->pgeno[chr]<<endl;
		cerr<<"There are "<<nsite<<" SNPs, while more markers in your genotype file!"<<endl;
		err_print("Pls check your genotype files!");
	}
	int nsample=(par->max_dis-par->min_dis)/par->bin_size;

	weight(n, w, wld, count, num, nsample);

	delete[]w;
	delete[]n;
}

void WldChrDataW2::calculation() throw()
/*
 * Calculation function for derived class WldChrDataW2 (Weight Square Model)
 *
 */
{
	int nbin;
	nbin=ceil((pos.back()-pos.front())/par->bin_size);
	int shift=0;
	while(1<<shift < nbin)
		++shift;
	++shift;
	int num=1<<shift;

	double *w2, *n;
	try
	{
		w2=new double[num];
		n=new double[num];
	}
	catch (std::bad_alloc &e)
	{
		err_print("Memory allocation failure!");
	}
	memset(w2,0,sizeof(double)*num);
	memset(n,0,sizeof(double)*num);

	ifstream fpi(par->pgeno[chr].c_str());
	string snp;
	int nsite=pos.size();
	int nind=par->pops.size();
	vector<int> & tpops(par->pops);
	for(int i=0;i<nsite;++i)
	{
		getline(fpi,snp);
		int bin=floor((pos[i]-pos.front())/par->bin_size);
		double count1(0),nallele1(0),count2(0),nallele2(0);
		if(snp.size()!=nind)
		{
			cerr<<par->pgeno[chr]<<"\nLine:"<<i+1<<endl;
			cerr<<"There are "<<snp.size()<<" SNPs, while your dataset contains "<<nind<<" indiviudlas!"<<endl;
			err_print("Pls check your genotype files!");
		}
		for(int j=0;j<nind;++j)
		{
			int tsnp=snp[j]-'0';
			if(tpops[j]==1 && tsnp!=9)
			{
				++count1;
				nallele1+=tsnp;
			}
			else if(tpops[j]==2 && tsnp!=9)
			{
				++count2;
				nallele2+=tsnp;
			}
		}
		if(count1 && count2)
		{
			double tmpw=(nallele1/count1/2-nallele2/count2/2);
			w2[bin]+=tmpw*tmpw;
			n[bin]++;
		}
	}
	getline(fpi,snp);
	if(fpi && snp!="")
	{
		cerr<<par->pgeno[chr]<<endl;
		cerr<<"There are "<<nsite<<" SNPs, while more markers in your genotype file!"<<endl;
		err_print("Pls check your genotype files!");
	}
	int nsample=(par->max_dis-par->min_dis)/par->bin_size;

	weight(n, w2, wld, count, num, nsample);

	delete[]w2;
	delete[]n;
}

void WldChrDataREF1::calculation() throw()
/*
 * Calculation function for derived class WldChrDataREF1 (1 Reference Model)
 *
 */
{
	int nbin=ceil((pos.back()-pos.front())/par->bin_size);
	int nsite=pos.size();
	int nad=par->popind[2];
	int shift=0;
	while(1<<shift < nbin)
		++shift;
	++shift;
	int num=1<<shift;

	double *n, *af, *sx, *sx2;
	char *geno;
	int *bin;
	bool *skip;
	try
	{
		n=new double[num];
		af=new double[nsite];
		sx=new double[nsite];
		sx2=new double[nsite];
		geno=new char[nad*nsite];
		bin=new int[nsite];
		skip=new bool[nsite];
	}
	catch (std::bad_alloc &e)
	{
		err_print("Memory allocation failure!");
	}
	memset(n,0,sizeof(double)*num);
	memset(af,0,sizeof(double)*nsite);
	memset(sx,0,sizeof(double)*nsite);
	memset(sx2,0,sizeof(double)*nsite);
	memset(geno,0,sizeof(char)*nad*nsite);
	memset(bin,0,sizeof(int)*nsite);
	memset(skip,0,sizeof(bool)*nsite);

	int nsample=(par->max_dis-par->min_dis)/par->bin_size;

	ifstream fpi(par->pgeno[chr].c_str());
	string snp;
	int p(0);
	int nind=par->pops.size();
	vector<int> & tpops(par->pops);
	for(int i=0;i<nsite;++i)
	{
		p=i;
		getline(fpi,snp);
		bin[i]=floor((pos[i]-pos.front())/par->bin_size);
		double count1(0),nallele1(0),tsx(0),tsx2(0),counta(0);
		bool miss(false);
		if(snp.size()!=nind)
		{
			cerr<<par->pgeno[chr]<<"\nLine:"<<i+1<<endl;
			cerr<<"There are "<<snp.size()<<" SNPs, while your dataset contains "<<nind<<" indiviudlas!"<<endl;
			err_print("Pls check your genotype files!");
		}
		for(int j=0;j<nind;++j)
		{
			int tsnp=snp[j]-'0';
			if(tpops[j]==1 && tsnp!=9)
			{
				++count1;
				nallele1+=tsnp;
			}
			else if(tpops[j]==2)
			{
				geno[p]=tsnp;
				if(tsnp!=9)
				{
					++counta;
					tsx+=tsnp;
					tsx2+=tsnp*tsnp;
				}
				else
				{
					miss=true;
					break;
				}
				p+=nsite;
			}
		}
		if(count1==0 || miss)
			skip[i]=true;
		else
		{
			af[i]=nallele1/count1/2;
			sx[i]=tsx;
			sx2[i]=tsx2;
			++n[bin[i]];
		}
	}
	getline(fpi,snp);
	if(fpi && snp!="")
	{
		cerr<<par->pgeno[chr]<<endl;
		cerr<<"There are "<<nsite<<" SNPs, while more markers in your genotype file!"<<endl;
		err_print("Pls check your genotype files!");
	}

	ref1(n, af, sx, sx2, bin, geno, skip, wld, count, num, nsite, nsample,par->popind[2]);

	delete[]n;
	delete[]af;
	delete[]sx;
	delete[]sx2;
	delete[]geno;
	delete[]bin;
	delete[]skip;
}

void WldChrDataREF2::calculation() throw()
/*
 * Calculation function for derived class WldChrDataREF2 (2 Reference Model)
 *
 */
{
	int nbin=ceil((pos.back()-pos.front())/par->bin_size);
	int nsite=pos.size();
	int nad=par->popind[3];
	int shift=0;
	while(1<<shift < nbin)
		++shift;
	++shift;
	int num=1<<shift;

	double *n, *miss2, *w, *s, *ka;
	char *geno;
	bool *skip;
	try
	{
		n=new double[num];
		miss2=new  double[num];
		w=new double[nsite];
		s=new double[nsite];
		ka=new double[nsite];
		geno=new char[nad*nsite];
		skip=new bool[nsite];
	}
	catch (std::bad_alloc &e)
	{
		err_print("Memory allocation failure!");
	}
	memset(n,0,sizeof(double)*num);
	memset(miss2,0,sizeof(double)*num);
	memset(w,0,sizeof(double)*nsite);
	memset(s,0,sizeof(double)*nsite);
	memset(ka,0,sizeof(double)*nsite);
	memset(geno,0,sizeof(char)*nad*nsite);
	memset(skip,0,sizeof(bool)*nsite);

	ifstream fpi(par->pgeno[chr].c_str());
	string snp;
	int p(0);
	vector<int> & tpops(par->pops);
	for(int i=0;i<nsite;++i)
	{
		p=i;
		getline(fpi,snp);
		int nind=par->pops.size();
		int bin=floor((pos[i]-pos.front())/par->bin_size);
		double count1(0),nallele1(0),count2(0),nallele2(0);
		if(snp.size()!=nind)
		{
			cerr<<par->pgeno[chr]<<"\nLine:"<<i+1<<endl;
			cerr<<"There are "<<snp.size()<<" SNPs, while your dataset contains "<<nind<<" indiviudlas!"<<endl;
			err_print("Pls check your genotype files!");
		}
		for(int j=0;j<nind;++j)
		{
			int tsnp=snp[j]-'0';
			if(tpops[j]==1 && tsnp!=9)
			{
				++count1;
				nallele1+=tsnp;
			}
			else if(tpops[j]==2 && tsnp!=9)
			{
				++count2;
				nallele2+=tsnp;
			}
			else if(tpops[j]==3)
			{
				geno[p]=tsnp;
				p+=nsite;
				if(tsnp==9)
					++ka[i];
				else
					s[i]+=tsnp;
			}
		}
		if(count1==0 || count2==0 || nad-ka[i]<=1)
			skip[i]=true;
		else
		{
			w[i]=nallele1/count1/2-nallele2/count2/2;
			++n[bin];
			if(ka[i])
				++miss2[bin];
		}
	}
	getline(fpi,snp);
	if(fpi && snp!="")
	{
		cerr<<par->pgeno[chr]<<endl;
		cerr<<"There are "<<nsite<<" SNPs, while more markers in your genotype file!"<<endl;
		err_print("Pls check your genotype files!");
	}

	int nsample=(par->max_dis-par->min_dis)/par->bin_size;

	ref2(n, miss2, w, s, ka, geno, skip, wld, count, num, nsite, nsample);

	delete[]n;
	delete[]miss2;
	delete[]w;
	delete[]s;
	delete[]ka;
	delete[]geno;
	delete[]skip;
}

void WldChrDataCOM::calculation() throw()
/*
 * Calculation function for derived class WldChrDataCOM (Combined LD Model)
 *
 */
{
	int nbin=ceil((pos.back()-pos.front())/par->bin_size);
	int nsite=pos.size();
	int nad=par->popind[3];
	int r1=par->popind[1];
	int r2=par->popind[2];
	int shift=0;
	while(1<<shift < nbin)
		++shift;
	++shift;
	int num=1<<shift;

	double *n1, *n2, *n3, *n4, *miss2, *af1, *sx11, *sx21, *af2, *sx12, *sx22, *w, *wref2, *s, *ka;
	char *geno1, *geno2, *geno3;
	int *bin;
	bool *skip1, *skip2, *skip3;
	try
	{
		n1=new double[num];
		n2=new double[num];
		n3=new double[num];
		n4=new double[num];
		miss2=new double[num];
		af1=new double[nsite];
		af2=new double[nsite];
		sx11=new double[nsite];
		sx21=new double[nsite];
		sx12=new double[nsite];
		sx22=new double[nsite];
		w=new double[num];
		wref2=new double[nsite];
		s=new double[nsite];
		ka=new double[nsite];
		geno1=new char[r1*nsite];
		geno2=new char[r2*nsite];
		geno3=new char[nad*nsite];
		bin=new int[nsite];
		skip1=new bool[nsite];
		skip2=new bool[nsite];
		skip3=new bool[nsite];
	}
	catch (std::bad_alloc &e)
	{
		err_print("Memory allocation failure!");
	}
	memset(n1,0,sizeof(double)*num);
	memset(n2,0,sizeof(double)*num);
	memset(n3,0,sizeof(double)*num);
	memset(n4,0,sizeof(double)*num);
	memset(miss2,0,sizeof(double)*num);
	memset(af1,0,sizeof(double)*nsite);
	memset(af2,0,sizeof(double)*nsite);
	memset(sx11,0,sizeof(double)*nsite);
	memset(sx21,0,sizeof(double)*nsite);
	memset(sx12,0,sizeof(double)*nsite);
	memset(sx22,0,sizeof(double)*nsite);
	memset(w,0,sizeof(double)*num);
	memset(wref2,0,sizeof(double)*nsite);
	memset(s,0,sizeof(double)*nsite);
	memset(ka,0,sizeof(double)*nsite);
	memset(geno1,0,sizeof(char)*r1*nsite);
	memset(geno2,0,sizeof(char)*r2*nsite);
	memset(geno3,0,sizeof(char)*nad*nsite);
	memset(bin,0,sizeof(int)*nsite);
	memset(skip1,0,sizeof(bool)*nsite);
	memset(skip2,0,sizeof(bool)*nsite);
	memset(skip3,0,sizeof(bool)*nsite);

	ifstream fpi(par->pgeno[chr].c_str());
	string snp;
	int p1(0),p2(0),p3(0);
	int nind=par->pops.size();
	double adp1(0), adp2(0);
	vector<int> & tpops(par->pops);
	for(int i=0;i<nsite;++i)
	{
		p1=i;
		p2=i;
		p3=i;
		getline(fpi,snp);
		bin[i]=floor((pos[i]-pos.front())/par->bin_size);
		double count1(0),count2(0),tsx11(0),tsx12(0),tsx21(0),tsx22(0),ts(0),tka(0);
		bool missing2ref(false),missing1to2(false),missing2to1(false),missingw(false);
		if(snp.size()!=nind)
		{
			cerr<<par->pgeno[chr]<<"\nLine:"<<i+1<<endl;
			cerr<<"There are "<<snp.size()<<" SNPs, while your dataset contains "<<nind<<" indiviudlas!"<<endl;
			err_print("Pls check your genotype files!");
		}
		for(int j=0;j<nind;++j)
		{
			int tsnp=snp[j]-'0';
			if(tpops[j]==1)
			{
				geno1[p1]=tsnp;
				if(tsnp!=9)
				{
					++count1;
					tsx11+=tsnp;
					tsx21+=tsnp*tsnp;
				}
				else
					missing2to1=true;
				p1+=nsite;
			}
			else if(tpops[j]==2)
			{
				geno2[p2]=tsnp;
				if(tsnp!=9)
				{
					++count2;
					tsx12+=tsnp;
					tsx22+=tsnp*tsnp;
				}
				else
					missing1to2=true;
				p2+=nsite;
			}
			else if(tpops[j]==3)
			{
				geno3[p3]=tsnp;
				if(tsnp==9)
					++tka;
				else
					ts+=tsnp;
				p3+=nsite;
			}
		}
		if(nad-tka<=1)
			missing2ref=true;
		if(count1==0 || count2==0)
		{
			missing2ref=true;
			missing1to2=true;
			missing2to1=true;
			missingw=true;
		}
		int tbin=bin[i];
		double taf1=tsx11/count1/2;
		double taf2=tsx12/count2/2;
		if(missing2ref)
			skip1[i]=true;
		else
		{
			++n1[tbin];
			wref2[i]=taf1-taf2;
			s[i]=ts;
			ka[i]=tka;
			if(tka)
				++miss2[tbin];
		}
		if(missing2to1)
			skip2[i]=true;
		else
		{
			++n2[tbin];
			af2[i]=taf2;
			sx11[i]=tsx11;
			sx21[i]=tsx21;
		}
		if(missing1to2)
			skip3[i]=true;
		else
		{
			++n3[tbin];
			af1[i]=taf1;
			sx12[i]=tsx12;
			sx22[i]=tsx22;
		}
		if(!missingw && nad!=tka)
		{
			++n4[tbin];
			w[tbin]+=(taf1-taf2)*(taf1-taf2);
			adp1+=(ts/(nad-tka)/2-taf2)*(taf1-taf2);
			adp2+=(taf1-taf2)*(taf1-taf2);
		}
	}
	getline(fpi,snp);
	if(fpi && snp!="")
	{
		cerr<<par->pgeno[chr]<<endl;
		cerr<<"There are "<<nsite<<" SNPs, while more markers in your genotype file!"<<endl;
		err_print("Pls check your genotype files!");
	}

	int nsample=(par->max_dis-par->min_dis)/par->bin_size;
	*adp=adp1;
	*(adp+1)=adp2;

	ref2(n1, miss2, wref2, s, ka, geno3, skip1, exwld, excount, num, nsite, nsample);
	ref1(n2, af2, sx11, sx21, bin, geno1, skip2, exwld+nsample, excount+nsample, num, nsite, nsample,par->popind[1]);
	ref1(n3, af1, sx12, sx22, bin, geno2, skip3, exwld+nsample*2, excount+nsample*2, num, nsite, nsample,par->popind[2]);
	weight(n4, w, exwld+nsample*3, excount+nsample*3, num, nsample);

	delete[]n1;
	delete[]n2;
	delete[]n3;
	delete[]n4;
	delete[]miss2;
	delete[]af1;
	delete[]af2;
	delete[]sx11;
	delete[]sx21;
	delete[]sx12;
	delete[]sx22;
	delete[]w;
	delete[]wref2;
	delete[]s;
	delete[]ka;
	delete[]geno1;
	delete[]geno2;
	delete[]geno3;
	delete[]bin;
	delete[]skip1;
	delete[]skip2;
	delete[]skip3;
}




