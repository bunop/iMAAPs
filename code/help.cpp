/*
 * help.cpp
 *
 *  Created on: 2015.4.17
 *      Author: Yuan Kai
 */

# include <cstdlib>
# include "help.h"

using namespace std;

void print_help()
{
	cout<<SOFT_NAME<<"\tV."<<VER<<endl;
	cout<<endl;
	cout<<SOFT_NAME<<" <options> <parameters>"<<endl<<endl;
	cout<<"Options:"<<endl;
	cout<<"     -h/-H                                       Show this help."<<endl;
	cout<<"     -p <parameter file>                         Run "<<SOFT_NAME<<" with parameter files."<<endl;
	cout<<"     -help parameter/-help para/-P               Show the help of parameter file."<<endl;
	cout<<"     -help version/-V                            Show the version of "<<SOFT_NAME<<"."<<endl;
	cout<<endl;
	exit(0);
}

void print_para()
{
	cout<<SOFT_NAME<<"\tV."<<VER<<endl;
	cout<<endl;
	cout<<"Parameter help"<<endl;
	cout<<endl;
	cout<<"Format:  <option>:  <parameter>"<<endl;
	cout<<endl;
	cout<<"Basic setting:"<<endl;
	cout<<"     Inputfilelist   List of input files.                           Default: -"<<endl;
	cout<<"     Indfile         Individual notation.                           Default: -"<<endl;
	cout<<"     Admpop          Label for admixed population.                  Default: -"<<endl;
	cout<<"     Refpops         Labels for reference populations.              Default: -"<<endl;
	cout<<"     Timefile        Time file for testing the admixture signals.   Default: -"<<endl;
	cout<<"     Outdir          Output directory.                              Default: ./"<<endl;
	cout<<endl;
	cout<<"Advanced setting:"<<endl;
	cout<<"     Jackknife       Whether applied jack knife test.               Default: off"<<endl;
	cout<<"                     \"1\" \"on\"  \"ON\"  stands for open."<<endl;
	cout<<"                     \"0\" \"off\" \"OFF\" stands for close."<<endl;
	cout<<"     Mindis          Minimum distance of two WLD bins(in Morgan).   Defalut: 0.0002"<<endl;
	cout<<"     Maxdis          Maximum distance of two WLD bins(in Morgan).   Default: 0.3"<<endl;
	cout<<"     Binsize         Bin size of WLD bin(in Morgan)                 Default: 0.0002"<<endl;
	cout<<"     Iter            Number of iteration to the final fitting results."<<endl;
	cout<<"                                                                    Default: 100000"<<endl;
	cout<<"     Num_threads     Number of threads used globally.               Default: 1"<<endl;
	cout<<"     Num_threads_wld Number threads used during WLD calculation.    Default: equal to num_threads"<<endl;
	exit(0);
}

void print_version()
{
	cout<<SOFT_NAME<<"\tV."<<VER<<endl;
	cout<<endl;
	exit(0);
}
