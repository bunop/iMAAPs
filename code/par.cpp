/*
 * par.cpp
 *
 *  Created on: 2015年3月31日
 *      Author: yuankai_local
 */

# include "par.h"
# include "err.h"
# include "tools.h"
# include <sstream>
# include <map>
# include <cmath>

using namespace std;

void WldPar::read(const std::string &path)
/*
 * Read parameters from file.
 *
 * path    string to store the path of parameter file
 *
 */
{
	if(!file_check(path))
		err_print("Cannot open parameter file!\n"+path);
	ifstream fpi(path.c_str());
	string line;
	while(getline(fpi,line))
		get_par(line);
	if(min_dis<0)
		min_dis=0.0002;
	if(max_dis<0)
		max_dis=0.3;
	if(bin_size<0)
		bin_size=0.0002;
	if(burnin<0)
		burnin=100000;
	if(n_thread<0)
		n_thread=1;
	if(n_thread_wld<0)
		n_thread_wld=n_thread;
	if(mode==Null)
		mode=COM;
	if(pout=="")
		pout=".";
	if(c_thd<0)
		c_thd=1.0/exp(-100.0*min_dis);
	pwld=pout+"/"+adname[0]+".rawld";
	pad=pout+"/"+adname[0]+".ad";
	//Set default value.

	pop_in();
	//Input individual information
}

void WldPar::check() const
/*
 * Check all of these parameters
 */
{
	if(!file_check(plist))
		err_print("Cannot open the input files list!\nInput list: "+plist);
	if(!file_check(ptime))
		err_print("Cannot open the generation file!Generation file: \n"+ptime);
	if(!file_check(pind))
		err_print("Cannot open the individual file!\nIndividual file: "+pind);
	if(pwld=="" || !file_open(pwld))
		err_print("Cannot open the Weighted LD file!\nWLD file: "+pwld);
	if(pad=="" || !file_open(pad))
		err_print("Cannot open the admixed file!\nAdmixed file: "+pad);
	if(pgeno.size()==0)
		err_print("Number of Genotype files cannot be 0!");
	if(psnp.size()==0)
		err_print("Number of SNP files cannot be 0!");
	if(pgeno.size()!=psnp.size())
		err_print("Number genotype files and that of SNP files should be equal!");
	if(pgeno.size()==1 && jackknife)
		err_print("Jackknife test cannot applied with only one chromosome!");
	if(min_dis<=0)
		err_print("Minimum distance should be larger than 0!");
	if(max_dis<=0)
		err_print("Maximum distance should be larger than 0!");
	if(min_dis>=max_dis)
		err_print("Maximum distance should be larger than minimum distance!");
	if(bin_size<=0)
		err_print("Bin size should be larger than 0!");
	if(burnin<=0)
		err_print("Burn in number should be larger than 0!");
	if(n_thread<=0)
		err_print("Number of thread should be larger than 0!");
	if(mode==W)
	{
		if(refname.size()!=2)
			err_print("Number of reference populations under Weighted Model should be 2!");
		if(popind[1]==0 || popind[2]==0)
			err_print("Number of individuals in reference population should not be 0!");
	}
	if(mode==W2)
	{
		if(refname.size()!=2)
			err_print("Number of reference populations under Weighted Square Model should be 2!");
		if(popind[1]==0 || popind[2]==0)
			err_print("Number of individuals in reference population should not be 0!");
	}
	if(mode==REF1)
	{
		if(refname.size()!=1)
			err_print("Number of reference populations under 1 Reference Model should be 1!");
		if(adname.size()!=1)
			err_print("Number of admixed populations under 1 Reference Model should be 1!");
		if(popind[1]==0)
			err_print("Number of individuals in reference population should not be 0!");
		if(popind[2]<=3)
			err_print("Number of individuals in admixed population should be larger than 3!");
	}
	if(mode==REF2)
	{
		if(refname.size()!=2)
			err_print("Number of reference populations under 2 Reference Model should be 2!");
		if(adname.size()!=1)
			err_print("Number of admixed populations under 2 Reference Model should be 1!");
		if(popind[1]==0 || popind[2]==0)
			err_print("Number of individuals in reference population should not be 0!");
		if(popind[3]<=1)
			err_print("Number of individuals in admixed population should be larger than 1!");
	}
	if(mode==COM)
	{
		if(refname.size()!=2)
			err_print("Number of reference populations under Combined LD Model should be 2!");
		if(adname.size()!=1)
			err_print("Number of admixed populations under Combined LD Model should be 1!");
		if(popind[1]==0 || popind[2]==0)
			err_print("Number of individuals in reference populations should not be 0!");
		if(popind[3]<=1)
			err_print("Number of individuals in admixed population should be larger than 1!");
		if(popind[1]<=3 || popind[2]<=3)
			err_print("Number of individuals in reference population should be larger than 3!");
	}
}

void WldPar::print() const
/*
 * Print parameters on screen
 */
{
	print_line();
	cout<<endl;
	cout<<"Input files:"<<endl;
	cout<<"    Input list:      "<<plist<<endl;
	int nchr=pgeno.size();
	for(int i=0;i<nchr;++i)
		cout<<"                     "<<pgeno[i]<<'\t'<<psnp[i]<<endl;
	cout<<"    Individual file: "<<pind<<endl;
	cout<<"    Generation file: "<<ptime<<endl;
	cout<<endl;
	cout<<"Output files:"<<endl;
	cout<<"    WLD file:        "<<pwld<<endl;
	cout<<"    Admixed file:    "<<pad<<endl;
	cout<<endl;
	cout<<"Reference populations:"<<endl;
	for(size_t i=0;i<refname.size();++i)
		cout<<"    "<<i+1<<'.'<<refname[i]<<endl;
	cout<<endl;
	cout<<"Admixed populations:"<<endl;
	for(size_t i=0;i<adname.size();++i)
		cout<<"    "<<i+1<<'.'<<adname[i]<<endl;
	cout<<endl;
	cout<<"WLD parameters:"<<endl;
	cout<<"    Minimum distance:"<<min_dis<<endl;
	cout<<"    Maximum distance:"<<max_dis<<endl;
	cout<<"    Bin size:        "<<bin_size<<endl;
	cout<<endl;
	cout<<"Fitting parameters:"<<endl;
	cout<<"    Iteration:       "<<burnin<<endl;
	cout<<"    C threshold:     "<<c_thd<<endl;
	cout<<endl;
	cout<<"Jackknife:  "<< (jackknife ? "on" : "off") <<endl;
	cout<<endl;
	cout<<"Mode:       ";
	if(mode==W)			cout<<"Weighted Model"<<endl;
	else if(mode==W2)	cout<<"Weighted Square Model"<<endl;
	else if(mode==REF1)	cout<<"1 Reference Model"<<endl;
	else if(mode==REF2)	cout<<"2 Reference Model"<<endl;
	else if(mode==COM)	cout<<"Combined LD Model"<<endl;
	cout<<endl;
}

void WldPar::get_par(const std::string & line)
/*
 * Function to abstract information parameter string.
 *
 * line     A string contains parameter label and its value
 *
 */
{
	istringstream cl(line);
	string lab;
	cl>>lab;
	string tmp;
	if(lab=="Inputfilelist:")
	{
		if(plist!="")
			err_print("Duplicate list path input!\nYou have already input list path "+plist+"!");
		cl>>plist;
		if(!file_check(plist))
			err_print("Cannot open list file!\n"+plist+" doesn't exist!");
		cl>>tmp;
		if(cl)
			err_print("List path parameter is too long!\n"+line);

		list_in();	//Input paths information of list file
	}
	else if(lab=="Timefile:")
	{
		if(ptime!="")
			err_print("Duplicate generation file path input!\nYou have already input path "+ptime+"!");
		cl>>ptime;
		if(!file_check(ptime))
			err_print("Cannot open generation file!\n"+ptime+" doesn't exist!");
		cl>>tmp;
		if(cl)
			err_print("Generation file path parameter is too long!\n"+line);

		time_in();	//Input generation information
	}
	else if(lab=="Outdir:")
	{
		if(pout!="")
			err_print("Duplicate output directory path input!\nYou have already input path "+pout+"!");
		cl>>pout;
		cl>>tmp;
		if(cl)
			err_print("Output directory path is too long!\n"+line);
	}
	else if(lab=="Indfile:")
	{
		if(pind!="")
			err_print("Duplicate individual file path input!\nYou have already input path "+pind+"!");
		cl>>pind;
		if(!file_check(pind))
			err_print("Cannot open individual file!\n"+pind+" doesn't exist!");
		cl>>tmp;
		if(cl)
			err_print("Individual file path parameter is too long!\n"+line);
	}
	else if(lab=="Admpop:")
	{
		if(adname.size())
			err_print("Duplicate admixed population labels input!");
		string subline;
		cl>>subline;
		split_str(adname,subline);
		cl>>tmp;
		if(cl)
			err_print("Admixed population labels parameter is too long!\n"+line);
		if(adname.size()==0)
			err_print("No admixed population label intput!");
	}
	else if(lab=="Refpops:")
	{
		if(refname.size())
			err_print("Duplicate reference population labels input!");
		string subline;
		cl>>subline;
		split_str(refname,subline);
		cl>>tmp;
		if(cl)
			err_print("Reference population labels parameter is too long!\n"+line);
	}
	else if(lab=="Mindis:")
	{
		if(min_dis>=0)
		{
			char fnum[64];
			sprintf(fnum,"%f",min_dis);
			string tfnum=fnum;
			err_print("Duplicate fitting minimum distance input!\nYou have already input "+tfnum);
		}
		cl>>min_dis;
		if(min_dis<=0)
			err_print("Minimum distance cannot be negative!");
		cl>>tmp;
		if(cl)
			err_print("Minimum distance parameter is too long!\n"+line);
	}
	else if(lab=="Maxdis:")
	{
		if(max_dis>=0)
		{
			char fnum[64];
			sprintf(fnum,"%f",max_dis);
			string tfnum=fnum;
			err_print("Duplicate fitting maximum distance input!\nYou have already input "+tfnum);
		}
		cl>>max_dis;
		if(max_dis<=0)
			err_print("Maximum distance cannot be negative!");
		cl>>tmp;
		if(cl)
			err_print("Maximum distance parameter is too long!\n"+line);
	}
	else if(lab=="Binsize:")
	{
		if(bin_size>=0)
		{
			char fnum[64];
			sprintf(fnum,"%f",bin_size);
			string tfnum=fnum;
			err_print("Duplicate WLD bin size input!\nYou have already input "+tfnum);
		}
		cl>>bin_size;
		if(bin_size<=0)
			err_print("Bin size should be positive!");
		cl>>tmp;
		if(cl)
			err_print("Bin size parameter is too long!\n"+line);
	}
	else if(lab=="Iter:")
	{
		if(burnin>=0)
		{
			char fnum[64];
			sprintf(fnum,"%d",burnin);
			string tfnum=fnum;
			err_print("Duplicate fitting burn in input!\nYou have already input "+tfnum);
		}
		cl>>burnin;
		if(burnin<=0)
			err_print("Fitting burn in should be positive!");
		cl>>tmp;
		if(cl)
			err_print("Fitting burn in parameter is too long!\n"+line);
	}
	else if(lab=="Num_threads:")
	{
		if(n_thread>=0)
		{
			char fnum[64];
			sprintf(fnum,"%d",n_thread);
			string tfnum=fnum;
			err_print("Duplicate number of thread input!\nYou have already input "+tfnum);
		}
		cl>>n_thread;
		if(n_thread<=0)
			err_print("Number of thread should be positive!");
		cl>>tmp;
		if(cl)
			err_print("Number of thread parameter is too long!\n"+line);
	}
	else if(lab=="Num_threads_wld:")
	{
		if(n_thread_wld>=0)
		{
			char fnum[64];
			sprintf(fnum,"%d",n_thread_wld);
			string tfnum=fnum;
			err_print("Duplicate input of number of thread for WLD calculation!\nYou have already input "+tfnum);
		}
		cl>>n_thread_wld;
		if(n_thread_wld<=0)
			err_print("Number of thread should be positive!");
		cl>>tmp;
		if(cl)
			err_print("Number of thread parameter is too long!\n"+line);
	}
	else if(lab=="Jackknife:")
	{
		if(jackknife)
		{
			err_print("Duplicate jackknife parameter input!");
		}
		cl>>tmp;
		if(tmp=="1" || tmp=="on" || tmp=="ON")
			jackknife=true;
		else if(tmp=="0" || tmp=="off" || tmp=="OFF")
			jackknife=false;
		else
			err_print("Cannot identify jackknife parameter "+tmp);
		cl>>tmp;
		if(cl)
			err_print("Jackknife parameter is too long!\n"+line);
	}
	else
		err_print("Cannot identify parameter "+line+"!");
}

void WldPar::split_str(vector< string > & names, const string & line)
/*
 * Split string into parts by semicolon
 */
{
	names.clear();
	size_t pre(0);
	for(size_t i=0;i<line.length();++i)
	{
		if(line[i]==';')
		{
			if(i!=pre)
				names.push_back(line.substr(pre,i-pre));
			pre=i+1;
		}
	}
	if(pre!=line.length())
		names.push_back(line.substr(pre));
}

void WldPar::list_in()
/*
 * Read paths of input files.
 */
{
	if(pgeno.size() || psnp.size() || chr.size())
		err_print("Duplicate input file paths input!");
	ifstream fpi(plist.c_str());
	string pg,ps;
	int ch;
	pgeno.clear();
	pgeno.push_back("");
	psnp.clear();
	psnp.push_back("");
	chr.clear();
	while(fpi>>ch>>pg>>ps)
	{
		chr.push_back(ch);
		pgeno.push_back(pg);
		psnp.push_back(ps);
		if(!file_check(pg))
			err_print("Cannot open genotype file "+pg);
		if(!file_check(ps))
			err_print("Cannot open snp file "+ps);
	}
}

void WldPar::pop_in()
/*
 * Read individual information
 */
{
	if(popind.size())
		err_print("Duplicate individual information input!");
	if(!refname.size() || !adname.size())
		err_print("Population labels of admixed or reference population cannot be 0!");
	ifstream fp(pind.c_str());
	string id,pop;
	std::map <std::string , int > mpop;
	mpop.clear();
	int p(0);
	for(size_t i=0;i<refname.size();++i)
		mpop.insert(make_pair(refname[i],++p));
	for(size_t i=0;i<adname.size();++i)
		mpop.insert(make_pair(adname[i],++p));
	popind.clear();
	popind.resize(p+1);
	for(int i=0;i<p;++i)
		popind[i]=0;
	pops.clear();
	string tmp;
	while(fp>>id>>tmp>>pop)
	{
		pops.push_back(mpop[pop]);
		++popind[mpop[pop]];
	}
}

void WldPar::time_in()
/*
 * Read generation file information
 */
{
	if(time.size())
		err_print("Duplicate generation information input!");
	ifstream fp(ptime.c_str());
	int nt;
	time.clear();
	while(fp>>nt)
		time.push_back(nt);
}
