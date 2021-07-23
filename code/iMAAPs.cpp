/*
 * iMAAPs.cpp
 *
 *  Created on: 2015.4.15
 *      Author: Yuan Kai
 */


# include <iostream>
# include "data.h"
# include "fftw3.h"
# include "tools.h"
# include "help.h"

using namespace std;

int main(int argc,char **argv)
{
	string ppar("");
	for(int i=1;i<argc;++i)
	{
		string targv=argv[i];
		if(targv=="-p" && i!=argc-1)
			ppar=argv[++i];
		else if(targv=="-h" || targv=="-H")
			print_help();
		else if(targv=="-help")
		{
			if(i==argc-1)
				break;
			targv=argv[++i];
			if(targv=="parameter" || targv=="para")
				print_para();
			else if(targv=="version")
				print_version();
		}
		else if(targv=="-P")
			print_para();
		else if(targv=="-V")
			print_version();
		else
		{
			cerr<<"Error parameter input"<<endl;
			print_para();
		}
	}
	if(ppar=="")
	{
		print_help();
	}
	timer t;
	t.set();

	WldPar p;
	p.read(ppar);
	p.check();
	p.print();

	WldData data(p);
	data();

	cout<<"\nTotal running time is "<<t.gettime()<<" s."<<endl;
}

