/*
 * tools.h
 *
 *  Created on: 2015.3.31
 *      Author: Yuan Kai
 *
 *  This file include lots of tool function used in this soft.
 */

#ifndef TOOLS_H_
#define TOOLS_H_

# include <iostream>
# include <fstream>
# include <sys/time.h>

void print_line();

bool file_check( const std::string & path );

bool file_open( const std::string & path );

class timer
{
public:
	struct timeval pre,add_pre;
	double mem;
	double gettime();
	void set();
	void add_set();
	void add();
	double add_get();
};

#endif /* TOOLS_H_ */
