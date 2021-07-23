/*
 * tools.cpp
 *
 *  Created on: 2015.3.31
 *      Author: Yuan Kai
 *
 *  This file include lots of tool function used in this soft.
 */

# include "tools.h"

void print_line()
/*
 * Print a line on screen
 */
{
	std::cout<<"------------------------------------------------------------------------------------"<<std::endl;
}

bool file_check( const std::string & path )
/*
 * Check whether file could be open as read only
 */
{
	std::ifstream fp(path.c_str());
	return fp;
}

bool file_open( const std::string & path )
/*
 * Check whether file could be open as write
 */
{
	std::ofstream fp(path.c_str());
	return fp;
}

double timer::gettime()
//Get time
{
	struct timeval now;
	gettimeofday(&now, NULL);
	return now.tv_sec-pre.tv_sec+(now.tv_usec-pre.tv_usec)*1e-6;
}

void timer::set()
//Set timer
{
	gettimeofday(&pre, NULL);
	mem=0;
}

void timer::add_set()
//Set add time
{
	gettimeofday(&add_pre, NULL);
}

void timer::add()
//Add time
{
	struct timeval now;
	gettimeofday(&now, NULL);
	mem+=now.tv_sec-add_pre.tv_sec+(now.tv_usec-add_pre.tv_usec)*1e-6;
}

double timer::add_get()
//Get add time
{
	return mem;
}
