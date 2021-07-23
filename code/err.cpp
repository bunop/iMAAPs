/*
 * err.cpp
 *
 *  Created on: 2015.3.31
 *      Author: Yuan Kai
 *
 */

# include "err.h"

void err_print( const std::string & err )
/*
 * Function used to print error information on the screen.
 *
 */
{
	std::cerr<<err<<std::endl;
	exit(1);
}
