/*
 * fitting.h
 *
 *  Created on: 2014.12.27
 *      Author: Joe
 *
 */

#ifndef FITTING_H_
#define FITTING_H_

# include <iostream>
# include <fstream>
# include <vector>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_sort_vector.h>
#include<gsl/gsl_blas.h>

int proj2splx(gsl_vector *y);

int proj2l1(gsl_vector *y, double C);

int apg(gsl_matrix *w_m,gsl_matrix *fit_m, gsl_matrix *A, gsl_vector *b, gsl_vector *w, double  eta0, size_t burnin, double C, size_t spl_num,size_t step);

void fitting(double *tmp_y, double *tmp_d, double *tmp_fit, const int & len_ld,const std::vector<int>& time_pt, double *time_out,const int& len_time,const int& n_burnin,const double& C_thd);


#endif /* FITTING_H_ */
