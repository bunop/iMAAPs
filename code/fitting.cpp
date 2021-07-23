/*
 * fitting.cpp
 *
 *  Created on: 2014.12.27
 *      Author: Joe
 *
 */

#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_sort_vector.h>
#include<gsl/gsl_blas.h>
#include <stdlib.h>
# include <vector>
# include <iostream>
# include <fstream>
# include <math.h>

using namespace std;

int proj2splx(gsl_vector *y){
	int m, bget=0, i;
	double tmpsum=0, tmax;
	gsl_vector *s =static_cast<gsl_vector*>(gsl_vector_alloc(y->size));
	gsl_vector_memcpy(s, y);
	m= s->size;
	gsl_sort_vector(s);
	gsl_vector_reverse(s);
//	gsl_vector_fprintf(stderr, s, "%lf");
	for(i=0;i<(m-1);i++){
		tmpsum=tmpsum+gsl_vector_get(s, i);
		tmax=(tmpsum-1)/(i+1);
//	fprintf(stderr, "tmax %lf\n", tmax);
		if(tmax>=gsl_vector_get(s, i+1)){
			bget = 1;
			break;
		}
	}
	if(!bget){
		tmax=(tmpsum + gsl_vector_get(s, m-1)-1)/m;
	}
//	gsl_vector_fprintf(stderr, y, "%lf");
//	fprintf(stderr, "tmax %lf\n", tmax);
	gsl_vector_add_constant(y, -tmax);
	for(i=0;i<m;i++){
		if(gsl_vector_get(y, i)<0)gsl_vector_set(y,i,0);
	}
	gsl_vector_free(s);
	return 0;
}


int proj2l1(gsl_vector *y, double C){
	int i,m=y->size;
	for(i=0;i<m;i++){if(gsl_vector_get(y, i)<0)gsl_vector_set(y,i,0);}
//	gsl_vector_fprintf(stderr, y, "0-%lf");
	if(gsl_blas_dasum(y)>C){
		gsl_vector_scale(y, 1/C);
		proj2splx(y);
		gsl_vector_scale(y, C);
	}
//	gsl_vector_fprintf(stderr, y, "1-%lf");
	return 0;

}

int apg(gsl_matrix *w_m,gsl_matrix *fit_m, gsl_matrix *A, gsl_vector *b, gsl_vector *w, double  eta0, int burnin, double C, int spl_num,int step){
	double t, eta, t_old;
	int maxiter;
	int jj, ii;
	maxiter= burnin+spl_num*step;
//	cout<<"fitting threshold: "<<C<<endl;
//	fprintf(stderr, "fitting threshold: %lf\n", C);
	proj2l1(w, C);
	gsl_vector *z =static_cast<gsl_vector*>(gsl_vector_alloc(w->size));
	gsl_vector *w_old =static_cast<gsl_vector*>(gsl_vector_alloc(w->size));
	gsl_vector *tmp1 =static_cast<gsl_vector*>(gsl_vector_alloc(w->size));
	gsl_vector *res =static_cast<gsl_vector*>(gsl_vector_alloc(b->size));
	gsl_vector *tmp2 =static_cast<gsl_vector*>(gsl_vector_alloc(b->size));
	gsl_vector *g =static_cast<gsl_vector*>(gsl_vector_alloc(w->size));
	gsl_vector_memcpy(z, w);
//	gsl_matrix *w_m=gsl_matrix_alloc(spl_num, w->size);
	gsl_vector *f =static_cast<gsl_vector*>(gsl_vector_alloc(maxiter));
	t = 1;
//	tol = 0.00000000001;
	eta = eta0;
	jj=0;
//	gsl_vector_scale(b, 1/C);
	for(ii=0;ii<maxiter;ii++){
//		for(i=0;i<w->size;i++)fprintf(p, "%lf ", gsl_vector_get(w,i));
//		fprintf(p, "\n");
		gsl_vector_memcpy(w_old, w);
		gsl_vector_memcpy(res, b);
		gsl_blas_dgemv(CblasNoTrans, 1,A, z, -1, res);
	//	res = A%*%z - b;
		gsl_vector_set(f,ii,gsl_blas_dnrm2(res));
	//      	fprintf(p, "%.15lf\n", gsl_vector_get(f, ii));

	//	f[ii] = norm(res, type="F")^2 / 2;
		gsl_blas_dgemv(CblasTrans, 1,A, res, 0, g);
	//	g = t(A)%*%res;
		gsl_vector_memcpy(w, z);
		gsl_blas_daxpy(-eta, g, w);
		proj2l1(w, C);

/*		while(1){
			w = proj2l1(z - eta*g, C);
			f_new = norm(A%*%w-b, type="F")^2 / 2;
			f_up = f[ii] + t(w-z)%*%g + 1/2/eta*norm(w-z,type="F")^2;
			if(f_new <= f_up + tol){
				eta = eta*2;
				break;}
			else{ eta[eta/1.5<eta0]<-eta0[eta/1.5<eta0];eta[eta/1.5>=eta0]<-eta[eta/1.5>=eta0]/1.5; }
		}
*/
		t_old = t;
		t = (1 + sqrt(1+4*t*t)) / 2;
		gsl_vector_memcpy(tmp1, w);
		gsl_vector_sub(tmp1, w_old);
		gsl_vector_memcpy(z, w);
		gsl_blas_daxpy((t_old -1)/t, tmp1, z);
//		z = w + (t_old -1) / t * (w - w_old);
		if(ii==(burnin+jj*step)){
			gsl_matrix_set_row(w_m, jj, w);
			gsl_blas_dgemv(CblasNoTrans, 1,A, w, 0, tmp2);
			gsl_matrix_set_row(fit_m, jj, tmp2);
//			fprintf(stderr, "sample: %d, total: %d\n",jj, spl_num);
			jj=jj+1;
		}
	}
	if(ii==maxiter){
//		cout<<"max iteration reached!!!"<<endl;
//		fprintf(stderr, "max iteration reached!!!\n");
	}

	/*debug~~~~~~~
	FILE *p=fopen("tmp", "w");

	 //     gsl_matrix_fprintf(p,z , "%.15lf");
	//      fprintf(stderr, "%d %d\n", A->size1, A->size2);
	//      gsl_vector_fprintf(p,z, "%.15lf");
//	fclose(p);
	*/


	gsl_vector_free(z);
	gsl_vector_free(w_old);
	gsl_vector_free(tmp1);
	gsl_vector_free(tmp2);
	gsl_vector_free(res);
	gsl_vector_free(g);
	gsl_vector_free(f);

	return 0;
}

//int main(int arg, char *argv[]){
void fitting(double *tmp_y, double *tmp_d, double *tmp_fit, const int & len_ld,const vector<int>& time_pt, double *time_out,const int& len_time,const int& n_burnin,const double & C_thd){


	int i,j;
//	double dist,d, LD,LD2, LD3;
//	char if_tmpt[256], of_nm[256],ofname[256], tmps[2000];
//	FILE *ifp, *ofp;

	double eta0;
	int burnin;
	double C;
	int spl_num=20;
	int step=2;
//	double C_thd;

	burnin=n_burnin;


//	C_thd=CC; //default value;

//	cout<<"There are "<<len_time<<" time points to test and "<<len_ld<<" samples points to fit"<<endl;
//	fprintf(stderr, "There are %d time points to test and %d samples points to fit\n", len_time, len_ld);

//	fprintf(stderr, "data reading finish!\nfitting!\n");
	//initiation of the coefficient of the exponential functions
	gsl_matrix *A =static_cast<gsl_matrix*>(gsl_matrix_calloc(len_ld, len_time));
	gsl_vector *w =static_cast<gsl_vector*>(gsl_vector_alloc(len_time));
	gsl_vector *y =static_cast<gsl_vector*>(gsl_vector_alloc(len_ld));
	gsl_matrix *w_m =static_cast<gsl_matrix*>(gsl_matrix_alloc(spl_num, w->size));
	gsl_matrix *fit_m =static_cast<gsl_matrix*>(gsl_matrix_alloc(spl_num, len_ld));
	//A<-t(exp(-(N1 %o% d)));
	double tmpd, tmp_sum=0;
	gsl_matrix_set_zero(w_m);
	for(i=0;i<len_ld;i++){
		gsl_vector_set(y,i,tmp_y[i]);
		//	gsl_vector_set(y,i,exp(-100*tmp_d[i]));
		for(j=0;j<len_time;j++){
			gsl_vector_set(w,j,drand48());
			tmpd=exp(-(tmp_d[i]/100*time_pt[j]));
//			tmpd=pow(1-tmp_d[i]/100, time_pt[j]);
			gsl_matrix_set(A,i, j,tmpd);
			tmp_sum=tmp_sum+tmpd*tmpd;
		}
	}
	//eta<-1/(norm(A, type="2"))^2;
	eta0=1/tmp_sum;
	//C
	C=gsl_vector_max(y)*C_thd;
	//cout<<C<<"~~~~~~~~~~~~~"<<endl;
	//fitting, save the results in w_m.
	/*debug~~~~~~~

	  FILE *p=fopen("tmp1", "w");
	//	gsl_matrix_fprintf(p,A , "%.15lf");
	//	fprintf(stderr, "%d %d\n", A->size1, A->size2);
	gsl_vector_fprintf(p,y, "%.15lf");
	fclose(p);
	//	return 0;
	*/
	apg(w_m, fit_m, A, y, w,eta0,burnin,C, spl_num,step);
	//sprintf(ofname,"%s_w", of_nm);
	//ofp=fopen(ofname, "w");
	int tmpw=w->size;
	for(j=0;j< tmpw; j++){
		//fprintf(ofp, "%d", time_pt[j]);
		for(i=0;i<spl_num;i++){
		//	fprintf(ofp, " %lf",gsl_matrix_get(w_m,i,j));
			time_out[j]=gsl_matrix_get(w_m,i,j);
		}
	//	fprintf(ofp, "\n");
	}
	//fclose(ofp);

	//sprintf(ofname,"%s_fit", of_nm);
	//ofp=fopen(ofname, "w");
	for(j=0;j< len_ld; j++){
		//fprintf(ofp, "%lf", tmp_d[j]);
		for(i=0;i<spl_num;i++){
			//fprintf(ofp, " %.10lf",gsl_matrix_get(fit_m,i,j));
			tmp_fit[j]=gsl_matrix_get(fit_m,i,j);
		}
	//	fprintf(ofp, "\n");
	}
	//fclose(ofp);



	gsl_vector_free(w);
	gsl_vector_free(y);
	gsl_matrix_free(A);
	gsl_matrix_free(w_m);
	gsl_matrix_free(fit_m);
}




