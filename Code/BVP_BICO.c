#include "BVP_header.h"

int BVP_BICO(void *params, double A[], double b[], double x[], double err[])
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int N = my_ptr->Nparam;
	const int imax = my_ptr->imaxparam;
	const double tol = my_ptr->LStolparam;
	int i,j,k;
	double temp;
	double temp2;
	double *r0hat;
	double *r0;
	double *ri;
	double *p0;
	double *pi;
	double *v0;
	double *vi;
	double *x0;
	double *xi;
	double *s;
	double *t;
	double *h;
	double rho0, rhoi;
	double omega0, omegai;
	double alpha, beta;
	double b_norm;
	double r_norm;
	
	r0hat = calloc(N, sizeof(double));
	if(!r0hat){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	r0 = calloc(N, sizeof(double));
	if(!r0){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	ri = calloc(N, sizeof(double));
	if(!ri){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	p0 = calloc(N, sizeof(double));
	if(!p0){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	pi = calloc(N, sizeof(double));
	if(!pi){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	v0 = calloc(N, sizeof(double));
	if(!v0){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	vi = calloc(N, sizeof(double));
	if(!vi){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	x0 = calloc(N, sizeof(double));
	if(!x0){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	xi = calloc(N, sizeof(double));
	if(!xi){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	s = calloc(N, sizeof(double));
	if(!s){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	t = calloc(N, sizeof(double));
	if(!t){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	h = calloc(N, sizeof(double));
	if(!h){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	
	b_norm = 0.0;
	r_norm = 0.0;

	//0.
	for (j = 0; j < N; ++j)
	{
		x0[j] = x[j];
	}
	//1. & 2.
	for (j = 0; j < N; ++j)
	{
		temp = 0.0;
		for (k = 0; k < N; ++k)
		{
			temp -= A[j*N+k]*x0[k];
		}
		b_norm += b[j]*b[j];
		r0[j] = b[j] + temp;
		r0hat[j] = r0[j];
		r_norm += r0[j]*r0[j];
	}
	b_norm = sqrt(b_norm);
	r_norm = sqrt(r_norm);
	err[0] = r_norm/b_norm;
	//3.
	rho0 = 1.0;
	alpha = 1.0;
	omega0 = 1.0;
	for (j = 0; j < N; ++j)
	{
		v0[j] = 0.0;
		p0[j] = 0.0;
	}
	for (i = 1; i < imax; ++i)
	{
		//1
		temp = 0.0;
		for (j = 0; j < N; ++j)
		{
			temp += r0hat[j]*r0[j];
		}
		rhoi = temp;
		//2
		beta = rhoi/rho0*alpha/omega0;
		//3
		for (j = 0; j < N; ++j)
		{
			pi[j] = r0[j] + beta*(p0[j] - omega0*v0[j]);
		}
		//4
		for (j = 0; j < N; ++j)
		{
			vi[j] = 0.0;
			for (k = 0; k < N; ++k)
			{
				vi[j] += A[j*N+k]*pi[k];
			}
		}
		//5
		temp = 0.0;
		for (j = 0; j < N; ++j)
		{
			temp += r0hat[j]*vi[j];
		}
		alpha = rhoi/temp;
		//6
		for (j = 0; j < N; ++j)
		{
			h[j] = x0[j] + alpha*pi[j];
		}
		//7 Check convergence of h
		r_norm = 0.0;
		for (j = 0; j < N; j++)
		{
			temp = 0.0;
			for (k = 0; k < N; k++)
			{
				temp += A[j*N+k]*h[k];
			}
			temp -= b[j];
			r_norm += temp*temp;
		}
		r_norm = sqrt(r_norm);
		//printf("error for iteration %03i = %22.15e.\n",i,r_norm);
		if (r_norm < tol)
		{
			printf("Success: error at iteration %04i = %11.4e < %.2e.\n",i,r_norm,tol);
			err[i] = r_norm;
			for (j = 0; j < N; ++j)
			{
				xi[j] = h[j];
			}
			break;
		}
		//8
		for (j = 0; j < N; ++j)
		{
			s[j] = r0[j] - alpha*vi[j];
		}
		//9
		for (j = 0; j < N; ++j)
		{
			t[j] = 0.0;
			for (k = 0; k < N; ++k)
			{
				t[j] += A[j*N+k]*s[k];
			}
		}
		//10
		temp = 0.0;
		temp2 = 0.0;
		for (j = 0; j < N; ++j)
		{
			temp += t[j]*s[j];
			temp2 += t[j]*t[j];
		}
		omegai = temp/temp2;
		//11
		for (j = 0; j < N; ++j)
		{
			xi[j] = h[j] + omegai*s[j];
		}
		//12 Check convergence of xi
		r_norm = 0.0;
		for (j = 0; j < N; j++)
		{
			temp = 0.0;
			for (k = 0; k < N; k++)
			{
				temp += A[j*N+k]*xi[k];
			}
			temp -= b[j];
			r_norm += temp*temp;
		}
		r_norm = sqrt(r_norm);
		//printf("error for iteration %03i = %22.15e.\n",i,r_norm);
		if (r_norm < tol)
		{
			printf("Success: error at iteration %04i = %11.4e < %.2e.\n",i,r_norm,tol);
			err[i] = r_norm;
			break;
		}
		//13
		for (j = 0; j < N; ++j)
		{
			ri[j] = s[j] - omegai*t[j];
		}
		//reset values
		for (j = 0; j < N; ++j)
		{
			r0[j] = ri[j];
			rho0 = rhoi;
			if (isnan(rho0))
			{
				printf("ERROR in BVP_BICO. rhoi is NaN. Last error = %11.4e.\n",b_norm);
				return -1;
			}
			p0[j] = pi[j];
			v0[j] = vi[j];
			x0[j] = xi[j];
		}
		omega0 = omegai;
		err[i] = r_norm;
	}
	if (i == imax)
	{
		printf("BICO FAILURE. Error after %i iterations = %23.16e.\n",i,err[i-1]);
		return -1;
	}
	for (j = 0; j < N; ++j)
	{
		x[j] = xi[j];
	}
	
	//Cleanup
	free(r0hat);
	free(r0);
	free(ri);
	free(p0);
	free(pi);
	free(v0);
	free(vi);
	free(x0);
	free(xi);
	free(s);
	free(t);
	free(h);
	
	return 0;
}