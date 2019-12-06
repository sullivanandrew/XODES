#include "BVP_header.h"

int BVP_LUdcmp(void *params, double A[], double b[], double x[], double *Det)
{
	struct param_type *my_ptr;
	my_ptr = params;
	int N = my_ptr->Nparam;
	int imax = my_ptr->imaxparam;
	double tol = my_ptr->LStolparam;
	int i,j,k;
	double *y;
	double *al;
	double *be;
	double temp;
	double Dettemp = 1.0;

	al = calloc(N*N, sizeof(double));
	if(!al){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	be = calloc(N*N, sizeof(double));
	if(!be){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	y = calloc(N, sizeof(double));
	if(!y){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }


	for (i = 0; i < N; i++)
	{
		al[i*N+i] = 1.0;
	}
	for (j = 0; j < N; j++)
	{
		for (i = 0; i < j+1; i++)
		{
			be[i*N+j] = A[i*N+j];
			for (k = 0; k < i; k++)
			{
				be[i*N+j] = be[i*N+j] - al[i*N+k]*be[k*N+j];
				//printf("al[%2i,%2i]*be[%2i,%2i] = %7.2f\n",i,k,k,j,al[i*N+k]*be[k*N+j]);
			}
			//printf("be[%2i,%2i] = %7.2f\n",i,j,be[i*N+j]);
		}
		for (i = j+1; i < N; i++)
		{
			al[i*N+j] = A[i*N+j]/be[j*N+j];
			for (k = 0; k < j; k++)
			{
				al[i*N+j] = al[i*N+j] - al[i*N+k]*be[k*N+j]/be[j*N+j];
				//printf("al[%2i,%2i]*be[%2i,%2i]/be[%2i,%2i] = %7.2f\n",i,k,k,j,j,j,al[i*N+k]*be[k*N+j]/be[j*N+j]);
			}
			//printf("al[%2i,%2i] = %7.2f\n",i,j,al[i*N+j]);
		}
	}

	for (i = 0; i < N; i++)
	{
		y[i] = b[i]/al[i*N+i];
		for (j = 0; j < i; j++)
		{
			y[i] = y[i] - al[i*N+j]*y[j]/al[i*N+i];
			//printf("al[%2i,%2i]*y[%2i]/al[%2i,%2i] = %7.2f\n",i,j,j,i,i,al[i*N+j]*y[j]/al[i*N+i]);
		}
		//printf("y[%2i] = %7.2f\n",i,y[i]);
	}
	for (i = N-1; i > -1; i--)
	{
		x[i] = y[i]/be[i*N+i];
		for (j = i+1; j < N; j++)
		{
			x[i] = x[i] - be[i*N+j]*x[j]/be[i*N+i];
			//printf("be[%2i,%2i]*x[%2i]/be[%2i,%2i] = %7.2f\n",i,j,j,i,i,be[i*N+j]*x[j]/be[i*N+i]);
		}
		//printf("x[%2i] = %7.2f\n",i,x[i]);
	}
	temp = be[0];
	for (i = 1; i < N; ++i)
	{
		temp = temp*be[i*N+i];
		//printf("beta[%3i,%3i] = %11.4e. temp = %11.4e.\n",i,i,be[i*N+i],temp);
	}
	*Det = temp;

	//Cleanup
	free(al);
	free(be);
	free(y);

	return 0;
}
