#include "BVP_header.h"

int arnoldi(void *params, double A[], double H[], double V[], int j);
int apply_givens_rot(void *params, double H[], double cs[], double sn[], int j);
int givens_rot(double r, double h, double *cs_j1, double *sn_j1);

int BVP_GMRES(void *params, double A[], double b[], double x[], double err[])
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int N = my_ptr->Nparam;
	const int imax = my_ptr->imaxparam;
	const double tol = my_ptr->LStolparam;
	int i,j,k,l,status;
	double temp;
	double b_norm;
	double r_norm;

	double *r0;
	double *betag;
	double *sn;
	double *cs;
	double *V;
	double *H;
	double *y;
	r0 = calloc(N, sizeof(double));
	if(!r0){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	betag = calloc(imax+1, sizeof(double));
	if(!betag){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	sn = calloc(imax+1, sizeof(double));
	if(!sn){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	cs = calloc(imax+1, sizeof(double));
	if(!cs){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	V = calloc(N*(imax+1), sizeof(double));
	if(!V){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	H = calloc((imax+1)*imax, sizeof(double));
	if(!H){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }


	b_norm = 0.0;
	r_norm = 0.0;
	for (i = 0; i < N; ++i)
	{
		temp = 0.0;
		for (j = 0; j < N; ++j)
		{
			temp -= A[i*N+j]*x[j];
		}
		b_norm += b[i]*b[i];
		r0[i] = b[i] + temp;
		r_norm += r0[i]*r0[i];
	}
	b_norm = sqrt(b_norm);
	r_norm = sqrt(r_norm);
	err[0] = r_norm/b_norm;

	betag[0] = r_norm;
	for (i = 0; i < N; ++i)
	{
		//V_0
		V[i*imax+0] = r0[i]/r_norm;
		//printf("V_%i[%i] = %22.15e.\n",0,i,V[i*imax+0]);
	}
	for (k = 0; k < imax; ++k)
	{
		status = arnoldi(params, A, H, V, k);
		//printf("Generated Arnoldi basis.\n");
		status = apply_givens_rot(params, H, cs, sn, k);
		//printf("Applied Givens rotation.\n");

		betag[k+1] = -sn[k]*betag[k];
		betag[k] = cs[k]*betag[k];
		err[k+1] = fabs(betag[k+1])/b_norm;
		//printf("Error [%4i] = %22.15e.\n",k+1,err[k+1]);

		if (err[k+1] < tol)
		{
			printf("Min error [%4i] = %25.18e.\n",k+1,err[k+1]);
			break;
		}
	}
	if (k == imax)
	{
		printf("GMRES FAILURE. Error after %i iterations = %23.16e.\n",k,err[k+1]);
		return -1;
	}

	//Solve system
	y = calloc(k+1, sizeof(double));

	//Direct solving

	for (i = k; i >= 0; --i)
	{
		temp = 0.0;
		for (j = i+1; j <= k; ++j)
		{
			temp += H[i*imax+j]*y[j];
		}
		y[i] = (betag[i]-temp)/H[i*imax+i];
		//printf("y[%i] = %22.15e.\n",i,y[i]);
	}

	for (i = 0; i < N; ++i)
	{
		temp = 0.0;
		for (l = 0; l <= k; ++l)
		{
			temp += V[i*imax+l]*y[l];
		}
		x[i] += temp;
		//printf("x[%i] = %22.15e.\n",i,x[i]);
	}

	//Cleanup
	free(r0);
	free(betag);
	free(sn);
	free(cs);
	free(V);
	free(H);
	free(y);

	return 0;
}



//------------------------------------------
//           Function Definitions
//------------------------------------------

//Arnoldi
int arnoldi(void *params, double A[], double H[], double V[], int j)
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int n = my_ptr->nparam;
	const int m = my_ptr->mparam;
	const int N = my_ptr->Nparam;
	const int imax = my_ptr->imaxparam;
	const double tol = my_ptr->LStolparam;
	int i,k,l;
	double v_norm;
	for (i = 0; i < N; ++i)
	{
		for (l = 0; l < N; ++l)
		{
			V[i*imax+j+1] += A[i*N+l]*V[l*imax+j];
		}
		//printf("q[%i] = %22.15e.\n",i,V[i*imax+j+1]);
	}
	for (i = 0; i <= j; ++i)
	{
		for (l = 0; l < N; ++l)
		{
			H[i*imax+j] += V[l*imax+j+1]*V[l*imax+i];
		}
		//H is an imax+1 by imax matrix
		//V is an N by imax+1 matrix
		//------------------------------------------------------------------------------------------------------------------------
		//printf("h_%i%i = %22.15e.\n",i,j,H[i*imax+j]);
		for (l = 0; l < N; ++l)
		{
			V[l*imax+j+1] -= H[i*imax+j]*V[l*imax+i];
			//printf("q[%i] = q[%i] - %22.15e = %22.15e.\n",i,i,H[i*imax+j]*V[l*imax+i],V[l*imax+j+1]);
		}
	}
	v_norm = 0.0;
	for (l = 0; l < N; ++l)
	{
		v_norm += V[l*imax+j+1]*V[l*imax+j+1];
	}
	v_norm = sqrt(v_norm);
	H[(j+1)*imax+j] = v_norm;
	//printf("h_%i%i = %22.15e.\n",j+1,j,H[(j+1)*imax+j]);
	for (l = 0; l < N; ++l)
	{
		//printf("V_%i[%i] = %22.15e = %22.15e/%22.15e.\n",j+1,l,V[l*imax+j+1]/H[(j+1)*imax+j],V[l*imax+j+1],H[(j+1)*imax+j]);
		V[l*imax+j+1] /= H[(j+1)*imax+j];
	}

	return 0;
}

//Givens Rotation
int apply_givens_rot(void *params, double H[], double cs[], double sn[], int j)
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int n = my_ptr->nparam;
	const int m = my_ptr->mparam;
	const int N = my_ptr->Nparam;
	const int imax = my_ptr->imaxparam;
	const double tol = my_ptr->LStolparam;
	int i,k,l,status;
	double cs_j;
	double sn_j;
	double temp;

	for (i = 0; i < j; ++i)
	{
		temp = cs[i]*H[i*imax+j]+sn[i]*H[(i+1)*imax+j];
		H[(i+1)*imax+j] = -sn[i]*H[i*imax+j]+cs[i]*H[(i+1)*imax+j];
		H[i*imax+j] = temp;
	}

	status = givens_rot(H[j*imax+j], H[(j+1)*imax+j], &cs_j, &sn_j);
	//printf("cs = %22.15e. sn = %22.15e.\n",cs_j,sn_j);
	cs[j]  = cs_j;
	sn[j] = sn_j;

	H[j*imax+j] = cs_j*H[j*imax+j] + sn_j*H[(j+1)*imax+j];
	H[(j+1)*imax+j] = 0.0;
	//for (i = 0; i < j; ++i)
	//{
		//printf("h_%i%i = %22.15e.\n",i,j,H[i*imax+j]);
	//}
	//printf("h_%i%i = %22.15e.\n",j,j,H[j*imax+j]);
	//printf("h_%i%i = %22.15e.\n",j+1,j,H[(j+1)*imax+j]);

	return 0;
}

int givens_rot(double r, double h, double *cs_j, double *sn_j)
{
	double temp;

	if (r == 0.0)
	{
		*cs_j = 0.0;
		*sn_j = 1.0;
	}
	else
	{
		temp = sqrt(r*r+h*h);
		*cs_j = fabs(r)/temp;
		*sn_j = *cs_j/r*h;
	}

	return 0;
}
