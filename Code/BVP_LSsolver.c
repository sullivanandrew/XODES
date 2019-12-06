#include "BVP_header.h"

int BVP_LSsolver(void *params, double x[], double Psi[], double Jac[], double dx[], double b[], double err[])
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int n = my_ptr->nparam;
	const int m = my_ptr->mparam;
	const int N = my_ptr->Nparam;
	const double r_H = my_ptr->r_Hparam;
	const int imax = my_ptr->imaxparam;
	int i,j,k,l,status;
	double Det;
	clock_t begin, end;
	double time_spent;

	begin = clock();
	status = BVP_GMRES(params, Jac, b, dx, err);
	//status = BVP_BICO(params, Jac, b, dx, err);
	//status = BVP_LUdcmp(params, Jac, b, dx, &Det);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Elapsed Time for GMRES Solver = %11.4f.\n",time_spent);
	//printf("Elapsed Time for BICO Solver = %11.4f.\n",time_spent);
	//printf("Elapsed Time for LU Decomp Solver = %11.4f.\n",time_spent);
	if (status == -1)
	{
		for (i = 0; i < N; ++i)
		{
			dx[i] = 0.0;
		}
		for (i = 0; i < imax; ++i)
		{
			//printf("err[%3i] = %11.4e.\n",i,err[i]);
			err[i] = 0.0;
		}
		printf("GMRES failed. Trying BICO.\n");
		//printf("BICO failed. Trying GMRES.\n");
		begin = clock();
		status = BVP_BICO(params, Jac, b, dx, err);
		//status = BVP_GMRES(params, Jac, b, dx, err);
		end = clock();
		time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
		printf("Elapsed Time for BICO Solver = %11.4f.\n",time_spent);
		//printf("Elapsed Time for GMRES Solver = %11.4f.\n",time_spent);
		if (status == -1)
		{
			printf("\n\nFAILURE: Linear solver has failed.\n");
			return -1;
		}
	}

	return 0;
}
