/*
Andrew Sullivan
Montana State University
Newton-Rapson method for solving PDE BVP in 1D as coupled equations
Procedure:
Notes:
double format
*/

#include "BVP_header.h"

int BVP_control(int r, int nstep, int p, double tol, double alpha, double beta, double kappa, double r_H)
{
	int i, j, k, l, status;
	double *tmp_ptr;

	//Size
	int n,m,N;
	n = nstep;
	N = n*p;
	double bnew, Dnew;
	int it;
	double temp;
	double btemp, Dtemp, Jtemp, utemp, dxtemp;
	double errorratio, duDxtemp;
	double dxold = (double) (1.0/(n-1));
	int nnew;
	double *xnew;
	double *psinew;


	double gamma = 1.0; //coupling parameter gamma
	const int bvpimax = 50; //maximum iterations
	double LStol = 1.0E-12; //iterative linear solver tolerance
	int LSimax = 2*N; //maximum iterations of iterative solver
	double err_norm;
	double eta_norm;
	double dxnew;
	double w = 1.0; //Start value of relaxed Newton-Raphson parameter
	double wmin = 1.0E-3;
	double wmax = 1.0;
	double wshrink = 0.5;
	double wgrow = 1.5;

	//----- Calloc -----
	double *x;
	double *dxo;
	double *dxn;
	double *R2;
	double *Rab2;
	double *Rabcd2;
	double *RGB;
	x = calloc(n, sizeof(double));
	if(!x){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	dxo = calloc(n, sizeof(double));//dx old
	if(!dxo){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	dxn = calloc(n, sizeof(double));//dx new
	if(!dxn){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	R2 = calloc(n, sizeof(double));//dx new
	if(!R2){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	Rab2 = calloc(n, sizeof(double));//dx new
	if(!Rab2){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	Rabcd2 = calloc(n, sizeof(double));//dx new
	if(!Rabcd2){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	RGB = calloc(n, sizeof(double));//dx new
	if(!RGB){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }


	double *Psi;
	double *dx;
	double *Ddx;
	double *b, *D;
	double *Jac, *err;
	double *b_norm, *D_norm;
	double *Psidx;
	Psi = calloc(N, sizeof(double));
	if(!Psi){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	dx = calloc(N, sizeof(double));
	if(!dx){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	Ddx = calloc(N, sizeof(double));
	if(!Ddx){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	b = calloc(N, sizeof(double));
	if(!b){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	D = calloc(N, sizeof(double));
	if(!D){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	Jac = calloc(N*N, sizeof(double));
	if(!Jac){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	err = calloc(LSimax, sizeof(double));
	if(!err){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	b_norm = calloc(bvpimax, sizeof(double));
	if(!b_norm){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	D_norm = calloc(bvpimax, sizeof(double));
	if(!D_norm){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	Psidx = calloc(N, sizeof(double));
	if(!Psidx){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	//----- Calloc -----


	//Setting struct
	struct param_type my_params;
	my_params.alphaparam = alpha;
	my_params.betaparam = beta;
	my_params.gammaparam = gamma;
	my_params.kappaparam = kappa;
	my_params.nparam = n;
	my_params.pparam = p;
	my_params.Nparam = N;
	my_params.rparam = r;
	my_params.r_Hparam = r_H;
	my_params.tolparam = tol;
	my_params.LStolparam = LStol;
	my_params.imaxparam = LSimax;

	//Initial Conditions
	status = BVP_IC(&my_params,x,Psi);

	for (it = 0; it < bvpimax; ++it)
	{
		printf("\n----- Iteration: %03i -----\n\n",it);

		//Linear system, J dx = -b
		status = BVP_sys(&my_params, r, x, Psi, Jac, b, D, RGB, R2, Rab2, Rabcd2);
		status = BVP_LSsolver(&my_params, x, Psi, Jac, dx, b, err);
		status = BVP_LSsolver(&my_params, x, Psi, Jac, Ddx, D, err);
		//status = BVP_out(&my_params, x, Psi, Jac, dx, b, D, it, RGB, R2, Rab2, Rabcd2);


		btemp = norm(&my_params, b);
		Dtemp = norm(&my_params, D);
		Jtemp = Jac_norm(&my_params, Jac);
		b_norm[it] = btemp;
		D_norm[it] = Dtemp;
		if ((btemp < tol) && (it != 0)) //Convergence condition 2
		{
			printf("\nSUCCESS!! ||b|| = %11.4e < 0.1*||Dx|| = %11.4e or tol = %11.4e.\n\n",btemp,Dtemp,tol);
			status = BVP_physics(&my_params, x,Psi, b, D, RGB, R2, Rab2, Rabcd2);
			//status = BVP_out(&my_params, x, Psi, Jac, dx, b, D, it, RGB, R2, Rab2, Rabcd2); //Save iteration in directory
			status = BVP_out(&my_params, x, Psi, Jac, dx, b, D, -3, RGB, R2, Rab2, Rabcd2); //Save solution in separate directory
			break;
		}

		utemp = norm(&my_params,Psi);
		dxtemp = norm(&my_params,dx);
		duDxtemp = norm(&my_params,Ddx);
		printf("||u||    = %11.4e.\n",utemp);
		printf("||b||    = %11.4e.\n",btemp);
		printf("||Dx||   = %11.4e.\n",Dtemp);
		printf("||J||    = %11.4e.\n",Jtemp);
		printf("||dx||   = %11.4e.\n",dxtemp);
		printf("||dDx||  = %11.4e.\n",duDxtemp);
		printf("Ratio comparison: %11.4e ~ %11.4e.\n",duDxtemp,tol*utemp);


		while (duDxtemp > tol*utemp)
		{
			//Resize
			for (j = 1; j < n-1; ++j)
			{
				dxo[j] = (x[j+1] - x[j-1])/2.0;
				duDxtemp = MAX(MAX(fabs(Ddx[j*p +0]),fabs(Ddx[j*p +1])),fabs(Ddx[j*p +2]));
				errorratio = (1.0/3.0*tol*utemp)/(duDxtemp);
				if (errorratio >= 1.0)
				{
					errorratio = 1.0;
				}
				dxnew = pow(errorratio,1.0/((double) r))*dxo[j];
				dxn[j] = dxnew;
				//printf("%3i: x = %6.3f. dxold = %11.4e. dxnew = %11.4e.\n",j,x[j],dxo[j],dxn[j]);
				//printf("Ddx_f = %11.4e. Ddx_m = %11.4e. Ddx_psi = %11.4e.\n",fabs(Ddx[j*p +0]),fabs(Ddx[j*p +1]),fabs(Ddx[j*p +2]));
			}
			j = 0;
			dxo[j] = (x[j+1] - x[j]);
			duDxtemp = MAX(MAX(fabs(Ddx[j*p +0]),fabs(Ddx[j*p +1])),fabs(Ddx[j*p +2]));
			errorratio = (1.0/3.0*tol*utemp)/(duDxtemp);
			dxnew = pow(errorratio,1.0/((double) r))*dxo[j];
			dxn[j] = dxnew;
			//printf("%3i: x = %6.3f. dxold = %11.4e. dxnew = %11.4e.\n",j,x[j],dxo[j],dxn[j]);
			j = n-1;
			dxo[j] = (x[j] - x[j-1]);
			errorratio = 1.0;
			dxnew = pow(errorratio,1.0/((double) r))*dxo[j];
			dxn[j] = dxnew;
			//printf("%3i: x = %6.3f. dxold = %11.4e. dxnew = %11.4e.\n",j,x[j],dxo[j],dxn[j]);

			nnew = BVP_gridsize(&my_params, r, x, Psi, dxo, dxn, &x, &Psi);
			printf("nnew = %i after nold = %i.\n",nnew,n);

			n = nnew;
			N = n*p;
			LSimax = 2*N;
			my_params.nparam = n;
			my_params.Nparam = N;
			my_params.imaxparam = LSimax;
			printf("N = %i.\n",N);


			//----- Realloc -----
			tmp_ptr = realloc(dxo, sizeof(double) * n);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				dxo = tmp_ptr;
			}
			tmp_ptr = realloc(dxn, sizeof(double) * n);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				dxn = tmp_ptr;
			}
			tmp_ptr = realloc(R2, sizeof(double) * n);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				R2 = tmp_ptr;
			}
			tmp_ptr = realloc(Rab2, sizeof(double) * n);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				Rab2 = tmp_ptr;
			}
			tmp_ptr = realloc(Rabcd2, sizeof(double) * n);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				Rabcd2 = tmp_ptr;
			}
			tmp_ptr = realloc(RGB, sizeof(double) * n);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				RGB = tmp_ptr;
			}
			tmp_ptr = realloc(Ddx, sizeof(double) * N);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				Ddx = tmp_ptr;
			}
			tmp_ptr = realloc(Psidx, sizeof(double) * N);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				Psidx = tmp_ptr;
			}
			tmp_ptr = realloc(b, sizeof(double) * N);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				b = tmp_ptr;
			}
			tmp_ptr = realloc(D, sizeof(double) * N);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				D = tmp_ptr;
			}
			tmp_ptr = realloc(Jac, sizeof(double) * N * N);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				Jac = tmp_ptr;
			}
			tmp_ptr = realloc(dx, sizeof(double) * N);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				dx = tmp_ptr;
			}
			tmp_ptr = realloc(err, sizeof(double) * LSimax);
			if(!tmp_ptr)
			{
				printf("Error! Memory not allocated. Exiting");
				exit(0);
			}
			else
			{
				err = tmp_ptr;
			}
			printf("Successful realloc.\n");

			//----- Realloc -----
			for (i = 0; i < N; ++i)
			{
				dx[i] = 0.0;
				Ddx[i] = 0.0;
			}
			status = BVP_sys(&my_params, r, x, Psi, Jac, b, D, RGB, R2, Rab2, Rabcd2);
			status = BVP_LSsolver(&my_params, x, Psi, Jac, dx, b, err);
			status = BVP_LSsolver(&my_params, x, Psi, Jac, Ddx, D, err);
			btemp = norm(&my_params, b);
			Dtemp = norm(&my_params, D);
			Jtemp = Jac_norm(&my_params, Jac);
			utemp = norm(&my_params, Psi);
			dxtemp = norm(&my_params,dx);
			duDxtemp = norm(&my_params,Ddx);
			printf("After resize ||u||    = %11.4e.\n",utemp);
			printf("After resize ||b||    = %11.4e.\n",btemp);
			printf("After resize ||Dx||   = %11.4e.\n",Dtemp);
			printf("After resize ||dx||   = %11.4e.\n",dxtemp);
			printf("After resize ||dDx||  = %11.4e.\n",duDxtemp);
			printf("After resize Ratio comparison: If %11.4e < %11.4e stepsize is fine.\n",duDxtemp,tol*utemp);
			status = BVP_out(&my_params, x, Psi, Jac, dx, b, D, it, RGB, R2, Rab2, Rabcd2);

		}
		printf("||dDx|| = %11.4e < tol ||u||.\n",duDxtemp);
		printf("Step accuracy is sufficient.\n");


		//Check convergence
		for (j = 0; j < N; ++j)
		{
			Psidx[j] = Psi[j] + w*dx[j];
		}
		status = BVP_sys(&my_params, r, x, Psidx, Jac, b, D, RGB, R2, Rab2, Rabcd2);
		bnew = norm(&my_params, b);
		Dnew = norm(&my_params, D);
		printf("New ||b||  = %11.4e, old ||b||  = %11.4e.\n",bnew,btemp);
		printf("New ||Dx|| = %11.4e, old ||Dx|| = %11.4e.\n",Dnew,Dtemp);
		while (bnew > btemp)
		{
			if (w < wmin)
			{
				printf("w below wmin = %6.3f.\n",wmin);
				break;
			}
			printf("Not converging, reducing w = %6.3f.\n",w);
			w *= wshrink;
			printf("new w = %6.3f.\n",w);
			for (j = 0; j < N; ++j)
			{
				Psidx[j] = Psi[j] + w*dx[j];
			}
			status = BVP_sys(&my_params, r, x, Psidx, Jac, b, D, RGB, R2, Rab2, Rabcd2);
			bnew = norm(&my_params, b);
			Dnew = norm(&my_params, D);
			printf("New ||b||  = %11.4e, old ||b||  = %11.4e.\n",bnew,btemp);
			printf("New ||Dx|| = %11.4e, old ||Dx|| = %11.4e.\n",Dnew,Dtemp);
		}
		if (w < wmin)
		{
			break;
		}

		//If converging
		for (j = 0; j < N; ++j)
		{
			Psi[j] += w*dx[j];
		}

	}

	//Cleanup
	free(x);
	free(dxo);
	free(dxn);
	free(Psi);
	free(dx);
	free(Ddx);
	free(b);
	free(D);
	free(Jac);
	free(err);
	free(b_norm);
	free(D_norm);
	free(Psidx);
	free(RGB);
	free(R2);
	free(Rab2);
	free(Rabcd2);

	return 0;
}





//------------------------------------------
//           Function Definitions
//------------------------------------------

double norm(void *params, double v[])
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int N = my_ptr->Nparam;
	int i;
	double temp = 0.0;
	//Euclidean Norm
	/*
	for (i = 0; i < N; ++i)
	{
		temp += v[i]*v[i];
	}
	temp = sqrt(temp);
	*/
	//Max Norm
	for (i = 0; i < N; ++i)
	{
		if (temp < fabs(v[i]))
		{
			temp = fabs(v[i]);
		}
	}

	return temp;
}

double Jac_norm(void *params, double v[])
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int p = my_ptr->pparam;
	const int N = my_ptr->Nparam;
	int i, j;
	double temp = 0.0;
	//Max Norm
	for (i = p; i < N-p; ++i)
	{
		for (j = p; j < N-p; ++j)
		{
			if (temp < fabs(v[i*N+j]))
			{
				temp = fabs(v[i*N+j]);
			}
		}
	}

	return temp;
}
