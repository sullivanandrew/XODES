#include "BVP_header.h"

int BVP_gridsize(void *params, int r, double x[], double Psi[], double dxo[], double dxn[], double **xp, double **pp)
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int n = my_ptr->nparam;
	const int p = my_ptr->pparam;
	const int N = my_ptr->Nparam;
	int i,j,k,l;
	int status;
	int pn;
	int nnew, Nnew;
	double xk;
	double *tmp_ptr;

	double dxtemp, vi;
	double x0, x1;
	double temp, temp2;
	int nmax;
	temp = 0.0;
	temp2 = 0.0;

	double *xtemp;
	double *psitemp;
	xtemp = calloc(n, sizeof(double));//x Right
	if(!xtemp){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	psitemp = calloc(N, sizeof(double));
	if(!psitemp){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	nmax = 0;
	for (i = 0; i < n; ++i)
	{
		if (temp2 < fabs(dxn[i]))
		{
			temp2 = fabs(dxn[i]);
			nmax = i;
		}
		//printf("%3i: x = %6.3f, dxo = %6.3f, dxn = %7.4f.\n",i, x[i], dxo[i], dxn[i]);
		if (i == 0 || i == n-1)
		{
			temp += dxn[i]/2.0;
		}
		else
		{
			temp += dxn[i];
		}
	}
	nnew = (int) (((double) n)*3.0/2.0/temp);
	printf("total step distance = %11.4e. Estimated new points = %i.\n",temp,nnew);
	printf("max stepsize is %7.4f at %i.\n",temp2,nmax);

	double *xL;
	double *dxL;
	double *xR;
	double *dxR;
	xL = calloc(nnew, sizeof(double));//x Left
	if(!xL){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	dxL = calloc(nnew, sizeof(double));//x Left
	if(!dxL){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	xR = calloc(nnew, sizeof(double));//x Right
	if(!xR){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	dxR = calloc(nnew, sizeof(double));//x Right
	if(!dxR){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }

	//Temp
	for (i = 0; i < n; ++i)
	{
		xtemp[i] = x[i];
		//printf("xtemp[%3i] = %7.4f.\n",i,xtemp[i]);
	}
	for (i = 0; i < N; ++i)
	{
		psitemp[i] = Psi[i];
		//printf("psitemp[%3i] = %11.4e.\n",i,psitemp[i]);
	}

	//Left to middle
	j = 0;
	xL[j] = xtemp[0];
	dxL[j] = dxn[0];

	for (i = 0; i < nmax; ++i)
	{
		x0 = xL[j] + dxL[j]/2.0;
		dxtemp = dxn[i];
		x1 = xtemp[i] + dxo[i]/2.0;
		while (x0 + dxtemp/2.0 < x1)
		{
			xL[j+1] = x0 + dxtemp/2.0;
			dxL[j+1] = dxtemp;
			x0 = xL[j+1] + dxtemp/2.0;
			++j;
		}

	}

	//Right to middle
	k = 0;
	xR[k] = xtemp[n-1];
	dxR[k] = dxn[n-1];
	for (i = n-1; i > nmax; --i)
	{
		x1 = xR[k] - dxR[k]/2.0;
		dxtemp = dxn[i];
		x0 = xtemp[i] - dxo[i]/2.0;
		while (x1 - dxtemp/2.0 > x0)
		{
			xR[k+1] = x1 - dxtemp/2.0;
			dxR[k+1] = dxtemp;
			x1 = xR[k+1] - dxtemp/2.0;
			++k;
		}

	}

	//Center
	dxtemp = dxn[nmax];
	x0 = xL[j] + dxL[j]/2.0;
	x1 = xR[k] - dxR[k]/2.0;
	vi = (x1 - x0)/dxtemp;
	while (vi > 2.0)
	{
		//Left
		xL[j+1] = x0 + dxtemp/2.0;
		dxL[j+1] = dxtemp;
		x0 = xL[j+1] + dxtemp/2.0;
		++j;
		//Right
		xR[k+1] = x1 - dxtemp/2.0;
		dxR[k+1] = dxtemp;
		x1 = xR[k+1] - dxtemp/2.0;
		++k;
		vi = vi - 2.0;
	}
	if (vi > 1.0)
	{
		dxtemp = (x1-x0)/2.0;
		//Left
		xL[j+1] = x0 + dxtemp/2.0;
		dxL[j+1] = dxtemp;
		x0 = xL[j+1] + dxtemp/2.0;
		++j;
		//Right
		xR[k+1] = x1 - dxtemp/2.0;
		dxR[k+1] = dxtemp;
		x1 = xR[k+1] - dxtemp/2.0;
		++k;
	}
	else
	{
		dxtemp = (x1 - x0);
		//Left
		xL[j+1] = x0 + dxtemp/2.0;
		dxL[j+1] = dxtemp;
		x0 = xL[j+1] + dxtemp/2.0;
		++j;
	}

	nnew = j+k+2;
	Nnew = nnew*p;
	printf("New points = %i.\n",nnew);
	tmp_ptr = realloc(x, sizeof(double) * nnew);//x Left
	if(!tmp_ptr)
	{
		printf("Error! Memory not allocated. Exiting");
		exit(0);
	}
	else
	{
		x = tmp_ptr;
	}
	tmp_ptr = realloc(Psi, sizeof(double) * Nnew);//x Left
	if(!tmp_ptr)
	{
		printf("Error! Memory not allocated. Exiting");
		exit(0);
	}
	else
	{
		Psi = tmp_ptr;
	}

	*xp = x;
	*pp = Psi;
	printf("xn pointer after realloc  = %p.\n",*xp);
	printf("pn pointer after realloc  = %p.\n",*pp);

	l = 0;
	for (i = 0; i <= j; ++i)
	{
		//printf("xL = %7.4f.\n",xL[i]);
		x[l] = xL[i];
		++l;
	}
	for (i = k; i >= 0; --i)
	{
		//printf("xR = %7.4f.\n",xR[i]);
		x[l] = xR[i];
		++l;
	}

	//Interpolate each Psi
	status = BVP_interp(params, r, x, Psi, nnew, xtemp, psitemp, n);

	free(xtemp);
	free(psitemp);
	free(xL);
	free(dxL);
	free(xR);
	free(dxR);

	return nnew;
}

int BVP_interp(void *params, int r, double x[], double Psi[], int nnew, double xtemp[], double Psitemp[], int ntemp)
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int p = my_ptr->pparam;
	int Nnew = nnew*p;
	int Ntemp = ntemp*p;
	int i,j,k,l,a,status;
	int pn;
	double xk;
	int I, oneside, twoside;
	int Itemp;

	double *a1_a;
	double *a1x_a;
	double *a1xx_a;
	a1_a = calloc(r+4, sizeof(double));
	a1x_a = calloc(r+4, sizeof(double));
	a1xx_a = calloc(r+4, sizeof(double));
	if(!a1_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a1x_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a1xx_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	double *a2_a;
	double *a2x_a;
	double *a2xx_a;
	a2_a = calloc(r+4, sizeof(double));
	a2x_a = calloc(r+4, sizeof(double));
	a2xx_a = calloc(r+4, sizeof(double));
	if(!a2_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a2x_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a2xx_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	double *a3_a;
	double *a3x_a;
	double *a3xx_a;
	a3_a = calloc(r+4, sizeof(double));
	a3x_a = calloc(r+4, sizeof(double));
	a3xx_a = calloc(r+4, sizeof(double));
	if(!a3_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a3x_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a3xx_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	int *s;
	s = (int*)calloc(r+3, sizeof(int));
	if(!s){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	double ps1k,dps1k,d2ps1k;
	double ps1dk,dps1dk,d2ps1dk;
	double ps2k,dps2k,d2ps2k;
	double ps2dk,dps2dk,d2ps2dk;
	double ps3k,dps3k,d2ps3k;
	double ps3dk,dps3dk,d2ps3dk;

	for (j = 0; j < nnew; ++j)
	{
		k = j*p;
		xk = x[j];
		//printf("xk = %11.4e.\n",xk);
		I = 0;
		for (i = 0; i < ntemp; ++i)
		{
			if (xtemp[i] >= xk)
			{
				I = i;
				//printf("I = %i.\n",I);
				break;
			}
		}
		//printf("xtemp[Itemp] = %11.4e.\n",xtemp[I]);

		ps1k = 0.0, dps1k = 0.0, d2ps1k = 0.0;
		ps1dk = 0.0, dps1dk = 0.0, d2ps1dk = 0.0;
		ps2k = 0.0, dps2k = 0.0, d2ps2k = 0.0;
		ps2dk = 0.0, dps2dk = 0.0, d2ps2dk = 0.0;
		ps3k = 0.0, dps3k = 0.0, d2ps3k = 0.0;
		ps3dk = 0.0, dps3dk = 0.0, d2ps3dk = 0.0;
		oneside = 0;
		twoside = 0;

		//Compute stencil
		status = steninit(s, &oneside, &twoside, I, r, ntemp);
		//Compute Newton polynomial representation for each function
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, xtemp, Psitemp, 0, a1_a, a1x_a, a1xx_a, &ps1dk, &dps1dk, &d2ps1dk, xk);
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, xtemp, Psitemp, 1, a2_a, a2x_a, a2xx_a, &ps2dk, &dps2dk, &d2ps2dk, xk);
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, xtemp, Psitemp, 2, a3_a, a3x_a, a3xx_a, &ps3dk, &dps3dk, &d2ps3dk, xk);

		for (a = 0; a <= r+oneside; ++a)
		{
			ps1k   += Psitemp[(I+s[a])*p +0]*a1_a[a];
			dps1k  += Psitemp[(I+s[a])*p +0]*a1x_a[a];
			d2ps1k += Psitemp[(I+s[a])*p +0]*a1xx_a[a];
			ps2k   += Psitemp[(I+s[a])*p +1]*a2_a[a];
			dps2k  += Psitemp[(I+s[a])*p +1]*a2x_a[a];
			d2ps2k += Psitemp[(I+s[a])*p +1]*a2xx_a[a];
			ps3k   += Psitemp[(I+s[a])*p +2]*a3_a[a];
			dps3k  += Psitemp[(I+s[a])*p +2]*a3x_a[a];
			d2ps3k += Psitemp[(I+s[a])*p +2]*a3xx_a[a];
		}

		pn = 0;
		Psi[k +pn] = ps1k;
		pn = 1;
		Psi[k +pn] = ps2k;
		pn = 2;
		Psi[k +pn] = ps3k;


	}

	//Cleanup
	free(a1_a);
	free(a1x_a);
	free(a1xx_a);
	free(a2_a);
	free(a2x_a);
	free(a2xx_a);
	free(a3_a);
	free(a3x_a);
	free(a3xx_a);
	free(s);

	return 0;
}
