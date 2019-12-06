#include "BVP_header.h"

int BVP_sys(void *params, const int r, double x[], double Psi[], double Jac[], double bvec[], double Dvec[], double RGB[], double R2[], double Rab2[], double Rabcd2[])
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int n = my_ptr->nparam;
	const int p = my_ptr->pparam;
	const int N = my_ptr->Nparam;
	const double r_H = my_ptr->r_Hparam;
	const double alpha = my_ptr->alphaparam;
	const double beta = my_ptr->betaparam;
	const double gamma = my_ptr->gammaparam;
	const double kappa = my_ptr->kappaparam;
	int i,j,k,l,a,status;
	int pn;

	const double sgn = -1.0;
	int I, oneside, twoside;
	double xhoriz, xinf;
	double psihoriz, psiinf;

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
	double xk,thk;

	double Tcoupl, Kcoupl;
	double Gpsicoupl, GBcoupl;

	double Gtt, Ttt, Ktt;
	double dGttdf2,dGttdf1,dGttdf0;
	double dGttdm2,dGttdm1,dGttdm0;
	double dGttdp2,dGttdp1,dGttdp0;
	double dTttdf2,dTttdf1,dTttdf0;
	double dTttdm2,dTttdm1,dTttdm0;
	double dTttdp2,dTttdp1,dTttdp0;
	double dKttdf2,dKttdf1,dKttdf0;
	double dKttdm2,dKttdm1,dKttdm0;
	double dKttdp2,dKttdp1,dKttdp0;

	double Grr, Trr, Krr;
	double dGrrdf2,dGrrdf1,dGrrdf0;
	double dGrrdm2,dGrrdm1,dGrrdm0;
	double dGrrdp2,dGrrdp1,dGrrdp0;
	double dTrrdf2,dTrrdf1,dTrrdf0;
	double dTrrdm2,dTrrdm1,dTrrdm0;
	double dTrrdp2,dTrrdp1,dTrrdp0;
	double dKrrdf2,dKrrdf1,dKrrdf0;
	double dKrrdm2,dKrrdm1,dKrrdm0;
	double dKrrdp2,dKrrdp1,dKrrdp0;

	double Gpsi, GB, GBcalc, R2calc, Rab2calc, Rabcd2calc;
	double dGpsidf2, dGpsidf1, dGpsidf0;
	double dGpsidm2, dGpsidm1, dGpsidm0;
	double dGpsidp2, dGpsidp1, dGpsidp0;
	double dGBdf2, dGBdf1, dGBdf0;
	double dGBdm2, dGBdm1, dGBdm0;
	double dGBdp2, dGBdp1, dGBdp0;


	Tcoupl = 1.0; //Stress energy tensor coupling
	Kcoupl = 1.0*alpha; //Interaction tensor coupling
	Gpsicoupl = 1.0;
	GBcoupl = alpha;

	//Initialize system to zero
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < N; ++j)
		{
			Jac[i*N+j] = 0.0;
		}
		bvec[i] = 0.0;
		Dvec[i] = 0.0;
	}

	//Boundary Conditions
	const double BC0 = 0.0;
	const double BC1 = 1.0;
	xhoriz = x[0];
	xinf = x[n-1];

	//Function 1
	pn = 0;
	psihoriz = BC0;
	psiinf = BC1;

	k = 0*p; //phi(x=0)
	bvec[k +pn] = sgn*(Psi[k +pn] - psihoriz);
	Jac[(k +pn)*N+ (k +pn)] = 1.0;
	k = (n-1)*p; //phi(x=1)
	bvec[k +pn] = sgn*(Psi[k +pn] - psiinf);
	Jac[(k +pn)*N +(k +pn)] = 1.0;

	//Function 2
	pn = 1;
	psihoriz = 16.0/pow(1.0+xhoriz,4.0);
	psiinf = BC1;

	k = 0*p; //phi(x=0)
	bvec[k +pn] = sgn*(Psi[k +pn] - psihoriz);
	Jac[(k +pn)*N+ (k +pn)] = 1.0;
	k = (n-1)*p; //phi(x=1)
	bvec[k +pn] = sgn*(Psi[k +pn] - psiinf);
	Jac[(k +pn)*N +(k +pn)] = 1.0;

	//Function 3
	pn = 2;
	psiinf = BC0; //G_psi bar and nonbar

	//First derivative BC
	k = 0*p; //phi(x=0)
	ps3k = 0.0, dps3k = 0.0, d2ps3k = 0.0;
	ps3dk = 0.0, dps3dk = 0.0, d2ps3dk = 0.0;
	I = 0;
	xk = x[I];
	oneside = 0;
	twoside = 0;
	//Compute stencil
	status = steninit(s, &oneside, &twoside, I, r, n);
	status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 2, a3_a, a3x_a, a3xx_a, &ps3dk, &dps3dk, &d2ps3dk, xk);

	for (a = 0; a <= r+oneside; ++a)
	{
		ps3k += Psi[(I+s[a])*p +pn]*a3_a[a];
		dps3k += Psi[(I+s[a])*p +pn]*a3x_a[a];
		d2ps3k += Psi[(I+s[a])*p +pn]*a3xx_a[a];
	}


	bvec[k +pn] = sgn*(dps3k - BC0);
	for (a = 0; a <= r+oneside; ++a)
	{
		Jac[(k +pn)*N+ (I+s[a])*p +pn] = a3x_a[a];
	}
	Dvec[k +pn] = dps3dk;

	//printf("IC: b = %11.4e. dps3k = %11.4e. D = %11.4e.\n",bvec[k +pn],dps3k,Dvec[k +pn]);

	k = (n-1)*p; //phi(x=1)
	bvec[k +pn] = sgn*(Psi[k +pn] - psiinf);
	Jac[(k +pn)*N +(k +pn)] = 1.0;


	//Inside grid
	for (j = 1; j < n-1; ++j)
	{
		k = j*p;
		I = j;
		xk = x[j];

		ps1k = 0.0, dps1k = 0.0, d2ps1k = 0.0;
		ps1dk = 0.0, dps1dk = 0.0, d2ps1dk = 0.0;
		ps2k = 0.0, dps2k = 0.0, d2ps2k = 0.0;
		ps2dk = 0.0, dps2dk = 0.0, d2ps2dk = 0.0;
		ps3k = 0.0, dps3k = 0.0, d2ps3k = 0.0;
		ps3dk = 0.0, dps3dk = 0.0, d2ps3dk = 0.0;
		oneside = 0;
		twoside = 0;

		//Compute stencil
		status = steninit(s, &oneside, &twoside, j, r, n);
		//Compute Newton polynomial representation for each function
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 0, a1_a, a1x_a, a1xx_a, &ps1dk, &dps1dk, &d2ps1dk, xk);
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 1, a2_a, a2x_a, a2xx_a, &ps2dk, &dps2dk, &d2ps2dk, xk);
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 2, a3_a, a3x_a, a3xx_a, &ps3dk, &dps3dk, &d2ps3dk, xk);

		for (a = 0; a <= r+oneside; ++a)
		{
			ps1k   += Psi[(I+s[a])*p +0]*a1_a[a];
			dps1k  += Psi[(I+s[a])*p +0]*a1x_a[a];
			d2ps1k += Psi[(I+s[a])*p +0]*a1xx_a[a];
			ps2k   += Psi[(I+s[a])*p +1]*a2_a[a];
			dps2k  += Psi[(I+s[a])*p +1]*a2x_a[a];
			d2ps2k += Psi[(I+s[a])*p +1]*a2xx_a[a];
			ps3k   += Psi[(I+s[a])*p +2]*a3_a[a];
			dps3k  += Psi[(I+s[a])*p +2]*a3x_a[a];
			d2ps3k += Psi[(I+s[a])*p +2]*a3xx_a[a];
		}


		//-------------------- GR FULL bar g_tt --------------------
		pn = 0;
		Gtt = Gtout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		Ttt = Ttout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		Ktt = Ktout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGttdf2 = dGtdf2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGttdf1 = dGtdf1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGttdf0 = dGtdf0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGttdm2 = dGtdm2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGttdm1 = dGtdm1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGttdm0 = dGtdm0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGttdp2 = dGtdp2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGttdp1 = dGtdp1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGttdp0 = dGtdp0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		dTttdf2 = dTtdf2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTttdf1 = dTtdf1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTttdf0 = dTtdf0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTttdm2 = dTtdm2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTttdm1 = dTtdm1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTttdm0 = dTtdm0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTttdp2 = dTtdp2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTttdp1 = dTtdp1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTttdp0 = dTtdp0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		dKttdf2 = dKtdf2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKttdf1 = dKtdf1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKttdf0 = dKtdf0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKttdm2 = dKtdm2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKttdm1 = dKtdm1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKttdm0 = dKtdm0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKttdp2 = dKtdp2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKttdp1 = dKtdp1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKttdp0 = dKtdp0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		//---------- b ----------
		bvec[k +pn] = sgn*(Gtt - Tcoupl*Ttt + Kcoupl*Ktt);
		//---------- Jac ----------
		for (a = 0; a <= r+oneside; ++a)
		{
			Jac[(k +pn)*N+ (I+s[a])*p +0] += dGttdf2*a1xx_a[a] +dGttdf1*a1x_a[a] +dGttdf0*a1_a[a] -Tcoupl*(dTttdf2*a1xx_a[a] +dTttdf1*a1x_a[a] +dTttdf0*a1_a[a]) +Kcoupl*(dKttdf2*a1xx_a[a] +dKttdf1*a1x_a[a] +dKttdf0*a1_a[a]);

			Jac[(k +pn)*N+ (I+s[a])*p +1] += dGttdm2*a2xx_a[a] +dGttdm1*a2x_a[a] +dGttdm0*a2_a[a] -Tcoupl*(dTttdm2*a2xx_a[a] +dTttdm1*a2x_a[a] +dTttdm0*a2_a[a]) +Kcoupl*(dKttdm2*a2xx_a[a] +dKttdm1*a2x_a[a] +dKttdm0*a2_a[a]);

			Jac[(k +pn)*N+ (I+s[a])*p +2] += dGttdp2*a3xx_a[a] +dGttdp1*a3x_a[a] +dGttdp0*a3_a[a] -Tcoupl*(dTttdp2*a3xx_a[a] +dTttdp1*a3x_a[a] +dTttdp0*a3_a[a]) +Kcoupl*(dKttdp2*a3xx_a[a] +dKttdp1*a3x_a[a] +dKttdp0*a3_a[a]);

		}

		//---------- D ----------
		Dvec[k +pn] = (dGttdf2 -Tcoupl*dTttdf2 +Kcoupl*dKttdf2)*d2ps1dk + (dGttdf1 -Tcoupl*dTttdf1 +Kcoupl*dKttdf1)*dps1dk + (dGttdm2 -Tcoupl*dTttdm2 +Kcoupl*dKttdm2)*d2ps2dk + (dGttdm1 -Tcoupl*dTttdm1 +Kcoupl*dKttdm1)*dps2dk + (dGttdp2 -Tcoupl*dTttdp2 +Kcoupl*dKttdp2)*d2ps3dk + (dGttdp1 -Tcoupl*dTttdp1 +Kcoupl*dKttdp1)*dps3dk; //Nonlinear



		//-------------------- GR FULL bar g_thth --------------------
		pn = 1;
		Grr = Grout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		Trr = Trout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		Krr = Krout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		dGrrdf2 = dGrdf2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGrrdf1 = dGrdf1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGrrdf0 = dGrdf0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGrrdm2 = dGrdm2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGrrdm1 = dGrdm1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGrrdm0 = dGrdm0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGrrdp2 = dGrdp2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGrrdp1 = dGrdp1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGrrdp0 = dGrdp0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		dTrrdf2 = dTrdf2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTrrdf1 = dTrdf1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTrrdf0 = dTrdf0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTrrdm2 = dTrdm2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTrrdm1 = dTrdm1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTrrdm0 = dTrdm0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTrrdp2 = dTrdp2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTrrdp1 = dTrdp1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dTrrdp0 = dTrdp0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		dKrrdf2 = dKrdf2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKrrdf1 = dKrdf1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKrrdf0 = dKrdf0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKrrdm2 = dKrdm2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKrrdm1 = dKrdm1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKrrdm0 = dKrdm0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKrrdp2 = dKrdp2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKrrdp1 = dKrdp1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dKrrdp0 = dKrdp0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		//---------- b ----------
		bvec[k +pn] = sgn*(Grr -Tcoupl*Trr +Kcoupl*Krr);
		//---------- Jac ----------
		for (a = 0; a <= r+oneside; ++a)
		{
			Jac[(k +pn)*N+(I+s[a])*p +0] += dGrrdf2*a1xx_a[a] +dGrrdf1*a1x_a[a] +dGrrdf0*a1_a[a] -Tcoupl*(dTrrdf2*a1xx_a[a] +dTrrdf1*a1x_a[a] +dTrrdf0*a1_a[a]) +Kcoupl*(dKrrdf2*a1xx_a[a] +dKrrdf1*a1x_a[a] +dKrrdf0*a1_a[a]);

			Jac[(k +pn)*N+(I+s[a])*p +1] += dGrrdm2*a2xx_a[a] +dGrrdm1*a2x_a[a] +dGrrdm0*a2_a[a] -Tcoupl*(dTrrdm2*a2xx_a[a] +dTrrdm1*a2x_a[a] +dTrrdm0*a2_a[a]) +Kcoupl*(dKrrdm2*a2xx_a[a] +dKrrdm1*a2x_a[a] +dKrrdm0*a2_a[a]);

			Jac[(k +pn)*N+(I+s[a])*p +2] += dGrrdp2*a3xx_a[a] +dGrrdp1*a3x_a[a] +dGrrdp0*a3_a[a] -Tcoupl*(dTrrdp2*a3xx_a[a] +dTrrdp1*a3x_a[a] +dTrrdp0*a3_a[a]) +Kcoupl*(dKrrdp2*a3xx_a[a] +dKrrdp1*a3x_a[a] +dKrrdp0*a3_a[a]);

		}

		//---------- D ----------
		Dvec[k +pn] = (dGrrdf2 -Tcoupl*dTrrdf2 +Kcoupl*dKrrdf2)*d2ps1dk + (dGrrdf1 -Tcoupl*dTrrdf1 +Kcoupl*dKrrdf1)*dps1dk + (dGrrdm2 -Tcoupl*dTrrdm2 +Kcoupl*dKrrdm2)*d2ps2dk + (dGrrdm1 -Tcoupl*dTrrdm1 +Kcoupl*dKrrdm1)*dps2dk + (dGrrdp2 -Tcoupl*dTrrdp2 +Kcoupl*dKrrdp2)*d2ps3dk + (dGrrdp1 -Tcoupl*dTrrdp1 +Kcoupl*dKrrdp1)*dps3dk; //Nonlinear



		//-------------------- Scalar Field on Fixed SCHW Background --------------------
		pn = 2;
		Gpsi = Gpsiout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		GB = GBout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		GBcalc = GBcalcout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		R2calc = R2calcout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		Rab2calc = Rab2calcout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		Rabcd2calc = Rabcd2calcout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		dGpsidf2 = dGpsidf2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGpsidf1 = dGpsidf1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGpsidf0 = dGpsidf0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGpsidm2 = dGpsidm2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGpsidm1 = dGpsidm1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGpsidm0 = dGpsidm0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGpsidp2 = dGpsidp2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGpsidp1 = dGpsidp1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGpsidp0 = dGpsidp0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		dGBdf2 = dGBdf2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGBdf1 = dGBdf1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGBdf0 = dGBdf0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGBdm2 = dGBdm2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGBdm1 = dGBdm1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGBdm0 = dGBdm0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGBdp2 = dGBdp2out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGBdp1 = dGBdp1out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dGBdp0 = dGBdp0out(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		//---------- b ----------
		bvec[k +pn] = sgn*(Gpsicoupl*Gpsi + GBcoupl*GB);
		//---------- Jac ----------
		for (a = 0; a <= r+oneside; ++a)
		{
			Jac[(k +pn)*N+(I+s[a])*p +0] += Gpsicoupl*(dGpsidf2*a1xx_a[a] + dGpsidf1*a1x_a[a] + dGpsidf0*a1_a[a]) + GBcoupl*(dGBdf2*a1xx_a[a] + dGBdf1*a1x_a[a] + dGBdf0*a1_a[a]);

			Jac[(k +pn)*N+(I+s[a])*p +1] += Gpsicoupl*(dGpsidm2*a2xx_a[a] + dGpsidm1*a2x_a[a] + dGpsidm0*a2_a[a]) + GBcoupl*(dGBdm2*a2xx_a[a] + dGBdm1*a2x_a[a] + dGBdm0*a2_a[a]);

			Jac[(k +pn)*N+(I+s[a])*p +2] += Gpsicoupl*(dGpsidp2*a3xx_a[a] + dGpsidp1*a3x_a[a] + dGpsidp0*a3_a[a]) + GBcoupl*(dGBdp2*a3xx_a[a] + dGBdp1*a3x_a[a] + dGBdp0*a3_a[a]);
		}
		RGB[j] = GBcalc;
		R2[j] = R2calc;
		Rab2[j] = Rab2calc;
		Rabcd2[j] = Rabcd2calc;

		//---------- D ----------
		Dvec[k +pn] = (Gpsicoupl*dGpsidf2 +GBcoupl*dGBdf2)*d2ps1dk + (Gpsicoupl*dGpsidf1 +GBcoupl*dGBdf1)*dps1dk + (Gpsicoupl*dGpsidm2 +GBcoupl*dGBdm2)*d2ps2dk + (Gpsicoupl*dGpsidm1 +GBcoupl*dGBdm1)*dps2dk + (Gpsicoupl*dGpsidp2 +GBcoupl*dGBdp2)*d2ps3dk + (Gpsicoupl*dGpsidp1 +GBcoupl*dGBdp1)*dps3dk;


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





//------------------------------------------
//           Function Definitions
//------------------------------------------
int BVP_NPcalc(void *params, int r, int rn, int I, int s[], double x[], double y[], int pn, double a_a[], double ax_a[], double axx_a[], double *ydk, double *dydk, double *d2ydk, double xk)
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int p = my_ptr->pparam;
	int k,t,q,a,b,c;
	double s_ba;
	double P_b;
	double Px_b;
	double Pxx_b;

	double temp;

	for (a = 0; a <= rn; ++a)
	{
		a_a[a] = 0.0;
		ax_a[a] = 0.0;
		axx_a[a] = 0.0;
		for (b = a; b <= rn; ++b)
		{
			P_b = 0.0;
			s_ba = 1.0;
			for (c = 0; c <= b; ++c)
			{
				if (c != a)
				{
					s_ba /= (x[I+s[a]]-x[I+s[c]]);
				}
			}
			temp = 1.0;
			for (c = 0; c <= b-1; ++c)
			{
				temp *= (xk - x[I+s[c]]);
			}
			P_b += temp;
			//P'
			Px_b = 0.0;
			for (q = 0; q <= b-1; ++q)
			{
				temp = 1.0;
				for (c = 0; c <= b-1; ++c)
				{
					if (c != q)
					{
						temp *= (xk - x[I+s[c]]);
					}
				}
				Px_b += temp;
			}
			//P''
			Pxx_b = 0.0;
			for (t = 0; t <= b-1; ++t)
			{
				for (q = 0; q <= b-1; ++q)
				{
					if (q != t)
					{
						temp = 1.0;
						for (c = 0; c <= b-1; ++c)
						{
							if (c != q && c != t)
							{
								temp *= (xk - x[I+s[c]]);
							}
						}
						Pxx_b += temp;
					}
				}
			}
			if (b > r)
			{
				*ydk += y[(I +s[a])*p +pn]*s_ba*P_b;
				*dydk += y[(I +s[a])*p +pn]*s_ba*Px_b;
				*d2ydk += y[(I +s[a])*p +pn]*s_ba*Pxx_b;
			}
			else
			{
				a_a[a] += s_ba*P_b;
				ax_a[a] += s_ba*Px_b;
				axx_a[a] += s_ba*Pxx_b;
			}
		}
	}

	return 0;
}

int steninit(int s[], int *oneside, int *twoside, int k, int r, int N)
{
	int l;
	s[0] = 0;
	s[(r+2)+1] = 0;
	for (l = 0; l < (r+2)/2; ++l)
	{
		s[2*l+1] = l+1;
		s[2*l+2] = -(l+1);
	}
	if (k < r/2)
	{
		*twoside = 0;
		*oneside = 1;
		for (l = 2*k+2; l <= (r+2)+1; ++l)
		{
			s[l] = s[l-1] + 1;
		}
	}
	else if (k < (r+2)/2)
	{
		*twoside = 1;
		*oneside = 0;
		for (l = 2*k+2; l <= (r+2)+1; ++l)
		{
			s[l] = s[l-1] + 1;
		}
	}
	else if ((N-1) - k < r/2)
	{
		*twoside = 0;
		*oneside = 1;
		for (l = 2*((N-1)-k)+1; l <= (r+2)+1; ++l)
		{
			s[l] = s[l-1] - 1;
		}
	}
	else if ((N-1) - k < (r+2)/2)
	{
		*twoside = 1;
		*oneside = 0;
		for (l = 2*((N-1)-k)+1; l <= (r+2)+1; ++l)
		{
			s[l] = s[l-1] - 1;
		}
	}
	/*
	printf("k = %03i. oneside = %i. twoside = %i. s[] = ",k,*oneside,*twoside);
	for (l = 0; l < (r+2)+1; ++l)
	{
		printf(" %2i,",s[l]);
	}
	printf(" %2i.\n",s[(r+2)+1]);
	*/

	return 0;
}
