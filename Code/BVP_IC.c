#include "BVP_header.h"

int BVP_IC(void *params, double x[], double Psi[])
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
	int i,j,k,l;
	int pn;
	double xk;
	double ICdiff = 1.0;
	double PsiICdiff = 1.0;

	//Generate IC from linearized solution

	for (j = 0; j < n; ++j)
	{
		k = j*p;
		x[j] = (double) j/(n-1);
		xk = x[j];

		//-------------------- Bar G_rr --------------------
		pn = 0;

		//Linear
		//Psi[k +pn] = ICdiff*pow(alpha,2.0)/(beta*kappa*pow(r_H,4.0))*(5460.0*pow(xk,14.0)-38416.0*pow(xk,12.0)+116050.0*pow(xk,10.0)-185504.0*pow(xk,8.0)+158004.0*pow(xk,6.0)-55440.0*pow(xk,4.0)+2383.0*pow(xk,3.0)-154.0*pow(xk,2.0)-2383.0*xk)/1155.0; //Arbitrary % GR GUESS


		//Nonlinear
		//Psi[k +pn] = pow(xk,2.0) + 2.0*ICdiff*pow(alpha,2.0)/(beta*kappa*pow(r_H,4.0))*(5460.0*pow(xk,14.0)-38416.0*pow(xk,12.0)+116050.0*pow(xk,10.0)-185504.0*pow(xk,8.0)+158004.0*pow(xk,6.0)-55440.0*pow(xk,4.0)+2383.0*pow(xk,3.0)-154.0*pow(xk,2.0)-2383.0*xk)/1155.0; //2nd order Lin SOL with x = (1-r_H/4r)/(1+r_H/4r);
		Psi[k +pn] = pow(xk,2.0) + 2.0*ICdiff*pow(alpha,2.0)/(beta*kappa*pow(r_H,4.0))*xk*(xk-1.0)*(xk+1.0)*(5460.0*pow(xk,11.0)-32956.0*pow(xk,9.0)+83094.0*pow(xk,7.0)-102410.0*pow(xk,5.0)+55594.0*pow(xk,3.0)+154.0*xk+2383.0)/1155.0; //Same as above but factored;


		//-------------------- Bar G_thth --------------------
		pn = 1;

		//Linear
		//Psi[k +pn] = -ICdiff*pow(alpha,2.0)/(beta*kappa*pow(r_H,4.0))*32.0/pow(1.0+xk,4.0)*(1610.0*pow(xk,12.0)-11704.0*pow(xk,10.0)+37125.0*pow(xk,8.0)-64064.0*pow(xk,6.0)+62370.0*pow(xk,4.0)-27720.0*pow(xk,2.0)+2383.0*xk)/1155.0; //Arbitrary % GR GUESS


		//Nonlinear
		//Old guess
		//Psi[k +pn] = 16.0/pow(1.0 +xk,4.0) - 2.0*ICdiff*pow(alpha,2.0)/(beta*kappa*pow(r_H,4.0))*32.0/pow(1.0+xk,4.0)*(1610.0*pow(xk,12.0)-11704.0*pow(xk,10.0)+37125.0*pow(xk,8.0)-64064.0*pow(xk,6.0)+62370.0*pow(xk,4.0)-27720.0*pow(xk,2.0)+2383.0*xk)/1155.0; //2nd order Lin SOL with x = (1-r_H/4r)/(1+r_H/4r);
		Psi[k +pn] = 16.0/pow(1.0 +xk,4.0) - ICdiff*pow(alpha,2.0)/(beta*kappa*pow(r_H,4.0))*64.0/pow(1.0+xk,4.0)*xk*(xk-1.0)*(1610.0*pow(xk,10.0)+1610.0*pow(xk,9.0)-10094.0*pow(xk,8.0)-10094.0*pow(xk,7.0)+27031.0*pow(xk,6.0)+27031.0*pow(xk,5.0)-37033.0*pow(xk,4.0)-37033.0*pow(xk,3.0)+25337.0*pow(xk,2.0)+25337.0*xk-2383.0)/1155.0; //same as above but factored

		//-------------------- Psi --------------------
		pn = 2;

		//Linear & Nonlinear
		Psi[k +pn] = PsiICdiff*alpha/beta*((-4.0*pow(xk,6.0) +18.0*pow(xk,4.0) -36.0*pow(xk,2.0) +22.0)/(3.0*pow(r_H,2.0))); //Arbitrary % 1st Order sol

	}


	return 0;
}
